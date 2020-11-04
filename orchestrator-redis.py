from pyanitools import anidataloader

from queue import Queue
import subprocess
import threading
import logging
import os

import h5py

import sshtunel
import redis

#server = sshtunnel.SSHTunnelForwarder('193.137.169.65', ssh_username='ssantos', ssh_password='passwd', remote_bind_address=('localhost', 6379))
#server.start()
#r = redis.StrictRedis(host='localhost', port=6379)
#r = redis.StrictRedis(host='localhost', port=server.local_bind_port)
#r.get('foo').decode('utf8')
#server.stop()


TEMPLATE = '''%nproc=2
%mem=2GB
#n B3LYP/6-31G(d,p) nmr

{{
 "smiles": "{smiles}",
 "path": "{path}",
 "index": {index}
}}

0 1
{conf}

'''

def write(species, coordinates, smiles, path, index):
    lines = []
    for s,(x,y,z) in zip(species, coordinates):
        lines.append('%-2s  %12.6f %12.6f %12.6f'%(s,x,y,z))
    conf = '\n'.join(lines)
    return TEMPLATE.format(smiles=smiles, index=index, path=path, conf=conf)

class Counter(object):
    def __init__(self, logger, slicer, logfreq=100, session=None, host=None):
        self.logfreq = logfreq
        self.logger = logger
        self.slice = slicer
        self._count = 0
        self._total = 0
        if session and host:
            h,p = host
            self.redis = redis.StrictRedis(host=h, port=p)
        else:
            self.redis = None
        self.session = session

    def visit(self, name, item):
        if isinstance(item, h5py.Dataset) and name.endswith('/coordinates'):
             item_len = item.len()
             n_calc = len(range(*self.slice.indices(item_len)))
             self._total += n_calc
             self.logger.info(f'  > {name} = [{n_calc}/{item_len}]')
             self._redis_hincrby('ntotal', n_calc)
             
    def update(self, n=1, refresh=False):
        self._count += n
        if refresh or (self._count%self.logfreq == 0):
            self.logger.info(f'deployed [{self._count}/{self._total}] jobs')
        self._redis_hincrby('ndone', n)
        
    def count(self):
        return self._count
    
    def _redis_hincrby(self, key, val):
        if self.redis is not None:
            self.redis.hincrby(self.session, key, val)
    
class ProducerThread(threading.Thread):
    
    def __init__(self, queue, src_file, redis_host, session, start=0, stop=None, step=1, batch_size=1, sort=False, *args, **kwargs):
        super(ProducerThread, self).__init__(*args, **kwargs)
        self.slice = slice(start, stop, step)
        self.batch_size = batch_size
        self.src_file = src_file
        self.queue = queue
        self.jobs = []
        self.logger = kwargs.get('logger', logging.getLogger(__name__))
        if not os.path.exists(src_file):
            error_msg = f'File {src_file} not found'
            self.logger.error(error_msg)
            raise IOError(error_msg)
        self.sort = sort
        self._stopevent = threading.Event()
        self.redis_host = redis_host
        self.session = session

    def stop(self):
        self._stopevent.set()

    def run(self):
        
        self.logger.info('starting ...')
        self.logger.info(f'Processing file {self.src_file}')
        self.logger.info(f'  >      start = {self.slice.start}')
        self.logger.info(f'  >       stop = {self.slice.stop}')
        self.logger.info(f'  >       step = {self.slice.step}')
        self.logger.info(f'  > batch-size = {self.batch_size}')
        
        #def visit(name, item):
        #    if isinstance(item, h5py.Dataset) and name.endswith('/coordinates'):
        #       self.logger.info(f'  > {name} = {item.len()}')
        
        #dataloader = anidataloader(self.src_file)
        #self.logger.info(f'Groups and dataset lengths:')
        #dataloader.store.visititems(visit)
        #self.logger.info(f'Done listing existing groups.')
        
        dataloader = anidataloader(self.src_file)
        counter = Counter(self.logger, self.slice, logfreq=50, session=self.session, host=self.redis_host)
        self.logger.info(f'Groups and dataset lengths:')
        dataloader.store.visititems(counter.visit)
        self.logger.info(f'Done listing existing groups.')
        counter.update(0, refresh=True)
        #job_counter = 0
        
        for dataset in dataloader:
            
            frames = dataset['coordinates']
            if self.sort:
                indices = dataset['energies'].argsort()
                indices = indices[self.slice]
            else:
                indices = range(*self.slice.indices(len(frames)))
            frames  = frames[self.slice]
            smiles  = ''.join(dataset['smiles'])
            species = dataset['species']
            path = dataset['path']
            self.logger.info(f'Processing group {path}')
            self.logger.info(f'  > Will produce {len(indices)} items')
            
            for index,frame in zip(indices,frames):
                if self._stopevent.isSet():
                    break
                job = (index, species, smiles, frame, path)
                self.append_job(job)
                #job_counter += 1
                #self.queue.put(job)

                #if job_counter%100 == 0:
                #    self.logger.info(f'deployed {job_counter} jobs so far')
                counter.update()
            if self._stopevent.isSet():
                self.logger.warning('received early termination signal')
                break
            
            self.logger.info(f'flushing jobs for group {path}')
            self.flush_job()
        counter.update(0, refresh=True)
        #self.logger.info(f'scheduled a total of {counter.count()} jobs')
        self.kill_workers()
        dataloader.cleanup()

        self.logger.info('completed!')

    def kill_workers(self):
        
        self.logger.info(f'deploying {self.queue.maxsize} poison pills')
        if self.queue.qsize() > 0:
            self.logger.warning('  > queue is not empty')
        for n in range(self.queue.maxsize):
            self.queue.put(None)
 
    def append_job(self, job):
        self.jobs.append(job)
        if len(self.jobs) == self.batch_size:
            self.flush_job()
    
    def flush_job(self):
        if self.jobs:
            self.queue.put(self.jobs)
            self.jobs = []

class ConsumerThread(threading.Thread):

    def __init__(self, queue, cmd, session, *args, **kwargs):
        super(ConsumerThread, self).__init__(*args, **kwargs)
        self.queue = queue
        self.session = session
        self.logger = kwargs.get('logger', logging.getLogger(__name__))
        self.cmd = cmd
        self._pass_counter = 0
        self._fail_counter = 0

    def run(self):
        self.logger.info('starting ...')

        dirname = f'{self.session}/{self.name}'
        if not os.path.exists(dirname):
            self.logger.info(f'  > creating directory {dirname}')
            os.makedirs(dirname)
        
        while True:
            job = self.queue.get()

            if job is None:
                self.queue.task_done()
                self.logger.info('received poison pill')
                break

            self.exec_job(job)
            self.queue.task_done()
        
        #ntot = self._pass_counter + self._fail_counter
        #self.logger.info(f'successfully executed {self._pass_counter}/{ntot} jobs')
        self.logger.info(f'stats:')
        self.logger.info(f'  >  fails = {self._fail_counter}')
        self.logger.info(f'  > passes = {self._pass_counter}')
        self.logger.info( 'completed!')
    

    def exec_job(self, job):
                
        index, species, smiles, frame, path = job[0]
        
        prefix = f'{self.session}/{self.name}/{os.path.basename(path)}__{index}'
         
        # write gaussian input file
        with open(prefix + '.com', 'wt') as fout:            
            ginput = '--link1--\n'.join(self.as_gin(*j) for j in job)
            fout.write(ginput)
        
        # execute command
        completed = subprocess.run(f'{self.cmd} {prefix}', shell=True)
        
        self.logger.debug(f'{prefix} - {completed.returncode}')
        if completed.returncode != 0:
            self.logger.error(f'Unable to successfully run job > {self.cmd} {prefix}')
            self._fail_counter += 1
        else:
            self._pass_counter += 1

    def as_gin(self, index, species, smiles, xyz, path):
        lines = []
        for s,(x,y,z) in zip(species, xyz):
            lines.append('%-2s  %12.6f %12.6f %12.6f'%(s,x,y,z))
        conf = '\n'.join(lines)
        return TEMPLATE.format(smiles=smiles, index=index, path=path, conf=conf)


def main(args):
    import signal
    import time
    
    queue = Queue(args.workers)
    workers = []
     
    try:
        tunnel = sshtunnel.SSHTunnelForwarder(
            '193.137.169.65',
            ssh_username='ssantos',
            ssh_password='passwd',
            remote_bind_address=('localhost', 6379)
        )
    except Exception as e:
        print(e)
        tunnel = None
        redis_host = None
    else:
        redis_host = ('localhost', tunnel.local_bind_port)

    for n in range(args.workers):
        c = ConsumerThread(queue, args.cmd, args.session, name=f'worker-{n}')
        workers.append(c)
    
    producer = ProducerThread(
        queue, args.file, redis_host, args.session,
        start=args.start, stop=args.stop, step=args.step,
        batch_size=args.batch_size,
        sort=args.sort,
        name=f'producer')
    workers.append(producer)
    
    for worker in workers: 
        worker.start()
        time.sleep(1)
    def sigterm_handler(*args, **kwargs):
       log = logging.getLogger()
       log.warning('Received a SIGTERM/SIGINT signal')
       log.warning('  > requesting early termination')
       producer.stop()
    
    signal.signal(signal.SIGTERM, sigterm_handler)
    signal.signal(signal.SIGINT, sigterm_handler)

    for worker in workers: 
        worker.join()
    
    if tunnel is not None:
        tunnel.close()

if __name__ == '__main__':
    from datetime import datetime
    import argparse
    import uuid
    import sys
    import os

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('file', help='h5df source file')
    parser.add_argument('-c','--cmd', required=True, type=str, help='command to execute: should take a G09.com file as input')
    parser.add_argument('--workers', help='number of workers', default=1, type=int)
    parser.add_argument('--start', help='slice start', default=0, type=int)
    parser.add_argument('--stop', help='slice stop', default=None, type=int)
    parser.add_argument('--step', help='slice step', default=1, type=int)
    parser.add_argument('--session', help='session identifier', default=None, type=str)
    parser.add_argument('--batch-size', help='configurations per file', default=10, type=int)
    parser.add_argument('--sort', help='sort configurations prior to selection, according to their "energy" attribute', action='store_true')

    args = parser.parse_args()
    #print(args)
    #exit(1)
    
    if args.session is None:
        args.session = uuid.uuid4().hex[-12:]
    
    # set logger options
    # ---------------------------------
    log = logging.getLogger()
    log.setLevel(logging.INFO)
    handler = logging.FileHandler(f'{args.session}.log')
    # handler = logging.StreamHandler(sys.stdout)
    # formatter = logging.Formatter('[%(levelname)-5s%(threadName)20s:%(lineno)-4i] %(message)s')
    formatter = logging.Formatter('[%(levelname)-5s %(threadName)s] %(message)s')
    handler.setFormatter(formatter)
    log.addHandler(handler)
    
    # MAIN 
    # ---------------------------------
    start_time = datetime.now() 
    log.info(f'Entering orchestrator on {start_time.strftime("%c")}')
    log.info(f'session name: {args.session}')
    
    main(args)
    
    end_time = datetime.now()
    log.info(f'All done {end_time.strftime("%c")}')
    log.info(f'  > Ellapsed time = {str(end_time-start_time)}')
    
    
