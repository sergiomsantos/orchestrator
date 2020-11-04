'''
This is my package description

(c) Sergio M. Santos, Univ Aveiro - PT, 2018

'''

__version__ = '0.2.1'


from pyanitools import anidataloader

from queue import Queue
import subprocess
import threading
import datetime
import logging
import os

import h5py
import json

TEMPLATE = '''%nproc=12
%mem={mem}
#n B3LYP/6-31G(d,p) nmr

{{
 "smiles": "{smiles}",
 "path": "{path}",
 "index": {index}
}}

0 1
{conf}

'''


class ProducerThread(threading.Thread):

    def __init__(self, queue, src_file, resume=0, start=0, stop=None, step=1, batch_size=1, sort=False, *args, **kwargs):
        super(ProducerThread, self).__init__(*args, **kwargs)
        self.slice = slice(start, stop, step)
        self.batch_size = batch_size
        self.queue = queue
        if not os.path.exists(src_file):
            error_msg = f'File {src_file} not found'
            self.logger.error(error_msg)
            raise IOError(error_msg)
        self.src_file = src_file
        self.sort = sort
        self.resume = resume
        self.logger = kwargs.get('logger', logging.getLogger(__name__))
        self.redis_conn = kwargs.get('redis_conn', None)
        if self.redis_conn is None:
            self.logger.info('No redis connection was provided')

        self._jobs = []
        self._stopevent = threading.Event()
        self._count = 0
        self._total = 0
        self._done = 0

    def stop(self):
        self._stopevent.set()
        self._redis_hset('status', 'killed-pending')

    def update(self, n=1, refresh=False):
        self._count += n
        if refresh or (self._count % 50 == 0):
            prefix = '' if self.is_deployable() else '(resuming) '
            if self._done == 0:
                self.logger.info(f'{prefix}deployed {self._count}/{self._total} jobs')
            else: 
                now = datetime.datetime.now()
                ellapsed = now - self._start_time
                # tpj = datetime.timedelta(seconds=time_per_job)
                time_per_job = ellapsed.total_seconds() / self._done
                eta = self._start_time + datetime.timedelta(seconds=(time_per_job * self._total))
                self.logger.info(f'{prefix}deployed {self._count}/{self._total} jobs (ETA: {eta.strftime("%c")})')
        self._redis_hincrby('ndone', n)
    
    # def count(self):
    #     return self._count

    def _redis_hincrby(self, key, val):
        if self.redis_conn is not None:
            self.redis_conn.hincrby(self.session, key, val)
    
    def _redis_hset(self, key, val):
        if self.redis_conn is not None:
            self.redis_conn.hset(self.session, key, val)

    def visit(self, name, item):
        if isinstance(item, h5py.Dataset) and name.endswith('/coordinates'):
             item_len = item.len()
             n_calc = len(range(*self.slice.indices(item_len)))
             self._total += n_calc
             self.logger.info(f'  > {name} = {n_calc}/{item_len}')
             self._redis_hincrby('ntotal', n_calc)

    def run(self):
        self._start_time = datetime.datetime.now()
        self.logger.info('starting ...')
        self.logger.info(f'Processing file {self.src_file}')
        self.logger.info(f'  > batch-size = {self.batch_size}')
        self.logger.info(f'  >     resume = {self.resume}')
        self.logger.info(f'  >      start = {self.slice.start}')
        self.logger.info(f'  >       stop = {self.slice.stop}')
        self.logger.info(f'  >       step = {self.slice.step}')
        self.logger.info(f'  >       sort = {self.sort}')
        
        self._redis_hset('status', 'running')
        self._redis_hset('details', json.dumps(dict(
            file=self.src_file,
            start=self.slice.start,
            stop=self.slice.stop,
            step=self.slice.step,
            batch_size=self.batch_size,
            sort=self.sort
        )))

        dataloader = anidataloader(self.src_file)
        self.logger.info(f'Groups and dataset lengths:')
        dataloader.store.visititems(self.visit)
        self.logger.info(f'Done listing existing groups.')
        self.update(0, refresh=True)

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
            
            if self.is_deployable():
                self.logger.info(f'Processing group {path}')
                self.logger.info(f'  > Will produce {len(indices)} items')

            for index,frame in zip(indices,frames):
                if self._stopevent.isSet():
                    break
                job = (index, species, smiles, frame, path)
                self.append_job(job)

            if self._stopevent.isSet():
                self.logger.warning('received early termination signal')
                break

            if self.is_deployable():
                self.logger.info(f'flushing jobs for group {path}')
            self.flush_job()
        self.update(0, refresh=True)

        self.kill_workers()
        dataloader.cleanup()

        st = 'killed' if self._stopevent.isSet() else 'done'
        self._redis_hset('status', st)
        self.logger.info('completed!')

    def kill_workers(self):
        self.logger.info(f'deploying {self.queue.maxsize} poison pills')
        if self.queue.qsize() > 0:
            self.logger.warning('  > queue is not empty')
        for n in range(self.queue.maxsize):
            self.queue.put(None)

    def append_job(self, job):
        self.update(1)
        if self.is_deployable():
            self._done += 1
        self._jobs.append(job)
        if len(self._jobs) == self.batch_size:
            self.flush_job()

    def flush_job(self):
        if self._jobs:
            if self.is_deployable():
                self.queue.put(self._jobs)
            self._jobs = []
    
    def is_deployable(self):
        return (self._count >= self.resume)


class ConsumerThread(threading.Thread):

    def __init__(self, queue, cmd, session, mem,  *args, **kwargs):
        super(ConsumerThread, self).__init__(*args, **kwargs)
        self.queue = queue
        self.session = session
        self.logger = kwargs.get('logger', logging.getLogger(__name__))
        self.cmd = cmd
        self.mem = mem
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
        return TEMPLATE.format(mem=self.mem, smiles=smiles, index=index, path=path, conf=conf)


def main(args):
    import signal
    import time

    log = logging.getLogger()

    try:
        import sshtunnel
        import redis
    except ImportError as e:
        log.error(str(e))

    # try:
    #     tunnel = sshtunnel.SSHTunnelForwarder(
    #         args.ssh_host
    #         ssh_username=args.ssh_username,
    #         ssh_password=args.ssh_password,
    #         remote_bind_address=('localhost', 6379)
    #     )
    #     tunnel.start()
    #     redis_conn = redis.StrictRedis(host='localhost', port=tunnel.local_bind_port)
    # except Exception as e:
    #     logger.error(str(e))
    #     redis_conn = None
    #     tunnel = None
    # else:
    #     logger.info(f'SSHTunnelForwarder - connected to localhost:{tunnel.local_bind_port}')


    queue = Queue(args.workers)
    workers = []

    for n in range(args.workers):
        c = ConsumerThread(queue, args.cmd, args.session, args.mem, name=f'worker-{n}')
        workers.append(c)

    producer = ProducerThread(
        queue, args.file, resume=args.resume,
        start=args.start, stop=args.stop, step=args.step,
        batch_size=args.batch_size,
        sort=args.sort,
        name=f'producer')
    workers.append(producer)

    for worker in workers:
        worker.start()
        time.sleep(1)

    def sigterm_handler(*args, **kwargs):
       log.warning('Received a SIGTERM/SIGINT signal')
       log.warning('  > requesting early termination')
       producer.stop()

    signal.signal(signal.SIGTERM, sigterm_handler)
    signal.signal(signal.SIGINT, sigterm_handler)

    for worker in workers:
        worker.join()
    
    # if tunnel is not None:
    #     logger.info('SSHTunnelForwarder - closing tunnel')
    #     tunnel.close()


if __name__ == '__main__':
    import argparse
    import uuid
    import sys
    import os

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', action='version', version=__version__)
    
    parser.add_argument('file', help='h5df source file')
    parser.add_argument('-c','--cmd', required=True, type=str, help='command to execute: should take a G09.com file as input')
    parser.add_argument('--workers', help='number of workers', default=1, type=int)
    parser.add_argument('--start', help='slice start', default=0, type=int)
    parser.add_argument('--stop', help='slice stop', default=None, type=int)
    parser.add_argument('--step', help='slice step', default=1, type=int)
    parser.add_argument('--session', help='session identifier', default=None, type=str)
    parser.add_argument('--batch-size', help='configurations per file', default=10, type=int)
    parser.add_argument('--sort', help='sort configurations prior to selection, according to their "energy" attribute', action='store_true')
    parser.add_argument('--resume', help='resume from job number <JOB>', type=int, default=0, metavar='JOB')
    parser.add_argument('--mem', help='memory per job', default='2GB', type=str)

    
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
    start_time = datetime.datetime.now()
    log.info(f'Entering orchestrator on {start_time.strftime("%c")}')
    log.info(f'session name: {args.session}')

    main(args)

    end_time = datetime.datetime.now()
    log.info(f'All done {end_time.strftime("%c")}')
    log.info(f'  > Ellapsed time = {str(end_time-start_time)}')

