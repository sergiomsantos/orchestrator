#!/usr/bin/env bash

NOW=$(date +"%s")
G09_SCRATCH=scratch_$NOW$RANDOM

#echo "Starting at `date`"
#echo "Running on hosts: $SLURM_NODELIST"
#echo "Running on $SLURM_NNODES nodes."
#echo "Running on $SLURM_NPROCS processors."
#echo "Current working directory is `pwd`"

mkdir -p $G09_SCRATCH
export g09root=/opt/scientific/g09
export GAUSS_SCRDIR=$G09_SCRATCH
source $g09root/g09/bsd/g09.profile
#export OMP_NUM_THREADS=2

g09 < $1.com > $1.log

gretval=$?

if [ $gretval -eq 0 ]; then
    gzip $1.log
fi

rm -rf $G09_SCRATCH

exit $gretval

