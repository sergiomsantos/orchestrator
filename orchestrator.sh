#!/bin/bash

ORCHESTRATOR_HOME="/opt/scientific/gromacs/orchestrator"

NPROC=`cat /proc/cpuinfo | grep processor | wc -l`
WORKERS=$(($NPROC / 2))

$ORCHESTRATOR_HOME/orchestrator "$@" --workers $WORKERS
#python orchestrator.py "$@" --workers $WORKERS

echo "Program finished with exit code $? at: `date`" 


