#!/bin/bash
#$ -N hbtrc_allfeat
#$ -cwd

chmod 775 hbtrc_allfeat.py

./hbtrc_allfeat.py > hbtrc_allfeat.log

mkdir -p /data/users/bs16b001/DDP/logs/20210707/hbtrc_allfeat/
cp *.py hbtrc_job.sh  hbtrc_allfeat* hbtrc_lime_primary.csv /data/users/bs16b001/DDP/logs/20210707/hbtrc_allfeat/
