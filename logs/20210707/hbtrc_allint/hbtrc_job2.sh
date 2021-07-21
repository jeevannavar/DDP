#!/bin/bash
#$ -N hbtrc_allint
#$ -cwd

chmod 775 hbtrc_allint.py

./hbtrc_allint.py > hbtrc_allint.log

mkdir -p /data/users/bs16b001/DDP/logs/20210707/hbtrc_allint/
cp *.py hbtrc_job2.sh  hbtrc_allint* hbtrc_lime_allint.csv /data/users/bs16b001/DDP/logs/20210707/hbtrc_allint/