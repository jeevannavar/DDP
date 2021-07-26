#!/bin/bash
#$ -N hbtrc_HD
#$ -cwd

chmod 775 hbtrc_HD.py

./hbtrc_HD.py > hbtrc_HD_allint.log

mkdir -p /data/users/bs16b001/DDP/logs/20210726/hbtrc_HD_allint/
cp *.py hbtrc_job1.sh  hbtrc_HD* hbtrc_shap_HD_allint.csv /data/users/bs16b001/DDP/logs/20210726/hbtrc_HD_allint/
