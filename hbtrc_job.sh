#!/bin/bash
#$ -N hbtrc_HD
#$ -cwd

chmod 775 hbtrc_HD.py

./hbtrc_HD.py > hbtrc_HD.log

mkdir -p /data/users/bs16b001/DDP/logs/20210721/hbtrc_HD/
cp *.py hbtrc_job.sh  hbtrc_HD* hbtrc_shap_HD.csv /data/users/bs16b001/DDP/logs/20210721/hbtrc_HD/
