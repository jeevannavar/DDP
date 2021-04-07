#!/bin/bash
#$ -N hbtrc_elasticnet
#$ -cwd

chmod 775 hbtrc.py

./hbtrc.py > hbtrc_elasticnet.log
mkdir -p /data/users/bs16b001/logs/20210101/hbtrc/elasticnet
cp *.py hbtrc_job.sh *.log hbtrc_elasticnet* *.png *.csv /data/users/bs16b001/logs/20210101/hbtrc/elasticnet
rm *.log *.png *.csv hbtrc_elasticnet*
