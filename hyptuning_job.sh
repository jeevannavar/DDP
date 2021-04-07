#!/bin/bash
#$ -N brca_hyptuning
#$ -cwd

chmod 775 hyp_tuner.py

./hyp_tuner.py > brca_hyptuning.log
mkdir -p /data/users/bs16b001/logs/20201230/hyp_tuner/brca_greedyanova
cp *.py hbtrc_job.sh *.log hbtrc_elasticnet* *.png *.csv /data/users/bs16b001/logs/20201230/hyp_tuner/brca_greedyanova
rm *.log *.png *.csv
