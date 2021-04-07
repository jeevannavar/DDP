#!/bin/bash
#$ -N brca_hyptuning
#$ -cwd

chmod 775 hyp_tuner.py

./hyp_tuner.py > brca_hyptuning.log
mkdir -p /data/users/bs16b001/logs/20210101/hyp_tuner/brca_elasticnet
cp *.py hyptuning_job.sh *.log brca_hyptuning* *.png *.csv /data/users/bs16b001/logs/20210101/hyp_tuner/brca_elasticnet
rm *.log *.png *.csv brca_hyptuning*
