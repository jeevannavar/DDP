#!/bin/bash
#$ -N hbtrc_vonsattel
#$ -cwd

chmod 775 hbtrc_vonsattel.py

./hbtrc_vonsattel.py > hbtrc_vonsattel_allint.log

mkdir -p /data/users/bs16b001/DDP/logs/20210726/hbtrc_vonsattel_allint/
cp *.py hbtrc_job2.sh  hbtrc_vonsattel* hbtrc_shap_vonsattel_allint.csv /data/users/bs16b001/DDP/logs/20210726/hbtrc_vonsattel_allint/
