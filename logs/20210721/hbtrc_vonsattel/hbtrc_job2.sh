#!/bin/bash
#$ -N hbtrc_vonsattel
#$ -cwd

chmod 775 hbtrc_vonsattel.py

./hbtrc_vonsattel.py > hbtrc_vonsattel.log

mkdir -p /data/users/bs16b001/DDP/logs/20210721/hbtrc_vonsattel/
cp *.py hbtrc_job2.sh  hbtrc_vonsattel* hbtrc_shap_vonsattel.csv /data/users/bs16b001/DDP/logs/20210721/hbtrc_vonsattel/
