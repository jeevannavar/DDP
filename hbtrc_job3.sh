#!/bin/bash
#$ -N hbtrc_AD
#$ -cwd

chmod 775 hbtrc_AD.py

./hbtrc_AD.py > hbtrc_AD_allint.log

mkdir -p /data/users/bs16b001/DDP/logs/20210726/hbtrc_allint/
cp *.py hbtrc_job3.sh  hbtrc_AD* hbtrc_shap_AD_allint.csv /data/users/bs16b001/DDP/logs/20210726/hbtrc_allint/
