#!/bin/bash
#$ -N hbtrc_AD
#$ -cwd

chmod 775 hbtrc_AD.py

./hbtrc_AD.py > hbtrc_AD.log

mkdir -p /data/users/bs16b001/DDP/logs/20210721/hbtrc_AD/
cp *.py hbtrc_job3.sh  hbtrc_AD* hbtrc_shap_AD.csv /data/users/bs16b001/DDP/logs/20210721/hbtrc_AD/
