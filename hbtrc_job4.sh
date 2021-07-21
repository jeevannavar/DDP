#!/bin/bash
#$ -N hbtrc_braak
#$ -cwd

chmod 775 hbtrc_braak.py

./hbtrc_braak.py > hbtrc_braak.log

mkdir -p /data/users/bs16b001/DDP/logs/20210721/hbtrc_braak/
cp *.py hbtrc_job4.sh  hbtrc_braak* hbtrc_shap_braak.csv /data/users/bs16b001/DDP/logs/20210721/hbtrc_braak/
