#!/bin/bash
#$ -N brca_many
#$ -cwd


./brca_self_cross.py > brca_self_cross.log
mkdir -p /data/users/bs16b001/logs/20210105/brca_all_int
cp *.py *.sh *.log *.csv /data/users/bs16b001/logs/20210105/brca_all_int
rm *.log *.csv
