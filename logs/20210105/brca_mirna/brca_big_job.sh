#!/bin/bash
#$ -N brca_many
#$ -cwd


./brca_mirna.py > brca_mirna.log
mkdir -p /data/users/bs16b001/logs/20210105/brca_mirna
cp *.py brca_big_job.sh *.log *.png *.csv /data/users/bs16b001/logs/20210105/brca_mirna
rm *.log *.png *.csv

./brca_meth.py > brca_meth.log
mkdir -p /data/users/bs16b001/logs/20210105/brca_meth
cp *.py *.sh *.log *.csv /data/users/bs16b001/logs/20210105/brca_meth
rm *.log *.csv

./brca_cross.py > brca_cross.log
mkdir -p /data/users/bs16b001/logs/20210105/brca_cross
cp *.py *.sh *.log *.csv /data/users/bs16b001/logs/20210105/brca_cross
rm *.log *.csv

./brca_self.py > brca_self.log
mkdir -p /data/users/bs16b001/logs/20210105/brca_self
cp *.py *.sh *.log *.csv /data/users/bs16b001/logs/20210105/brca_self
rm *.log *.csv

./brca_self_cross.py > brca_self_cross.log
mkdir -p /data/users/bs16b001/logs/20210105/brca_all_int
cp *.py *.sh *.log *.csv /data/users/bs16b001/logs/20210105/brca_al_int
rm *.log *.csv
