#!/bin/bash
#$ -N brca_mrna
#$ -cwd

chmod 775 brca.py

./brca.py > brca_mrna.log
DIR=/data/users/bs16b001/logs/20210104/brca_mrna
mkdir -p DIR
cp *.py brca_job.sh *.log brca_mrna* *.png *.csv DIR
rm *.log *.png *.csv brca_mrna*
