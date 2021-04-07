#!/bin/bash
#$ -N brca_mirna
#$ -cwd

chmod 775 main_brca.py

./main_brca.py > logs.log
cp *.py *.sh *.log *.png *.csv /data/users/bs16b001/logs/20201217/mirna/1
rm *.log *.png *.csv

./main_brca.py > logs.log
cp *.py *.sh *.log *.png *.csv /data/users/bs16b001/logs/20201217/mirna/2
rm *.log *.png *.csv

./main_brca.py > logs.log
cp *.py *.sh *.log brca_mirna* *.png *.csv /data/users/bs16b001/logs/20201217/mirna/3
rm *.log *.png *.csv brca_mirna*
