#!/bin/bash
#$ -N brca_meth
#$ -cwd

chmod 775 main_brca.py

./main_brca.py > logs.log
cp *.py *.sh *.log *.png *.csv /data/users/bs16b001/logs/20201217/meth/1
rm *.log *.png *.csv

./main_brca.py > logs.log
cp *.py *.sh *.log *.png *.csv /data/users/bs16b001/logs/20201217/meth/2
rm *.log *.png *.csv

./main_brca.py > logs.log
cp *.py *.sh *.log brca_meth* *.png *.csv /data/users/bs16b001/logs/20201217/meth/3
rm *.log *.png *.csv brca_meth*
