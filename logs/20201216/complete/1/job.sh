#!/bin/bash
#$ -N brca_complete
#$ -cwd

chmod 775 main_brca.py

./main_brca.py > complete.log
cp *.py *.sh *.log *.png *.csv /data/users/bs16b001/logs/20201216/complete/1
rm *.log *.png *.csv

./main_brca.py > complete.log
cp *.py *.sh *.log *.png *.csv /data/users/bs16b001/logs/20201216/complete/2
rm *.log *.png *.csv

./main_brca.py > complete.log
cp *.py *.sh *.log brca_complete* *.png *.csv /data/users/bs16b001/logs/20201216/complete/3
rm *.log *.png *.csv brca_complete*
