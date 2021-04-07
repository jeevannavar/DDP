#!/bin/bash
#$ -N brca_anova1000
#$ -cwd

chmod 775 main_brca.py

./main_brca.py > anova1000.log
cp *.py *.sh *.log brca_anova1000* *.png *.csv /data/users/bs16b001/logs/20201216/anova/1
rm *.log *.png *.csv

./main_brca.py > anova1000.log
cp *.py *.sh *.log brca_anova1000* *.png *.csv /data/users/bs16b001/logs/20201216/anova/2
rm *.log *.png *.csv

./main_brca.py > anova1000.log
cp *.py *.sh *.log brca_anova1000* *.png *.csv /data/users/bs16b001/logs/20201216/anova/3
rm *.log *.png *.csv
