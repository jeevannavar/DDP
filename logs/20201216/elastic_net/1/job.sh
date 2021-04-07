#!/bin/bash
#$ -N brca_elasticnet
#$ -cwd

chmod 775 main_brca.py

./main_brca.py > elastic.log
cp *.py *.sh *.log *.png *.csv /data/users/bs16b001/logs/20201216/elastic_net/1
rm *.log *.png *.csv

./main_brca.py > elastic.log
cp *.py *.sh *.log *.png *.csv /data/users/bs16b001/logs/20201216/elastic_net/2
rm *.log *.png *.csv

./main_brca.py > elastic.log
cp *.py *.sh *.log brca_elastic* *.png *.csv /data/users/bs16b001/logs/20201216/elastic_net/3
rm *.log *.png *.csv brca_elastic*
