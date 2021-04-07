#!/bin/bash
#$ -N solo_mrna
#$ -cwd

chmod 775 main_brca.py

./main_brca.py > output_mrna.log
cp *.log *.csv *.png *.py *.sh solo* /data/users/bs16b001/logs/20201209/mrna/3
rm *.log *.png *.csv solo*
