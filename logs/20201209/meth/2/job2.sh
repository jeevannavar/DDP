#!/bin/bash
#$ -N solo_meth
#$ -cwd

chmod 775 main_brca.py

./main_brca.py > output_meth.log
cp *.log *.csv *.png *.py *.sh solo* /data/users/bs16b001/logs/20201209/meth/2
rm *.log *.png *.csv solo*
