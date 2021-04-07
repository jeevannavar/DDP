#!/bin/bash
#$ -N solo_mirna
#$ -cwd

chmod 775 main_brca.py

./main_brca.py > output_mirna.log
