#!/bin/bash
#$ -N complete_data
#$ -cwd

chmod 775 main_brca.py

./main_brca.py > output_complete.log
