#!/bin/bash
#$ -N brca_complete
#$ -cwd

chmod 775 main_brca.py

./main_brca.py > output_complete.log
