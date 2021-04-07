#!/bin/bash
#$ -N six_gcns
#$ -cwd

chmod 775 main_brca.py

./main_brca.py > output_3GCNs_3naiveGCNs.log
