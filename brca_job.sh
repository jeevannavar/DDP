#!/bin/bash
#$ -N brca_allint_extract
#$ -cwd

chmod 775 brca.py

./brca.py > feature_extraction_all_int.log
