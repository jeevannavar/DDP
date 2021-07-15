#!/bin/bash
#$ -N citeseq_imp400
#$ -cwd

/opt/anaconda3/bin/python3 citeseq_imputation.py > citeseq_imp400.log
