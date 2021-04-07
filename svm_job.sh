#!/bin/bash
#$ -N svm_brca
#$ -cwd

/opt/anaconda3/bin/jupyter nbconvert --to notebook --execute svm.ipynb --output svm_run.ipynb --ExecutePreprocessor.timeout=None > svm.log
