#!/bin/bash
#$ -N ind_validation
#$ -cwd

/opt/anaconda3/bin/jupyter nbconvert --to notebook --execute Independent_Validation.ipynb --output output.ipynb --ExecutePreprocessor.timeout=None