#!/bin/bash
#$ -N hbtrc
#$ -cwd

chmod 775 hbtrc.py

/opt/anaconda3/bin/jupyter nbconvert --to notebook --execute HBTRC_Multi_Tissue_Analysis.ipynb --output Multi-Tissue-Analysis.ipynb --ExecutePreprocessor.timeout=None

mkdir -p /data/users/bs16b001/DDP/logs/20210622/hbtrc/
cp *.py hbtrc_job.sh *Analysis.ipynb  hbtrc* *.csv /data/users/bs16b001/DDP/logs/20210622/hbtrc/
rm *.csv
