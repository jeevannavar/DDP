#!/bin/bash
#$ -N python_test
#$ -cwd

chmod 775 main_brca.py

./main_brca.py > output_naivecross.log
