#!/opt/anaconda3/bin/python3.7
#!/usr/bin/python3

from main import train_model
import pandas as pd
import numpy as np
import random

import cufflinks as cf
import pathlib
#%matplotlib inline

import shap
shap.initjs()
import lime
import re

# SEED can be "random" or integer, if integer, it will be used as the seed for random, numpy, torch, and cuda
SEED = 42

# pre-processed data
cerebellum = "R/HBTRC/cerebellum_top.csv"
visualCortex = "R/HBTRC/visualCortex_top.csv"
prefrontalCortex = "R/HBTRC/prefrontalCortex_top.csv"
cerebellum_visual = "R/HBTRC/cere_visual.csv"
visual_prefrontal = "R/HBTRC/visual_prefrontal.csv"
prefrontal_cerebellum = "R/HBTRC/prefrontal_cere.csv"
cerebellum_cerebellum = "R/HBTRC/cere_cere.csv"
visual_visual = "R/HBTRC/visual_visual.csv"
prefrontal_prefrontal = "R/HBTRC/prefrontal_prefrontal.csv"
meta_csv = "R/HBTRC/disease_class.csv"
trte_partition_file = "R/HBTRC/trte_partition.txt"

# change label from text to integer
label_dict = {"Normal   ":0, "Alzheimer's":1, "Huntington's":2}

doSMOTE = True # Boolean

# training parameters
num_epoch = 900
test_interval = 50
lr = 5e-4
weight_decay = 1e-3
dropout = 0.25
adj_parameter = 8 # average number of edge per node in adj matrix

# VERBOSE setting for print results
VERBOSE = 1 #0, only print final result; 1, only testing result; 2, training and testing result
OUTPUT_FILES = False #Boolean to determine whether to output loss and metrics as csv files
MAKE_PLOTS = False #Boolean to determine whether to output loss and metrics as plots in png format
feature_extract = [] #["lime", "shap"]
num_gcn = 2

RUN_TITLE = "This run includes 3 GCNs - Cerebellum, Primary Visual Cortex, and Prefrontal Cortex"
RUN_TITLE_SHORT = "primary_data"

# load preprocessed data from csv
#load_list - list of csv files to laod. The -2 position should be meta_csv and -1 position should be trte_partition_file
load_list = [cerebellum, visualCortex, prefrontalCortex, meta_csv, trte_partition_file]
GCN_names = ["cerebellum", "visualCortex", "prefrontalCortex"]

COMBINER = "VCDN"
feature_extract = ["lime"]

losses_df, metrics_df, feature_imp, _, _ = train_model(load_list=load_list, label_dict=label_dict, GCN_names=GCN_names, COMBINER=COMBINER, SEED=SEED, num_epoch=num_epoch, test_interval=test_interval, lr=lr, weight_decay=weight_decay,     dropout=dropout, adj_parameter=adj_parameter, VERBOSE=VERBOSE, doSMOTE = doSMOTE, RUN_TITLE=RUN_TITLE, RUN_TITLE_SHORT=RUN_TITLE_SHORT, OUTPUT_FILES=OUTPUT_FILES, MAKE_PLOTS=MAKE_PLOTS, feature_extract=feature_extract, num_gcn=num_gcn)

#losses_df.to_csv("losses.csv")
#metrics_df.to_csv("metrics.csv")
feature_imp["lime"].to_csv("hbtrc_lime_primary.csv", index_label="features")
