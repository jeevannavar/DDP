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
cerebellum = "R/HBTRC/Huntington/cerebellum_HD.csv"
visualCortex = "R/HBTRC/Huntington/visualCortex_HD.csv"
prefrontalCortex = "R/HBTRC/Huntington/prefrontalCortex_HD.csv"
meta_csv = "R/HBTRC/Huntington/labels_HD.csv"
trte_partition_file = "R/HBTRC/Huntington/trte_partition_HD.txt"

# change label from text to integer
label_dict = {"control":0, "affected":1}

doSMOTE = True # Boolean

# training parameters
num_epoch = 800
test_interval = 25
lr = 5e-4
weight_decay = 5e-3
dropout = 0.5
adj_parameter = 8 # average number of edge per node in adj matrix

# VERBOSE setting for print results
VERBOSE = 1 #0, only print final result; 1, only testing result; 2, training and testing result
OUTPUT_FILES = False #Boolean to determine whether to output loss and metrics as csv files
MAKE_PLOTS = False #Boolean to determine whether to output loss and metrics as plots in png format
feature_extract = ["shap"] #["lime", "shap"]
num_gcn = 2

RUN_TITLE = "This run includes 3 GCNs - Cerebellum, Primary Visual Cortex, and Prefrontal Cortex"
RUN_TITLE_SHORT = "primary_data"

RUN_TITLE = "This run includes only 3 GCNs - Cerebellum, Visual Cortex, and Prefrontal Cortex \nClassifying between Huntington's Affected and Controls"
RUN_TITLE_SHORT = "HD"

# load preprocessed data from csv
#load_list - list of csv files to laod. The -2 position should be meta_csv and -1 position should be trte_partition_file
load_list = [cerebellum, visualCortex, prefrontalCortex, meta_csv, trte_partition_file]
GCN_names = ["cerebellum", "visualCortex", "prefrontalCortex"]

COMBINER = "VCDN"

losses_df, metrics_df, feature_imp, _, _ = train_model(load_list=load_list, label_dict=label_dict, GCN_names=GCN_names, COMBINER=COMBINER, SEED=SEED, num_epoch=num_epoch, test_interval=test_interval, lr=lr, weight_decay=weight_decay,     dropout=dropout, adj_parameter=adj_parameter, VERBOSE=VERBOSE, doSMOTE = doSMOTE, RUN_TITLE=RUN_TITLE, RUN_TITLE_SHORT=RUN_TITLE_SHORT, OUTPUT_FILES=OUTPUT_FILES, MAKE_PLOTS=MAKE_PLOTS, feature_extract=feature_extract, num_gcn=num_gcn)

#losses_df.to_csv("losses.csv")
#metrics_df.to_csv("metrics.csv")
feature_imp["shap"].to_csv("hbtrc_shap_HD.csv", index_label="features")
