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
cerebellum = "R/HBTRC/Huntington/cerebellum_vonsattel.csv"
visualCortex = "R/HBTRC/Huntington/visualCortex_vonsattel.csv"
prefrontalCortex = "R/HBTRC/Huntington/prefrontalCortex_vonsattel.csv"
cerebellum_visual = "R/HBTRC/Huntington/cere_visual_vonsattel.csv"
visual_prefrontal = "R/HBTRC/Huntington/visual_prefrontal_vonsattel.csv"
prefrontal_cerebellum = "R/HBTRC/Huntington/prefrontal_cere_vonsattel.csv"
cerebellum_cerebellum = "R/HBTRC/Huntington/cere_cere_vonsattel.csv"
visual_visual = "R/HBTRC/Huntington/visual_visual_vonsattel.csv"
prefrontal_prefrontal = "R/HBTRC/Huntington/prefrontal_prefrontal.csv"
meta_csv = "R/HBTRC/Huntington/labels_vonsattel.csv"
trte_partition_file = "R/HBTRC/Huntington/trte_partition_vonsattel.txt"

# change label from text to integer
#label_dict = {i:i for i in [0,2,3,4]}
label_dict = {0:0, 2:1, 3:2, 4:3}

doSMOTE = True # Boolean

# training parameters
num_epoch = 425
test_interval = 25
lr = 5e-4
weight_decay = 1e-3
dropout = 0.5
adj_parameter = 8 # average number of edge per node in adj matrix

# VERBOSE setting for print results
VERBOSE = 1 #0, only print final result; 1, only testing result; 2, training and testing result
OUTPUT_FILES = False #Boolean to determine whether to output loss and metrics as csv files
MAKE_PLOTS = False #Boolean to determine whether to output loss and metrics as plots in png format
feature_extract = ["shap"] #["lime", "shap"]
num_gcn = 2

RUN_TITLE = "This run includes 9 GCNs - Primary features and all the Interaction features \nPredicting Vonsattel Scores"
RUN_TITLE_SHORT = "vonsattel"

# load preprocessed data from csv
#load_list - list of csv files to laod. The -2 position should be meta_csv and -1 position should be trte_partition_file
load_list = [cerebellum, visualCortex, prefrontalCortex, cerebellum_visual, visual_prefrontal, prefrontal_cerebellum, cerebellum_cerebellum, visual_visual, prefrontal_prefrontal, meta_csv, trte_partition_file]
GCN_names = ["cerebellum", "visualCortex", "prefrontalCortex", "cere_vis", "vis_pre", "pre_cere", "cere_cere", "vis_vis", "pre_pre"]

COMBINER = "FullyConnected"

losses_df, metrics_df, feature_imp, _, _ = train_model(load_list=load_list, label_dict=label_dict, GCN_names=GCN_names, COMBINER=COMBINER, SEED=SEED, num_epoch=num_epoch, test_interval=test_interval, lr=lr, weight_decay=weight_decay,     dropout=dropout, adj_parameter=adj_parameter, VERBOSE=VERBOSE, doSMOTE = doSMOTE, RUN_TITLE=RUN_TITLE, RUN_TITLE_SHORT=RUN_TITLE_SHORT, OUTPUT_FILES=OUTPUT_FILES, MAKE_PLOTS=MAKE_PLOTS, feature_extract=feature_extract, num_gcn=num_gcn)

#losses_df.to_csv("losses.csv")
#metrics_df.to_csv("metrics.csv")
feature_imp["shap"].to_csv("hbtrc_shap_vonsattel_allint.csv", index=False)
