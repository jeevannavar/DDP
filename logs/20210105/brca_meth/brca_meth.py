#!/opt/anaconda3/bin/python3.7

from main import train_model
import pandas as pd

RUN_TITLE = "ANOVA top 1000 features selected greedily with correlation below 0.8 \nThis run includes only one GCN - DNA Methylation"
RUN_TITLE_SHORT = "brca_meth"
# SEED can be "random" or integer, if integer, it will be used as the seed for random, numpy, torch, and cuda
SEED = "random" 

# pre-processed data
mrna_preproc_csv = "/data/users/bs16b001/R/TCGA BRCA/mrna_top1000.csv"
meth_preproc_csv = "/data/users/bs16b001/R/TCGA BRCA/meth_top1000.csv"
mirna_preproc_csv = "/data/users/bs16b001/R/TCGA BRCA/mirna_anova.csv"
meta_csv = "/data/users/bs16b001/R/TCGA BRCA/PAM50_subtype.csv"
trte_partition_file = "/data/users/bs16b001/R/TCGA BRCA/trte_partition.txt"

# change label from text to integer
label_dict = {'Normal':0, 'Basal':1, 'Her2':2, 'LumA':3, 'LumB':4}

# load preprocessed data from csv
#load_list - list of csv files to laod. The -2 position should be meta_csv and -1 position should be trte_partition_file
load_list = [meth_preproc_csv, meta_csv, trte_partition_file]
GCN_names = ["methylation"]#,"methylation","miRNA"]

COMBINER = False # Can take values "VCDN", "FullyConnected", or False.
#Use False when only a sinngle GCN output is to be generated and

doSMOTE = True # Boolean

# MORONET training parameters
num_epoch = 1000
test_interval = 50
lr = 5e-4
weight_decay = 5e-4
dropout = 0.25
adj_parameter = 8 # average number of edge per node in adj matrix

# VERBOSE setting for print results
VERBOSE = 2 #0, only print final result; 1, only testing result; 2, training and testing result
OUTPUT_FILES = False #Boolean to determine whether to output loss and metrics as csv files
MAKE_PLOTS = False #Boolean to determine whether to output loss and metrics as plots in png format
REPEATS = 3 #Integer, how many times to independently train the model

for i in range(1,REPEATS+1):
    print("\n"+"-"*100)
    print("Run {}".format(i))
    print("-"*100)
    losses_df, metrics_df = train_model(load_list=load_list, label_dict=label_dict, GCN_names=GCN_names, COMBINER=COMBINER,
        SEED=SEED, num_epoch=num_epoch, test_interval=test_interval, lr=lr, weight_decay=weight_decay, 
        dropout=dropout, adj_parameter=adj_parameter, VERBOSE=VERBOSE, doSMOTE = doSMOTE,
        RUN_TITLE=RUN_TITLE, RUN_TITLE_SHORT=RUN_TITLE_SHORT,
        OUTPUT_FILES=OUTPUT_FILES, MAKE_PLOTS=MAKE_PLOTS)
    losses_df.to_csv("losses_{}.csv".format(i))
    metrics_df.to_csv("metrics_{}.csv".format(i))
