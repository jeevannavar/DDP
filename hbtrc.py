#!/opt/anaconda3/bin/python3.7

from main import train_model
import pandas as pd

RUN_TITLE = "HBTRC Huntington's data. Feature selected using elastic net. \nThis run includes 3 tissues - cerebellum, primary visual cortex, and prefrontal cortex"
RUN_TITLE_SHORT = "hbtrc_elasticnet"
# SEED can be "random" or integer, if integer, it will be used as the seed for random, numpy, torch, and cuda
SEED = "random" 

# pre-processed data
cerebellum_csv = "HBTRC/cerebellum_elasticnet.csv"#"/data/users/bs16b001/R/HBTRC/cerebellum_elasticnet.csv"
visualCortex_csv = "HBTRC/visualCortex_elasticnet.csv"#"/data/users/bs16b001/R/HBTRC/visualCortex_elasticnet.csv"
prefrontalCortex_csv = "HBTRC/prefrontalCortex_elasticnet.csv"#"/data/users/bs16b001/R/HBTRC/prefrontalCortex_elasticnet.csv"
meta_csv = "HBTRC/disease_class.csv"#"/data/users/bs16b001/R/HBTRC/disease_class.csv"
trte_partition_file = "HBTRC/trte_partition.txt"#"/data/users/bs16b001/R/HBTRC/trte_partition.txt"

# change label from text to integer
label_dict = {"Normal   ":0, "Alzheimer's":1, "Huntington's":2}

# load preprocessed data from csv
#load_list - list of csv files to laod. The -2 position should be meta_csv and -1 position should be trte_partition_file
load_list = [cerebellum_csv, visualCortex_csv, prefrontalCortex_csv, meta_csv, trte_partition_file]
GCN_names = ["cerebellum","visualCortex","prefrontalCortex_csv"]

COMBINER = "VCDN" # Can take values "VCDN", "FullyConnected", or False.
#Use False when only a sinngle GCN output is to be generated and

doSMOTE = True # Boolean

# MORONET training parameters
num_epoch = 30
test_interval = 5
lr = 1e-4
weight_decay = 3e-3
dropout = 0.5
adj_parameter = 5 # average number of edge per node in adj matrix

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