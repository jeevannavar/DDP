#!/opt/anaconda3/bin/python3.7

from main import train_model
import pandas as pd

RUN_TITLE = "Greedy ANOVA top 1000 \nThis run includes three GCNs - mRNA, DNA Methylation, miRNA"
RUN_TITLE_SHORT = "greedy_anova"
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
load_list = [mrna_preproc_csv, meth_preproc_csv, mirna_preproc_csv, meta_csv, trte_partition_file]
GCN_names = ["mRNA","methylation","miRNA"]

COMBINER = "VCDN" # Can take values "VCDN", "FullyConnected", or False.
#Use False when only a sinngle GCN output is to be generated and

doSMOTE = True # Boolean

# MORONET training parameters
num_epoch = 1000
test_interval = 50
lr = [1e-4, 5e-4, 1e-3, 5e-3, 1e-2]
weight_decay = [5e-4, 1e-3, 3e-3, 5e-3]
dropout = [0.1, 0.25, 0.5, 0.75]
adj_parameter = [2, 5, 8] # average number of edge per node in adj matrix

# VERBOSE setting for print results
VERBOSE = 1 #0, only print final result; 1, only testing result; 2, training and testing result
OUTPUT_FILES = False #Boolean to determine whether to output loss and metrics as csv files
MAKE_PLOTS = True #Boolean to determine whether to output loss and metrics as plots in png format


combinations = [(lr_,weight_decay_,dropout_,adj_parameter_) for lr_ in lr for weight_decay_ in weight_decay for dropout_ in dropout for adj_parameter_ in adj_parameter]

tracker = []
for (lr, weight_decay, dropout, adj_parameter) in combinations:
    print("\n"+"-"*100)
    print("Learning rate : {}\tWeight_decay : {}\tDropout : {}\tAdjacency Parameter : {}".format(lr, weight_decay, dropout, adj_parameter))
    print("-"*100)
    loss_df, metrics_df = train_model(load_list=load_list, label_dict=label_dict, GCN_names=GCN_names, COMBINER=COMBINER,
        SEED=SEED, num_epoch=num_epoch, test_interval=test_interval, lr=lr, weight_decay=weight_decay, 
        dropout=dropout, adj_parameter=adj_parameter, VERBOSE=VERBOSE, doSMOTE = doSMOTE,
        RUN_TITLE=RUN_TITLE, RUN_TITLE_SHORT=RUN_TITLE_SHORT,
        OUTPUT_FILES=OUTPUT_FILES, MAKE_PLOTS=MAKE_PLOTS)
    tracker.append([lr, weight_decay, dropout, adj_parameter, metrics_df.iloc[-1,:]["train_accuracy"], metrics_df.iloc[-1,:]["train_F1"], metrics_df.iloc[-1,:]["test_accuracy"], metrics_df.iloc[-1,:]["test_F1"]])

tracker_df = pd.DataFrame(tracker, columns=["lr","weight_decay","dropout", "adj_parameter", "tr_acc", "tr_F1", "test_acc", "test_F1"])
print(tracker_df)
tracker_df.to_csv("hyperparameters.csv")
