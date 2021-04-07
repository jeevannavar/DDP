#!/opt/anaconda3/bin/python3.7

""" Example for the BRCA dataset
"""
import numpy as np
import random
from sklearn.metrics import accuracy_score, f1_score
import matplotlib.pyplot as plt
import datetime
import torch
import csv

from load_data import load_preproc_data, load_trte_partition
from utils import one_hot_tensor, cal_sample_weight 
from models import init_model_dict, init_optim
from train_test import prepare_trte_data, gen_trte_adj_mat
from train_test import train_1_epoch, test_VCDN

from utils import print_dict
from temporaries import findInteractionsSelf, findInteractionsCross

# General information about the run to be saved into the log file
print(datetime.datetime.now(),"\n")
print("This run includes GCNs - mRNA, DNA methylation, and miRNA - each with non-preprocessed data","\n")

# assign GPU number
CUDA_DEVICE = 0
cuda = True if torch.cuda.is_available() else False

if cuda:
    Print("CUDA Device in use")

# Set random seed
#SEED = 123
#random.seed(SEED)
#np.random.seed(SEED)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False
#torch.manual_seed(SEED)
if cuda:
    #torch.cuda.manual_seed(SEED)
    torch.cuda.set_device(CUDA_DEVICE)


if __name__ == "__main__":    
    # pre-processed data
    mrna_preproc_csv = "/data/users/bs16b001/R/TCGA BRCA/mrna_complete.csv"
    meth_preproc_csv = "/data/users/bs16b001/R/TCGA BRCA/meth_complete.csv"
    mirna_preproc_csv = "/data/users/bs16b001/R/TCGA BRCA/mirna_complete.csv"
    mirna_mrna_csv = "BRCA/BRCA_mrna_mirna_interaction.csv"
    meta_csv = "/data/users/bs16b001/R/TCGA BRCA/PAM50_subtype.csv"
    trte_partition_file = "/data/users/bs16b001/R/TCGA BRCA/trte_partition.txt"
    
    # change label from text to integer
    label_dict = {'Normal':0, 'Basal':1, 'Her2':2, 'LumA':3, 'LumB':4}
    num_class = len(label_dict)
    
    # load preprocessed data from csv
    #load_list - list of csv files to laod other than the meta_csv
    load_list = [mrna_preproc_csv, meth_preproc_csv, mirna_preproc_csv]
    data_list, labels, patient_id, feat_name_list = load_preproc_data(load_list, meta_csv, label_dict)
    # load tr/te partition
    tr_idx, te_idx = load_trte_partition(trte_partition_file, patient_id)
    
    # MORONET training parameters
    num_epoch = 2000
    test_inverval = 50
    lr = 1e-3
    GCN_names = ["mRNA","methylation","miRNA"]
    #dim_list = [200,200,200,40000,40000,40000] # Input dimension of GCN
    #dim_he_list = [200,200,200,40000,40000,40000] # Hidden layer dimension of GCN
    #dim_hvcdn = num_class**len(dim_list) # Input dimension of VCDN
    adj_parameter = 5 # average number of edge per node in adj matrix
    
    # VERBOSE setting for print results
    VERBOSE = 2 #0, only print final result; 1, only testing result; 2, training and testing result
    
    """ Prepare data and model for training and testing
    """
    # data to tensor
    data_tensor_list = []
    for i in range(len(data_list)):
        data_tensor_list.append(torch.FloatTensor(data_list[i]))
        if cuda:
            data_tensor_list[i] = data_tensor_list[i].cuda()
    
    # generate training and testing data
    data_tr_list, data_trte_list, trte_idx, labels_trte = prepare_trte_data(data_tensor_list, labels, tr_idx, te_idx)
    '''
    #Mods
    #For training data
    concat = torch.cat((data_tr_list[0], data_tr_list[1], data_tr_list[2]),1)
    data_tr_list.append(concat)
    #For testing data
    concat = torch.cat((data_trte_list[0], data_trte_list[1], data_trte_list[2]),1)
    data_trte_list.append(concat)
    
    # Mods
    # mRNA - SELF
    # For training data
    mrna_interaction = findInteractions(data_tr_list[0])
    if cuda:    
        mrna_interaction.cuda()
    data_tr_list.append(mrna_interaction)

    # For training and testing data
    mrna_interaction = findInteractions(data_trte_list[0])
    if cuda:    
        mrna_interaction.cuda()
    data_trte_list.append(mrna_interaction)
    #
    
    # Cross Interactions
    mrna_meth= findInteractionsCross(data_tr_list[0], data_tr_list[1])
    data_tr_list.append(mrna_meth)
    mrna_meth= findInteractionsCross(data_trte_list[0], data_trte_list[1])
    data_trte_list.append(mrna_meth)
    
    mrna_mirna= findInteractionsCross(data_tr_list[0], data_tr_list[2])
    data_tr_list.append(mrna_mirna)
    mrna_mirna= findInteractionsCross(data_trte_list[0], data_trte_list[2])
    data_trte_list.append(mrna_mirna)

    mirna_meth= findInteractionsCross(data_tr_list[1], data_tr_list[2])
    data_tr_list.append(mirna_meth)
    mirna_meth= findInteractionsCross(data_trte_list[1], data_trte_list[2])
    data_trte_list.append(mirna_meth)
    

    # End Mods
    '''
    # calculating the dimensions for the GCN and VCDN inputs
    dim_list = dim_he_list = [each.shape[1] for each in data_tr_list]
    dim_hvcdn = num_class**len(dim_list)

    labels_tr_tensor = torch.LongTensor(labels_trte[trte_idx["tr"]])
    onehot_labels_tr_tensor = one_hot_tensor(labels_tr_tensor, num_class)
    sample_weight_tr = cal_sample_weight(labels_trte[trte_idx["tr"]], num_class)
    sample_weight_tr =  torch.FloatTensor(sample_weight_tr)
    if cuda:
        labels_tr_tensor = labels_tr_tensor.cuda()
        onehot_labels_tr_tensor = onehot_labels_tr_tensor.cuda()
        sample_weight_tr = sample_weight_tr.cuda()
    
    # calculate adjacency matrix
    adj_tr_list, adj_te_list = gen_trte_adj_mat(data_tr_list, data_trte_list, trte_idx, adj_parameter)
        
    # model and optimization
    model_dict = init_model_dict(num_class, dim_list, dim_he_list, dim_hvcdn)

    for m in model_dict:
        if cuda:
            model_dict[m].cuda()
    optim_dict = init_optim(model_dict, dim_list, lr)  
    
    
    """ Training and testing
    """
    losses = [[] for _ in range(len(dim_list)+1)] 
    test_ACC = []
    test_F1 = [] 
    train_ACC = []
    train_F1 = []
    print("\nTraining...")
    for epoch in range(num_epoch+1):
        """ Train one epoch
        """
        tr_loss_dict, train_prob = train_1_epoch(data_tr_list, adj_tr_list, labels_tr_tensor, 
                                           onehot_labels_tr_tensor, sample_weight_tr, model_dict, optim_dict)
        if VERBOSE >= 2:
            print("Epoch {:d}, Train loss: {:}".format(epoch, print_dict(tr_loss_dict, print_format=".2e")))
        
        classifiers = ["C%d" % i for i in range(1,len(data_tr_list)+1)]
        for i in range(len(classifiers)):
            losses[i].append(tr_loss_dict[classifiers[i]])
        losses[-1].append(tr_loss_dict["C"])
        
        """ Testing
        """
        if epoch % test_inverval == 0:
            te_prob = test_VCDN(data_trte_list, adj_te_list, trte_idx["te"], model_dict)
            te_acc = accuracy_score(labels_trte[trte_idx["te"]], te_prob.argmax(1))
            te_f1 = f1_score(labels_trte[trte_idx["te"]], te_prob.argmax(1), average='weighted')
            train_acc = accuracy_score(labels_trte[trte_idx["tr"]], train_prob.argmax(1))
            train_f1 = f1_score(labels_trte[trte_idx["tr"]], train_prob.argmax(1), average='weighted')
            if VERBOSE >= 1 or (epoch == num_epoch):
                print("\nTest: Epoch {:d}".format(epoch))
                print("Train Accuracy: {:.4f}   Test ACC: {:.4f}".format(train_acc, te_acc))
                print("Train F1: {:.4f}         Test F1: {:.4f}".format(train_f1, te_f1))
                print()
            test_ACC.append(te_acc)
            test_F1.append(te_f1)
            train_ACC.append(train_acc)
            train_F1.append(train_f1)

    # Output the necessary data as a csv files
    x = list(range(num_epoch+1))
    loss_with_index = [x,*losses]
    transpose_loss = list(map(list, zip(*loss_with_index)))
    transpose_loss_labelled = [["index"]+GCN_names+["VCDN"], *transpose_loss]
    with open('loss.csv', 'w', newline='') as f:
        losswriter = csv.writer(f, delimiter=",")
        losswriter.writerows(transpose_loss_labelled)

    x = list(range(0, num_epoch+1, test_inverval))
    acc_with_index = [x,train_ACC, train_F1, test_ACC, test_F1]
    transpose_acc = list(map(list, zip(*acc_with_index)))
    transpose_acc_labelled = [["index","train_accuracy","train_F1","test_accuracy","test_F1"], *transpose_acc]
    with open('metrics.csv', 'w', newline='') as f:
        metricwriter = csv.writer(f, delimiter=",")
        metricwriter.writerows(transpose_acc_labelled)

    # Plotting Figures
    loss = plt.figure(1)
    x = list(range(num_epoch+1))
    plt.plot(x,losses[0], color='b', label="GCN_mRNA")
    plt.plot(x,losses[1], color='g', label="GCN_methylation")
    plt.plot(x,losses[2], color='r', label="GCN_miRNA")
    #plt.plot(x,losses[3], color='m', label="GCN_mRNA-meth_Naive")
    #plt.plot(x,losses[4], color='y', label="GCN_mRNA-miRNA_Naive")
    #plt.plot(x,losses[5], color='k', label="GCN_miRNA-meth_Naive")
    plt.plot(x,losses[-1], color='c', label="VCDN")
    plt.yscale("log")
    plt.xlabel("Epochs")
    plt.ylabel("Loss (on log scale)")
    plt.title("MORONET Losses - Complete Data")
    plt.legend()
    #plt.show()
    time = str(datetime.datetime.now())
    plt.savefig("Losses_complete_"+time+".png")

    test = plt.figure(2)
    x = list(range(0, num_epoch+1, test_inverval))
    plt.plot(x, train_ACC, color="g", label="Train Accuracy")
    plt.plot(x, train_F1, color="m", label="Train F1 Score")
    plt.plot(x, test_ACC, color='b', label="Test Accuracy")
    plt.plot(x, test_F1, color="r", label="Test F1 Score")
    plt.ylim(0,1)
    plt.xlabel("Epochs")
    plt.ylabel("Test Scores")
    plt.title("MORONET Test Evaluation - BRCA Dataset - Complete Data")
    plt.legend()
    #plt.show()
    plt.savefig("Test_Scores_complete_"+time+".png")
