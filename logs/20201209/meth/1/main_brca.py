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
from collections import Counter

from load_data import load_preproc_data, load_trte_partition
from utils import one_hot_tensor, cal_sample_weight 
from models import init_model_dict, init_optim
from train_test import prepare_trte_data, gen_trte_adj_mat
from train_test import train_1_epoch, test_VCDN

from utils import print_dict
from temporaries import findInteractionsSelf, findInteractionsCross, label_specific_acc

# General information about the run to be saved into the log file
print(datetime.datetime.now(),"\n")
RUN_TITLE = "Complete data - DNA Methylation"
RUN_TITLE_SHORT = "meth"
print(RUN_TITLE)
print("This run includes only one GCN - DNA Methylation","\n")
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
    meta_csv = "/data/users/bs16b001/R/TCGA BRCA/PAM50_subtype.csv"
    trte_partition_file = "/data/users/bs16b001/R/TCGA BRCA/trte_partition.txt"
    
    # change label from text to integer
    label_dict = {'Normal':0, 'Basal':1, 'Her2':2, 'LumA':3, 'LumB':4}
    num_class = len(label_dict)
    
    # load preprocessed data from csv
    #load_list - list of csv files to laod other than the meta_csv
    load_list = [meth_preproc_csv] #, meth_preproc_csv, mirna_preproc_csv]
    data_list, labels, patient_id, feat_name_list = load_preproc_data(load_list, meta_csv, label_dict)
    # load tr/te partition
    tr_idx, te_idx = load_trte_partition(trte_partition_file, patient_id)

    # MORONET training parameters
    num_epoch = 2000
    test_inverval = 50
    lr = 1e-3
    GCN_names = ["methylation"] #,"methylation","miRNA"]
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
    '''
    # Cross Interactions
    '''
    mrna_meth= findInteractionsCross(data_tr_list[0], data_tr_list[1])
    data_tr_list.append(mrna_meth[:,:400])
    mrna_meth= findInteractionsCross(data_trte_list[0], data_trte_list[1])
    data_trte_list.append(mrna_meth[:,:400])
    
    mrna_mirna= findInteractionsCross(data_tr_list[0], data_tr_list[2])
    data_tr_list.append(mrna_mirna[:,:400])
    mrna_mirna= findInteractionsCross(data_trte_list[0], data_trte_list[2])
    data_trte_list.append(mrna_mirna[:,:400])

    mirna_meth= findInteractionsCross(data_tr_list[1], data_tr_list[2])
    data_tr_list.append(mirna_meth[:,:400])
    mirna_meth= findInteractionsCross(data_trte_list[1], data_trte_list[2])
    data_trte_list.append(mirna_meth[:,:400])
    '''
    # End Mods

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
    
    # Calculating the dimensions for the GCN and VCDN inputs
    dim_list = [each.shape[1] for each in data_tr_list] # Input dimension of GCN
    dim_he_list = dim_list                              # Hidden layer  dimension of GCN
    dim_hvcdn = num_class**len(dim_list)                # Input dimension of VCDN
           
    # model and optimization
    model_dict = init_model_dict(num_class, dim_list, dim_he_list, dim_hvcdn, gcn_dropout=0.5)

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
            if VERBOSE >= 1 or (epoch == num_epoch):                    #Second condition to print values at the end of run
                print("\nTest: Epoch {:d}".format(epoch))
                print("Train Accuracy: {:.4f}   Test ACC: {:.4f}".format(train_acc, te_acc))
                print("Train F1: {:.4f}         Test F1: {:.4f}".format(train_f1, te_f1))
                print()
            test_ACC.append(te_acc)
            test_F1.append(te_f1)
            train_ACC.append(train_acc)
            train_F1.append(train_f1)

    # Checking label specific accuracies
    # in order to determine if there is an issue due to an unbalanced label distribution
    ground_truth_test = labels_trte[trte_idx["te"]]
    ground_truth_train = labels_trte[trte_idx["tr"]]
    test_accuracies = label_specific_acc(ground_truth_test, te_prob.argmax(1))
    train_accuracies = label_specific_acc(ground_truth_train, train_prob.argmax(1))
    counter_test = Counter(ground_truth_test)
    counter_train = Counter(ground_truth_train)
    print("Here, we have the label specific accuracies:\n")
    print("Label", "\t  ", "Occurrence", "\t", "Training Accuracy", "\t", "Prediction Accuracy")
    print("-"*75)
    for each in label_dict:
        print(each, "\t\t", "{:.4f}".format(counter_test[label_dict[each]]/len(ground_truth_test)), "\t", "{:.4f}".format(train_accuracies[label_dict[each]]), "\t\t", "{:.4f}".format(test_accuracies[label_dict[each]]))

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
    colour = plt.cm.rainbow(np.linspace(0,1,len(GCN_names)+1))
    for i,c in zip(range(len(GCN_names)),colour):
        plt.plot(x, losses[i], c = c, label = "GCN_"+GCN_names[i])
    plt.plot(x,losses[-1], color=colour[-1], label="VCDN")
    plt.yscale("log")
    plt.xlabel("Epochs")
    plt.ylabel("Loss (on log scale)")
    plt.title("MORONET Losses - " + RUN_TITLE_SHORT)
    plt.legend()
    #plt.show()
    time = str(datetime.datetime.now())
    plt.savefig("Losses_"+RUN_TITLE_SHORT+"_"+time+".png")

    test = plt.figure(2)
    x = list(range(0, num_epoch+1, test_inverval))
    plt.plot(x, train_ACC, color="g", label="Train Accuracy")
    plt.plot(x, train_F1, color="m", label="Train F1 Score")
    plt.plot(x, test_ACC, color='b', label="Test Accuracy")
    plt.plot(x, test_F1, color="r", label="Test F1 Score")
    plt.ylim(0,1)
    plt.xlabel("Epochs")
    plt.ylabel("Test Scores")
    plt.title("MORONET Test Evaluation - BRCA Dataset - "+RUN_TITLE_SHORT)
    plt.legend()
    #plt.show()
    plt.savefig("Test_Scores_"+RUN_TITLE_SHORT+"_"+time+".png")
