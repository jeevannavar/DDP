#!/opt/anaconda3/bin/python3.7

"""
MAIN FUNCTION
"""
import numpy as np
import random
from sklearn.metrics import accuracy_score, f1_score
import matplotlib.pyplot as plt
import datetime
import torch
import csv
from imblearn.over_sampling import SMOTE
from collections import Counter
import pandas as pd

from load_data import load_preproc_data, load_trte_partition
from utils import one_hot_tensor, cal_sample_weight 
from models import init_model_dict, init_optim
from train_test import prepare_trte_data, gen_trte_adj_mat
from train_test import train_1_epoch, test_VCDN

from utils import print_dict
from temporaries import findInteractionsSelf, findInteractionsCross, label_specific_acc

def train_model(load_list, label_dict, GCN_names, COMBINER, SEED="random", num_epoch=1000, test_interval=50, lr=1e-4, weight_decay=3e-3, dropout=0.5, adj_parameter=5, VERBOSE=2, doSMOTE=True, RUN_TITLE="", RUN_TITLE_SHORT="", OUTPUT_FILES=False, MAKE_PLOTS=False):
    # General information about the run to be saved into the log file
    print(datetime.datetime.now(),"\n")
    print(RUN_TITLE)
    if COMBINER: GCN_names.append(COMBINER)
    
    # assign GPU number
    CUDA_DEVICE = 0
    cuda = True if torch.cuda.is_available() else False

    if cuda:
        Print("CUDA Device in use")

    # Set random seed
    if SEED == "random":
        SEED = random.randint(0,100000)
    print("SEED = ", SEED)
    random.seed(SEED)
    np.random.seed(SEED)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    torch.manual_seed(SEED)
    if cuda:
        torch.cuda.manual_seed(SEED)
        torch.cuda.set_device(CUDA_DEVICE)
    
    # load preprocessed data from csv

    num_class = len(label_dict)
    data_list, labels, patient_id, feat_name_list = load_preproc_data(load_list[:-2], load_list[-2], label_dict)
    # load tr/te partition
    tr_idx, te_idx = load_trte_partition(load_list[-1], patient_id)

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


    labels_tr_tensor = torch.LongTensor(labels_trte[trte_idx["tr"]])

    # Here we modify the training data set using SMOTE. We have only changed the training data.
    # While testing, we use training data, but only the real ones, not the synthesised ones
    if doSMOTE:
        sm = SMOTE(random_state = SEED)
        for i in range(len(data_tr_list)):
            data_tr_list[i], labels = sm.fit_sample(data_tr_list[i], labels_tr_tensor)
            data_tr_list[i] = torch.FloatTensor(data_tr_list[i])
        labels_tr_tensor = torch.LongTensor(labels)
        onehot_labels_tr_tensor = one_hot_tensor(labels_tr_tensor, num_class)

        sample_weight_tr = cal_sample_weight(labels_tr_tensor.numpy(), num_class)
        sample_weight_tr =  torch.FloatTensor(sample_weight_tr)
    else:
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
    model_dict = init_model_dict(num_class, dim_list, dim_he_list, dim_hvcdn, GCN_names, COMBINER, gcn_dropout=dropout)

    for m in model_dict:
        if cuda:
            model_dict[m].cuda()
    optim_dict = init_optim(model_dict, dim_list, GCN_names, COMBINER, lr, weight_decay)  
    
    #Some necessary lists
    ground_truth_test = labels_trte[trte_idx["te"]]
    ground_truth_train = labels_tr_tensor.numpy()
    counter_test = Counter(ground_truth_test)
    counter_train = Counter(ground_truth_train)
    train_accuracies_dict = {key:[] for key in label_dict}
    test_accuracies_dict = {key:[] for key in label_dict}
    
    
    """ Training and testing
    """
    losses = [[] for _ in range(len(GCN_names))] 
    test_ACC = []
    test_F1 = [] 
    train_ACC = []
    train_F1 = []
    print("\nTraining...")
    for epoch in range(num_epoch+1):
        """ Train one epoch
        """
        tr_loss_dict, train_prob = train_1_epoch(data_tr_list, adj_tr_list, labels_tr_tensor, 
                                           onehot_labels_tr_tensor, sample_weight_tr, model_dict, optim_dict,
                                           GCN_names, COMBINER)
        if VERBOSE >= 2:
            print("Epoch {:d}, Train loss: {:}".format(epoch, print_dict(tr_loss_dict, print_format=".2e")))
        
        for i in range(len(GCN_names)):
            losses[i].append(tr_loss_dict[GCN_names[i]])
        
        """ Testing
        """
        if epoch % test_interval == 0:
            te_prob = test_VCDN(data_trte_list, adj_te_list, trte_idx["te"], model_dict, GCN_names, COMBINER)
            te_acc = accuracy_score(labels_trte[trte_idx["te"]], te_prob.argmax(1))
            te_f1 = f1_score(labels_trte[trte_idx["te"]], te_prob.argmax(1), average='weighted')
            train_acc = accuracy_score(labels_tr_tensor, train_prob.argmax(1))
            train_f1 = f1_score(labels_tr_tensor, train_prob.argmax(1), average='weighted')
            if VERBOSE >= 1 or (epoch == num_epoch):                    #Second condition to print values at the end of run
                print("\nTest: Epoch {:d}".format(epoch))
                print("Train Accuracy: {:.4f}   Test ACC: {:.4f}".format(train_acc, te_acc))
                print("Train F1: {:.4f}         Test F1: {:.4f}".format(train_f1, te_f1))
                print()
            
                # Checking label specific accuracies
                # in order to determine if there is an issue due to an unbalanced label distribution
                accuracies_test = label_specific_acc(ground_truth_test, te_prob.argmax(1))
                accuracies_train = label_specific_acc(ground_truth_train, train_prob.argmax(1))
                #print("Here, we have the label specific accuracies:\n")
                print("Label", "      ", "Train Distribution", "   ", "Train Accuracy", "   ", "Test Distribution", "   ", "Test Accuracy")
                print("-"*90)
                for each in label_dict:
                    dist_train = "{:.4f}".format(counter_train[label_dict[each]]/len(ground_truth_train))
                    dist_test = "{:.4f}".format(counter_test[label_dict[each]]/len(ground_truth_test))
                    acc_train = "{:.4f}".format(accuracies_train[label_dict[each]])
                    acc_test = "{:.4f}".format(accuracies_test[label_dict[each]])
                    print(each, "\t          ", dist_train, "            ", acc_train, "            ", dist_test, "            ", acc_test)
                print() #Just a spacer. Need the extra blank line

            #Storing all metrics at test intervals
            accuracies_test = label_specific_acc(ground_truth_test, te_prob.argmax(1))
            accuracies_train = label_specific_acc(ground_truth_train, train_prob.argmax(1))
            for each in label_dict:
                train_accuracies_dict[each].append(accuracies_train[label_dict[each]])
                test_accuracies_dict[each].append(accuracies_train[label_dict[each]])
            test_ACC.append(te_acc)
            test_F1.append(te_f1)
            train_ACC.append(train_acc)
            train_F1.append(train_f1)


    # Output the necessary data as a csv files
    losses_df = pd.DataFrame(zip(*losses), columns=GCN_names)
    if OUTPUT_FILES:
        losses_df.to_csv("loss.csv")


    metrics = [train_ACC, train_F1, test_ACC, test_F1]+ list(train_accuracies_dict.values()) + list(test_accuracies_dict.values())
    columns = ["train_accuracy","train_F1","test_accuracy","test_F1"] + ["train_acc_"+key for key in label_dict] + ["test_acc_"+key for key in label_dict]
    metrics_df = pd.DataFrame(zip(*metrics), columns=columns, index=range(0,num_epoch+1,test_interval))
    if OUTPUT_FILES:
        metrics_df.to_csv("metrics.csv")

    # Plotting Figures
    if MAKE_PLOTS:
        loss = plt.figure(1)
        x = list(range(num_epoch+1))
        colour = plt.cm.rainbow(np.linspace(0,1,len(GCN_names)+1))
        for i,c in zip(range(len(GCN_names)),colour):
            plt.plot(x, losses_df.iloc[:,i], c = c, label = GCN_names[i])
        plt.yscale("log")
        plt.xlabel("Epochs")
        plt.ylabel("Loss (on log scale)")
        plt.title("MORONET Losses - " + RUN_TITLE_SHORT)
        plt.legend()
        #plt.show()
        time = str(datetime.datetime.now())
        plt.savefig("Losses_"+RUN_TITLE_SHORT+"_"+time+".png")

        test = plt.figure(2)
        x = list(range(0, num_epoch+1, test_interval))
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

    return losses_df, metrics_df