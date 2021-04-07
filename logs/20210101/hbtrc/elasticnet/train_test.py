""" Training and testing of the model
"""
import torch
import torch.nn.functional as F
from utils import gen_adj_mat_tensor, gen_test_adj_mat_tensor, cal_adj_mat_parameter


def prepare_trte_data(data_tensor_list, labels, tr_idx, te_idx):
    """ 
        Re-orgnize train/test data list and scale the data to [0,1]
        Training data is at the top and testing data is among the bottom rows

        Inputs:
            data_tensor_list -> list of torch tensors each of which are inputs for the GCNs
            labels           -> numpy list of ground truth labels
            tr_idx           -> list of indexes for training
            te_idx           -> list of indexes for testing

        Outputs:
            data_tr_list     -> list of torch tensors for GCN inputs for TRAINING ONLY
                                list of len num_input data types. each tensor of shape (N,D)
            data_trte_list   -> list of torch tensors for GCN inputs for trainging and testing
                                list of len num_input data types. each tensor of shape (N,D)
            trte_idx         -> dictionary {"tr":[0 ... num_train]; "te":[num_train ... num_total]}
            labels_new       -> numpy list where initial values are training indices and the later values are test indices
    """
    
    data_tr_list = []
    data_trte_list = []
    for i in range(len(data_tensor_list)):
        data_tr_list.append(data_tensor_list[i][tr_idx].clone())
        data_trte_list.append(torch.cat((data_tensor_list[i][tr_idx].clone(),
                                         data_tensor_list[i][te_idx].clone()),0))
    # scale to [0,1]
    for i in range(len(data_tensor_list)):
        min_tmp = data_tr_list[i].min()
        max_tmp = data_tr_list[i].max()
        data_tr_list[i] = (data_tr_list[i] - min_tmp)/(max_tmp - min_tmp)
        data_trte_list[i] = (data_trte_list[i] - min_tmp)/(max_tmp - min_tmp)
        
    trte_idx = {}
    trte_idx["tr"] = list(range(len(tr_idx)))
    trte_idx["te"] = list(range(len(tr_idx), len(tr_idx)+len(te_idx)))
    
    labels_new = labels.copy()
    labels_new[trte_idx["tr"]] = labels[tr_idx]
    labels_new[trte_idx["te"]] = labels[te_idx]
    
    return data_tr_list, data_trte_list, trte_idx, labels_new         


def gen_trte_adj_mat(data_tr_list, data_trte_list, trte_idx, adj_parameter):
    """ 
        Generate adjacency matrix for GCN training and testing
        Inputs:
            data_tr_list   -> list of torch tensors for GCN inputs for TRAINING ONLY
                                list of len num_input data types. each tensor of shape (N,D)
            data_trte_list -> list of torch tensors for GCN inputs for trainging and testing
                                list of len num_input data types. each tensor of shape (N,D)
            trte_idx       -> dictionary {"tr":[0 ... num_train]; "te":[num_train ... num_total]}
            adj_parameter  -> average number of edges per node in adjacency matrix

        Outputs:
            adj_train_list  -> list of adjacency matrices (as torch tensors) for the training data
            adj_test_list   -> list of adjacency matrices (as torch tensors) for the training and testing data
    """
    adj_metric = "cosine" # cosine distance
    
    adj_train_list = []
    adj_test_list = []
    for i in range(len(data_tr_list)):
        adj_parameter_adaptive = cal_adj_mat_parameter(adj_parameter, data_tr_list[i], adj_metric)
        adj_train_list.append(gen_adj_mat_tensor(data_tr_list[i], adj_parameter_adaptive, adj_metric))
        adj_test_list.append(gen_test_adj_mat_tensor(data_trte_list[i], trte_idx, adj_parameter_adaptive, adj_metric))
    
    return adj_train_list, adj_test_list


def train_1_epoch(data_list, adj_list, label, one_hot_label, sample_weight, model_dict, optim_dict, gcn_names, combiner):
    """
    Train one epoch of GCNs and one epoch of VCDN
    Input:
        data_list       -> list of GCN input tensors, for ex: data_list[0] is view 1 data
        adj_list        -> list of adjacency matrix tensors
        one_hot_label   -> ground truth values as one hot labels
        sample_weight   -> tensor of length num_labels indicating proportion of each label in training data
        model_dict      -> dictionary containing all the individual GCNs and the VCDN
        optim_dict      -> Adam optimiser for each of the GCN and VCDN layers
                            model_dict and optim_dict contain keys in the format:
                                "Ei" for GCN encoder layer i
                                gcn_name for GCN classifier layer
                                combiner for combiner        
        gcn_names       -> list of names for the gcns and the combiner, if it is being used
        combiner        -> name of combiner, false if no combiner is used
    Output:
        loss_dict       -> dict with each layer's name for keys and losses for values
        prob            -> softmax output, i.e., the probability of each label for the train data
                            numpy matrix of shape (N,len(Classes))
    """

    loss_dict = {}
    criterion = torch.nn.CrossEntropyLoss(reduction='none')
    if combiner: classifiers = gcn_names[:-1]
    if not combiner: classifier = gcn_names
    encoders = ["E%d" % i for i in range(1,len(data_list)+1)]

    # Regularization parameter between classifiers
    # The weights are 1 by default. Change here by hardcoding or make this an argument that can be passed to the function
    C_loss_weight = {classifier: 1 for classifier in classifiers}
    if combiner: C_loss_weight[combiner] = 1
    
    modules = model_dict.keys()
    for m in modules:
        model_dict[m].train()    

    ##  Train the GCNs
    ##  View specific encoder and classifier loss

    for classifier in classifiers:
        optim_dict[classifier].zero_grad()
    
    classifier_loss = [None for _ in range(len(data_list))]
    
    for i in range(len(data_list)):
        classifier = model_dict[classifiers[i]](model_dict[encoders[i]](data_list[i],adj_list[i]))
        classifier_loss[i] = torch.mean(torch.mul(torch.mean(torch.mul(classifier-one_hot_label, classifier-one_hot_label),1), sample_weight))*C_loss_weight[classifiers[i]]
        classifier_loss[i].backward()
        optim_dict[classifiers[i]].step()
        loss_dict[classifiers[i]] = classifier_loss[i].detach().cpu().numpy().item()
    
    ##  Train the VCDN
    ##  View integration classifier loss

    if combiner:
        optim_dict[combiner].zero_grad()
        c_loss = 0
        views = []
        for i in range(len(data_list)):
            views.append(model_dict[classifiers[i]](model_dict[encoders[i]](data_list[i],adj_list[i])))

        c = model_dict[combiner](views)    
        c_loss = torch.mean(torch.mul(criterion(c, label),sample_weight))*C_loss_weight[combiner]
        c_loss.backward()
        optim_dict[combiner].step()
        loss_dict[combiner] = c_loss.detach().cpu().numpy().item()

    else:
        c = classifier
    # predictionprobability N by C matrix
    prob = F.softmax(c,dim=1).data.cpu().numpy()
    
    return loss_dict, prob
    

def test_VCDN(data_list, adj_list, te_idx, model_dict, gcn_names, combiner):
    """ 
    Test the model
    Input:
        data_list   -> list of GCN input tensors, for ex: data_list[0] is view 1 data
        adj_lis     -> list of adjacency matrix tensors
        te_idx      -> data idx for testing
        model_dict  -> dictionary containing all the individual GCNs and the VCDN
        gcn_names   -> list of namesfor the gcns and the combiner, if it is being used
        combiner    -> name of combiner, false if no combiner is used
    Output:
        prob    -> softmax output, i.e., the probability of each label for the test data
                    numpy matrix of shape (N,len(Classes))
    """
    modules = model_dict.keys()
    for m in modules:
        model_dict[m].eval()
    
    if combiner: classifiers = gcn_names[:-1]
    if not combiner: classifiers = gcn_names
    encoders = ["E%d" % i for i in range(1,len(data_list)+1)]

    gcns = [None for _ in range(len(data_list))]
    for i in range(len(data_list)):
        gcns[i] = model_dict[classifiers[i]](model_dict[encoders[i]](data_list[i],adj_list[i]))
    if combiner:
        c = model_dict[combiner](gcns)
    else:
        c = gcns[0]
    c = c[te_idx,:]
    
    # prediction probability N by C matrix
    prob = F.softmax(c, dim=1).data.cpu().numpy()
    
    return prob