"""
Dividing up the "train_model" function into two chunks. 
One for data processing and another for model training.
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
from utils import one_hot_tensor, cal_sample_weight, cal_adj_mat_parameter, gen_adj_mat_tensor
from models import init_model_dict, init_optim
from train_test import prepare_trte_data, gen_trte_adj_mat
from train_test import train_1_epoch, test_VCDN

from utils import print_dict
from temporaries import findInteractionsSelf, findInteractionsCross, label_specific_acc

def load_and_process(load_list, label_dict, doSMOTE=True, SEED="random"):
    # assign GPU number
    CUDA_DEVICE = 0
    cuda = True if torch.cuda.is_available() else False

    if cuda:
        print("CUDA Device in use")
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

    return (data_tr_list, data_trte_list, trte_idx, labels_tr_tensor, onehot_labels_tr_tensor, labels_trte, sample_weight_tr)



def train(data, label_dict, GCN_names, COMBINER, SEED="random", num_epoch=1000, test_interval=50, lr=1e-4, weight_decay=3e-3, dropout=0.5, adj_parameter=5, VERBOSE=2, RUN_TITLE="", RUN_TITLE_SHORT="", OUTPUT_FILES=False, feature_extract=[]):
    
    # loading processed data
    (data_tr_list, data_trte_list, trte_idx, labels_tr_tensor, onehot_labels_tr_tensor, labels_trte, sample_weight_tr) = data

    # General information about the run to be saved into the log file or just on the console
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
    
    num_class = onehot_labels_tr_tensor.shape[1]

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
                test_accuracies_dict[each].append(accuracies_test[label_dict[each]])
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


    # Feature Extraction
    feature_imp = {}
    
    # SHAP
    if "shap" in feature_extract:
        print("Performing feature extraction using SHAP ...\n")
        import shap

        def pred_funct(sample_numpy):
            '''
            Required for SHAP. SHAP alters feature values and uses this to check model behaviour
            This function inherits many of the variables from the outer function train_model()
            Input:
                sample_numpy = numpy 2d array containing perturbed samples of size (samples, features)
            Output:
                predictions  = softmax output, i.e., the probability of each label for the test data numpy matrix of shape (N,len(Classes))
            '''
            sample_torch = torch.FloatTensor(sample_numpy)
            indices = list(range(sample_torch.shape[0]))
            
            lengths = [each.shape[1] for each in data_tr_list]
            sample_list = []
            adj_list = []
            end = 0
            for i in range(len(lengths)):
                start = end
                end = start + lengths[i]
                sample_list.append(sample_torch[:,start:end])
                adj_parameter_adaptive = cal_adj_mat_parameter(adj_parameter, data_tr_list[i])
                adj_list.append(gen_adj_mat_tensor(sample_list[i], adj_parameter_adaptive))
            
            predictions = test_VCDN(sample_list, adj_list, indices, model_dict, GCN_names, COMBINER)
            
            return predictions

        train_data_full = np.array(torch.cat(tuple(each for each in data_tr_list), dim=1))
        kmeans_data = shap.kmeans(train_data_full, k=1)
        explainer = shap.KernelExplainer(pred_funct, kmeans_data)
        data_test = np.array(torch.cat(tuple(each for each in data_trte_list), dim=1))
        shap_values = explainer.shap_values(data_test)
        feature_names = [name for sublist in feat_name_list for name in sublist]
        #shap.summary_plot(shap_values = shap_values, feature_names=feature_names, max_display=40, plot_type='bar', class_names=list(label_dict.keys()))

        shap_list = np.sum(np.array([np.mean(abs(each), axis=0) for each in shap_values]), axis=0)
        shap_df = pd.DataFrame({"features":feature_names, "shapley_values":shap_list})
        shap_df.sort_values(by="features", ascending=True).set_index("features")
        #This is a pandas dataframe of all features with their shapley values averaged across samples and summed across classes
        feature_imp["shap"] = shap_df
        shap_df.to_csv("shapley_values.csv")

        raw_shapley_values = np.array(shap_values)
        classes = label_dict.keys()
        samples = patient_id
        features = feature_names
        names = ['classes', 'samples', 'features']
        index = pd.MultiIndex.from_product([range(s)for s in raw_shapley_values.shape], names=names)
        #index.set_names(names, inplace=True)
        time = str(datetime.datetime.now())
        raw_shapley_values = pd.DataFrame({'Shap_values': raw_shapley_values.flatten()}, index=index)['Shap_values']
        raw_shapley_values.index = pd.MultiIndex.from_tuples([(x,y,z) for x in classes for y in samples for z in features])
        raw_shapley_values.to_csv("raw_shapley_values_"+RUN_TITLE_SHORT+"_"+time+".csv", index_label=names)

    ## LIME
    if "lime" in feature_extract:
        print("Performing feature extraction using LIME ...\n")
        from lime import lime_tabular
        from lime import submodular_pick
        import re

        def pred_funct(sample_numpy):
            '''
            Required for LIME. LIME perturbs feature values and uses this to check model behaviour
            This function inherits many of the variables from the outer function train_model()
            Input:
                sample_numpy = numpy 2d array containing perturbed samples of size (samples, features)
            Output:
                predictions  = softmax output, i.e., the probability of each label for the test data numpy matrix of shape (N,len(Classes))
            '''
            sample_torch = torch.FloatTensor(sample_numpy)
            if len(sample_torch)==1: sample_torch = sample_torch.unsqueeze(0)
            
            data_full = torch.cat(tuple(each for each in data_tr_list), dim=1)
            data_combined = torch.cat((data_full, sample_torch), dim=0)
            indices = list(range(data_full.shape[0], data_combined.shape[0]))
            
            lengths = [each.shape[1] for each in data_tr_list]
            data_combined_list = []
            adj_list = []
            end = 0
            for i in range(len(lengths)):
                start = end
                end = start + lengths[i]
                data_combined_list.append(data_combined[:,start:end])
                adj_parameter_adaptive = cal_adj_mat_parameter(adj_parameter, data_tr_list[i])
                adj_list.append(gen_adj_mat_tensor(data_combined_list[i], adj_parameter_adaptive))
            
            return test_VCDN(data_combined_list, adj_list, indices, model_dict, GCN_names, COMBINER)

        data_full = np.array(torch.cat(tuple(each for each in data_trte_list), dim=1))
        feature_names = [name for sublist in feat_name_list for name in sublist]
        explainer = lime_tabular.LimeTabularExplainer(data_full, feature_names=feature_names, class_names=label_dict.values(), discretize_continuous=True)
        sp_obj = submodular_pick.SubmodularPick(explainer, data_full, pred_funct, method="full", num_features=min(2000, data_full.shape[1]), num_exps_desired=100)
        # The SubmodularPick function takes two arguments for method: "full" and "sample"
        # If "sample" method is used, the sample_size argument has to be mentioned which defaults to 1000. More samples == more computation time

        #Making a dataframe of all the explanations of sampled points
        W = pd.DataFrame([dict(each.as_list(each.available_labels()[0])) for each in sp_obj.explanations]).fillna(0)
        W['prediction'] = [each.available_labels()[0] for each in sp_obj.explanations]

        importance = W.groupby("prediction").mean().T   #Mean of lime coefficients by class
        importance["aggregate"] = abs(importance).sum(axis=1) #Calculating aggregate importance

        new_index = [re.search("(.*)[\>\<]\d*", name).group(1).strip().split(' ')[-1] for name in importance.index]

        importance = importance.rename(index=dict(zip(importance.index,new_index)))["aggregate"]
        importance = importance.groupby(level=0).sum() 
        #This is a pandas series of all features with their importance score according to lime
        feature_imp["lime"] = importance

        #LIME score quantifies how well the lime was able to 
        average_lime_score = np.mean([each.score for each in sp_obj.explanations])
        std_dev = np.std([each.score for each in sp_obj.explanations])
        print("Average lime score = ", average_lime_score, "\n", "Standard deviation = ", std_dev, "\n")

    return losses_df, metrics_df, feature_imp, model_dict


def process(data_list, labels, indices, SEED="random", doSMOTE=True):
    """
    Inputs:
        data_list = list of pandas DataFrames representing individual omics
        labels = list of labels in the form of ints
        indices = (training_indices, test_indices), tuple of lists

    """
    if isinstance(data_list[0], pd.DataFrame):
        patient_id = data_list[0].index.tolist()
        feat_name_list = [None]*len(data_list)
        for i in range(len(data_list)):
            feat_name_list[i] = data_list[i].columns.to_numpy()
            data_list[i] = data_list[i].to_numpy()
    else:
        sys.exit("Wrond data type. Supply list of 'pandas.DataFrame's.")
    
    (tr_idx, te_idx) = indices
    
    num_class = len(set(labels))
    
    # assign GPU number
    CUDA_DEVICE = 0
    cuda = True if torch.cuda.is_available() else False

    if cuda:
        print("CUDA Device in use")
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

    return (data_tr_list, data_trte_list, trte_idx, labels_tr_tensor, onehot_labels_tr_tensor, labels_trte, sample_weight_tr)



if __name__ == '__main__':
    mirna = "R/TCGA BRCA/mirna_anova.csv"
    mrna = "R/TCGA BRCA/mrna_top1000.csv"
    meta_csv = "R/TCGA BRCA/PAM50_subtype.csv"
    trte_partition_file = "R/TCGA BRCA/trte_partition.txt"

    label_dict = {'Normal':0, 'Basal':1, 'Her2':2, 'LumA':3, 'LumB':4}

    load_list = [mirna, mrna, meta_csv, trte_partition_file]
    GCN_names = ["mirna", "mrna"]

    COMBINER = False

    #data = load_and_process(load_list, label_dict)
    #loss, metrics, imp, model = train(data, GCN_names, COMBINER, num_epoch=20)

    tr, te = process(load_list, label_dict)
    print(len(tr), len(te))
    # data_list = list of numpy.ndarrays of shape (num_samples, num_features)
    # labels = list of ints indicating label; len = num_samples
    # patient_id = list of patient ids; len = num_samples
    # feat_name_list = list of 1d numpy arrays with feature names
    # tr_idx = list of training indices
    # te_idx = list of test indices
