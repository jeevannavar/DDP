""" Load and pre-process data
"""
import numpy as np
import pandas as pd
from sklearn.feature_selection import SelectKBest, f_classif


def load_data_csv(csv_file, meta_patient_id, take_log=False):
    """ csv_file: columns: patient_id, features
        take_log: whether take log of the data
        meta_patient_id: patient_id order from metadata
    """
    df = pd.read_csv(csv_file, delimiter=',')
    patient_id = df["patient_id"].values.tolist()
    # check if patient_id in the same order as meta_patient_id
    assert len(patient_id) == len(meta_patient_id), "patient number different from metadata"
    for p, meta_p in zip(patient_id, list(meta_patient_id)):
        assert p == meta_p, "patient order different from metadata"
    feat_name = np.array(list(df.columns)[1:])
    data = df.iloc[:,1:].values
    if take_log:
        data = np.log2(data+1)
    
    return data, feat_name


def load_metadata(meta_csv, label_dict):
    """ meta_csv: metadata csv file
        label_dict: label text to integer
    """
    meta = pd.read_csv(meta_csv, delimiter=',')
    patient_id = meta.iloc[:,0].values
    labels = meta.iloc[:,1].values
    # change label from text to integer
    new_labels = np.zeros(len(labels)).astype(int)
    for k in label_dict:
        idx = np.where(labels == k)[0]
        new_labels[idx] = label_dict[k]
    labels = new_labels.copy()
    
    return labels, patient_id
    

def load_data(mrna_csv, meth_csv, mirna_csv, meta_csv, label_dict):
    """ label_dict: label text to integer 
    """
    labels, meta_patient_id = load_metadata(meta_csv, label_dict)
    mrna_data, mrna_feat_name = load_data_csv(mrna_csv, meta_patient_id, take_log=True)
    meth_data, meth_feat_name = load_data_csv(meth_csv, meta_patient_id, take_log=False)
    mirna_data, mirna_feat_name = load_data_csv(mirna_csv, meta_patient_id, take_log=True)
    
    data_list = [mrna_data, meth_data, mirna_data]
    feat_name_list = [mrna_feat_name, meth_feat_name, mirna_feat_name]
    
    return data_list, labels, meta_patient_id, feat_name_list


def load_trte_partition(trte_partition_file, patient_id):
    """ trte_partition_file
        lines: "training patient_id"; training patient_ids; "testing patient_id"; test patient_ids
        patient_id: list of patient_id indicating order
        return tr_idx, te_idx
    """
    f = open(trte_partition_file, 'r') 
    lines = f.readlines()
    f.close()
    tr_patient_id_list = lines[1].strip().split(',')
    te_patient_id_list = lines[3].strip().split(',')
    
    tr_idx = [np.where(patient_id == pid)[0][0] for pid in tr_patient_id_list]
    te_idx = [np.where(patient_id == pid)[0][0] for pid in te_patient_id_list]
    
    return tr_idx, te_idx
    

def pre_proc_mean_var(data, mean_min=0, mean_max=float("inf"), var_min=0, var_max=float("inf")):
    """ features pre-selection based on mean and variance
        data: N by P
    """
    data_mean = np.mean(data, axis=0)
    data_var = np.var(data, axis=0)
    idx_mean_min = np.where(data_mean>mean_min)[0]
    idx_mean_max = np.where(data_mean<mean_max)[0]
    idx_mean = np.intersect1d(idx_mean_min,idx_mean_max)
    idx_var_min = np.where(data_var>var_min)[0]
    idx_var_max = np.where(data_var<var_max)[0]
    idx_var = np.intersect1d(idx_var_min,idx_var_max)
    idx = np.intersect1d(idx_mean,idx_var)
    
    return idx


def pre_processing(train_data, labels, mean_min, mean_max, var_min, var_max, retain_dim):
    """ preprocessing and feature pre-selection
        train_data: N by P
        labels
        mean_min, mean_max, var_min, var_max: for features pre-selection based on mean and variance
        retain_dim: retained number of features after pre-processing
        return
        selected feature idx
    """
    idx_meanvar = pre_proc_mean_var(train_data, mean_min, mean_max, var_min, var_max)
    idx_feat_pre = SelectKBest(f_classif,k=retain_dim).fit(train_data[:,idx_meanvar],labels).get_support(indices=True)
    idx_pre = idx_meanvar[idx_feat_pre]
    
    return idx_pre

    
def save_data_to_csv(csv_file, patient_id, data, feat_name):
    df = pd.DataFrame(data, columns=list(feat_name))
    df.insert(0,"patient_id",patient_id)
    df.to_csv(csv_file, index=False, header=True)


def load_preproc_data(load_list, meta_csv, label_dict):
    """ label_dict: label text to integer 
    """
    labels, meta_patient_id = load_metadata(meta_csv, label_dict)
    data_list = []
    feat_name_list = []
    for each in load_list:
        data, feat_name = load_data_csv(each, meta_patient_id)
        data_list.append(data)
        feat_name_list.append(feat_name)
    
    return data_list, labels, meta_patient_id, feat_name_list


if __name__ == "__main__":
    # original data
    mrna_csv = "BRCA/BRCA_mrna.csv"
    meth_csv = "BRCA/BRCA_meth.csv"
    mirna_csv = "BRCA/BRCA_mirna.csv"
    meta_csv = "BRCA/BRCA_PAM50_subtype.csv"
    trte_partition_file = "BRCA/trte_partition.txt"
    
    # pre-processed data
    mrna_preproc_csv = "BRCA/BRCA_mrna_preproc.csv"
    meth_preproc_csv = "BRCA/BRCA_meth_preproc.csv"
    mirna_preproc_csv = "BRCA/BRCA_mirna_preproc.csv"
    
    # change label from text to integer
    label_dict = {'BRCA.Normal':0, 'BRCA.Basal':1, 'BRCA.Her2':2, 'BRCA.LumA':3, 'BRCA.LumB':4}
    
    # pre-processing parameters
    mean_min_list = [0,0,0]
    mean_max_list = [float("inf"),float("inf"),float("inf")]
    var_min_list = [0.1,0.001,0.0]
    var_max_list = [float("inf"),float("inf"),float("inf")]
    retain_dim_list = [200,200,200]

    data_list, labels, patient_id, feat_name_list = load_data(mrna_csv, meth_csv, mirna_csv, meta_csv, label_dict)
    tr_idx, te_idx = load_trte_partition(trte_partition_file, patient_id)
    preproc_feat_idx = []
    for i in range(len(data_list)):
        idx = pre_processing(data_list[i][tr_idx,:], labels[tr_idx], 
                             mean_min_list[i], mean_max_list[i], var_min_list[i], var_max_list[i], retain_dim_list[i])
        preproc_feat_idx.append(idx.copy())
    
    # save preprocessed data as csv
    preproc_csv_file_list = [mrna_preproc_csv, meth_preproc_csv, mirna_preproc_csv]
    for i in range(len(data_list)):
        save_data_to_csv(preproc_csv_file_list[i], patient_id, 
                         data_list[i][:,preproc_feat_idx[i]], feat_name_list[i][preproc_feat_idx[i]])