# Loading necessary libraries and functions
import sys
sys.path.append('../')
import pandas as pd
import numpy as np

# File locations
rna_file = "../../CITEseq_BMNC/rna.csv"
adt_file = "../R/CITE_Seq_BMNC/adt.csv"
adt_xgboost_file = "../R/CITE_Seq_BMNC/adt_xgboost.csv"
adt_ctpnet_file = "../R/CITE_Seq_BMNC/adt_ctpnet.csv"
meta_file = "../R/CITE_Seq_BMNC/cell_type.csv"
trte_partition_file = "../R/CITE_Seq_BMNC/trte_partition.txt"

labels = pd.read_csv(meta_file, index_col="cell_id")
label_dict = dict(zip(set(labels["cell_type"]), range(len(set(labels["cell_type"])))))

rna = pd.read_csv(rna_file, index_col="cell_id")
adt = pd.read_csv(adt_file, index_col="cell_id")
meta = pd.read_csv(meta_file, index_col="cell_id")

# Getting training and testing indices
patient_id = rna.index.to_numpy()
with open(trte_partition_file, 'r') as f:
    lines = f.readlines()

tr_patient_id_list = lines[1].strip().split(',')
te_patient_id_list = lines[3].strip().split(',')
tr_idx = [np.where(patient_id == pid)[0][0] for pid in tr_patient_id_list]
te_idx = [np.where(patient_id == pid)[0][0] for pid in te_patient_id_list]
te_idx = sorted(te_idx)

rna_train = rna.iloc[tr_idx,:]
rna_test = rna.iloc[te_idx,:]

adt_train = adt.iloc[tr_idx,:]
adt_test = adt.iloc[te_idx,:]

meta_train = meta.iloc[tr_idx,:]
meta_test = meta.iloc[te_idx,:]

# Merging data together
bmnc_train = pd.merge(rna_train, adt_train,  left_index=True, right_index=True)
bmnc_test = pd.merge(rna_test, adt_test,  left_index=True, right_index=True)
datatypes = ["rna"]*rna_train.shape[1] + ["adt"]*adt_train.shape[1]

# Remove test set values
mask = [x == "adt" for x in datatypes]
bmnc_test.loc[:,mask] = np.nan

# Vertically joining tcga and metabric data; missing features have NaN values
merged_data = pd.concat([bmnc_train, bmnc_test], keys=["train", "test"])

# Imputing the missing features in metabric data set
from sklearn.impute import KNNImputer

imputer = KNNImputer(n_neighbors=400)
imputer.fit(merged_data.loc["train"])

bmnc_test_imputed = imputer.transform(merged_data.loc["test"])
bmnc_test_imputed = pd.DataFrame(bmnc_test_imputed, columns = merged_data.columns, index = merged_data.loc["test"].index)
bmnc_test_imputed.head()

mask = [x == "adt" for x in datatypes]

test_set = pd.merge(rna_test, adt_test,  left_index=True, right_index=True)
truth = test_set.loc[:,mask]
knn50_imp = bmnc_test_imputed.loc[:,mask]
print(truth.shape)
print(knn50_imp.shape)
print("n_neighbors=400")

sp_score = knn50_imp.corrwith(truth)
print("\n\nSpearman Scores:")
print(sp_score)
print("")
sp_score.to_csv("citeseq_sp_scores_n_neighbors400.csv")