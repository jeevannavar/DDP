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

# Train Validation split
from sklearn.model_selection import train_test_split
bmnc_train, bmnc_val, meta_train, meta_val = train_test_split(bmnc_train, meta_train, test_size = 0.2, random_state = 42, stratify = meta_train)

# SMOTE for slightly adjusting for imbalanced classes

from imblearn.over_sampling import SMOTE
SEED=42
target = meta_train["cell_type"].value_counts().to_dict()
for key in target:
    target[key] = max(target[key], 400)
smote = SMOTE(random_state=SEED, sampling_strategy=target)
bmnc_train, meta_train = smote.fit_resample(bmnc_train, meta_train["cell_type"])
# Largest classes still have >1k samples, but even the smallest ones have min of 400 now


# Remove validation set values after saving the ground truth values 
mask = [x == "adt" for x in datatypes]
truth = bmnc_val.loc[:,mask]
bmnc_val.loc[:,mask] = np.nan

# Vertically joining tcga and metabric data; missing features have NaN values
merged_data = pd.concat([bmnc_train, bmnc_val], keys=["train", "val"])

# Imputing the missing features in metabric data set
from sklearn.impute import KNNImputer
neighbors = [2500, 3000, 4000, 5000] #[1300, 1500, 1750, 2000] #[1, 100, 300, 500, 700, 900, 1100]
for n in neighbors:
    imputer = KNNImputer(n_neighbors=n)
    imputer.fit(merged_data.loc["train"])

    bmnc_val_imputed = imputer.transform(merged_data.loc["val"])
    bmnc_val_imputed = pd.DataFrame(bmnc_val_imputed, columns = merged_data.columns, index = merged_data.loc["val"].index)

    imputed_values = bmnc_val_imputed.loc[:,mask]
    print("n_neighbors=", n)

    sp_score = imputed_values.corrwith(truth, method="spearman")
    pe_score = imputed_values.corrwith(truth, method="pearson")
    print("Spearman Scores:")
    print(sp_score)
    print("\n\n")
    print("Pearson Scores:")
    print(pe_score)
    print("\n\n")
    print("-"*50)
    
    filename = "citeseq_spearman_scores_n_" + str(n) + ".csv"
    sp_score.to_csv(filename)
    filename = "citeseq_pearson_scores_n_" + str(n) + ".csv"
    pe_score.to_csv(filename)
