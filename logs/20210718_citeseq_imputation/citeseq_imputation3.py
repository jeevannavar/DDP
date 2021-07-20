# Loading necessary libraries and functions
import sys
sys.path.append('../')
import pandas as pd
import numpy as np

# File locations
rna_file = "../../CITEseq_BMNC/rna.csv"
adt_file = "../R/CITE_Seq_BMNC/adt.csv"
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
# from sklearn.model_selection import train_test_split
# bmnc_train, bmnc_val, meta_train, meta_val = train_test_split(bmnc_train, meta_train, test_size = 0.2, random_state = 42, stratify = meta_train)

# SMOTE for slightly adjusting for imbalanced classes
# A mixture of oversampling and undersampling
from imblearn.over_sampling import SMOTE
from imblearn.under_sampling import RandomUnderSampler
SEED=42
# Setting the target so that the data is undersampled so that every cell_type that has more than 200 samples is reduced to 200\n",
target = meta_train["cell_type"].value_counts().to_dict()
for key in target:
    target[key] = min(target[key], 1000)

under = RandomUnderSampler(random_state=SEED, sampling_strategy=target)
bmnc_train, meta_train = under.fit_resample(bmnc_train, meta_train["cell_type"])

# Setting the target so that cell_types are all over-sampled to max amount
for key in target:
    target[key] = 1000
smote = SMOTE(random_state=SEED, sampling_strategy=target)
bmnc_train, meta_train = smote.fit_resample(bmnc_train, meta_train)
# Largest classes still have >1k samples, but even the smallest ones have min of 500 now

print(bmnc_train.shape)
print(bmnc_test.shape)
print("\n")

# Remove test set values after saving the ground truth values 
mask = [x == "adt" for x in datatypes]
truth = bmnc_test.loc[:,mask]
bmnc_test.loc[:,mask] = np.nan

# Vertically joining tcga and metabric data; missing features have NaN values
merged_data = pd.concat([bmnc_train, bmnc_test], keys=["train", "test"])

# Imputing the missing features in metabric data set
from sklearn.impute import KNNImputer
neighbors = [7500, 10000]
for n in neighbors:
    imputer = KNNImputer(n_neighbors=n)
    imputer.fit(merged_data.loc["train"])

    bmnc_test_imputed = imputer.transform(merged_data.loc["test"])
    bmnc_test_imputed = pd.DataFrame(bmnc_test_imputed, columns = merged_data.columns, index = merged_data.loc["test"].index)

    imputed_values = bmnc_test_imputed.loc[:,mask]
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
    
    filename = "smote1000_spearman_n_" + str(n) + ".csv"
    sp_score.to_csv(filename)
    filename = "smote1000_pearson_n_" + str(n) + ".csv"
    pe_score.to_csv(filename)
