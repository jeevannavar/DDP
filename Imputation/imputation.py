import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from copy import deepcopy
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import SimpleImputer, IterativeImputer, KNNImputer
from sklearn.linear_model import ElasticNet
from sklearn.neighbors import KNeighborsRegressor
from sklearn.ensemble import RandomForestRegressor
from metrics import *


NUMBER = 42

mrna = pd.read_csv("../R/TCGA BRCA/mrna_top1000.csv", index_col=0)
meth = pd.read_csv("../R/TCGA BRCA/meth_top1000.csv", index_col=0)
mirna = pd.read_csv("../R/TCGA BRCA/mirna_anova.csv", index_col=0)
labels = pd.read_csv("../R/TCGA BRCA/PAM50_subtype.csv", index_col=0)

all_data = pd.merge(pd.merge(mrna, meth, left_index=True, right_index=True), mirna,  left_index=True, right_index=True)
all_data = (all_data - all_data.mean())/all_data.std()
datatypes = ["mrna"]*mrna.shape[1] + ["meth"]*meth.shape[1] + ["mirna"]*mirna.shape[1]

X_train, X_test, y_train, y_test = train_test_split(all_data, labels, test_size = 0.2, random_state = NUMBER, stratify = labels)
X_train, X_val, y_train, y_val = train_test_split(X_train, y_train, test_size = 0.25, random_state = NUMBER, stratify = y_train)

X_test_truth = deepcopy(X_test)
X_val_truth = deepcopy(X_val)

methods = ["mean", "median", "enet", "knn_iter", "rf3", "rf7", "knn50", "knn75"]
imputed = {}
imputers = {
    "mean": SimpleImputer(missing_values=np.nan, strategy="mean"),
    "median": SimpleImputer(missing_values=np.nan, strategy="median"),
    "enet": IterativeImputer(estimator = ElasticNet(l1_ratio=0.75), initial_strategy = "mean", imputation_order = "random", random_state = NUMBER, n_nearest_features = 75),
    "knn_iter": IterativeImputer(estimator = KNeighborsRegressor(n_neighbors=75), initial_strategy="mean", imputation_order="random", random_state = NUMBER, n_nearest_features = 75),
    "rf": IterativeImputer(estimator = RandomForestRegressor(max_depth=10), initial_strategy="mean", imputation_order="random", random_state = NUMBER, n_nearest_features = 200),
    "knn50": KNNImputer(n_neighbors=50),
    "knn75": KNNImputer(n_neighbors=75),
    }

missing_types = [["mirna"], ["meth"], ["mrna"], ["mirna", "meth"], ["meth", "mrna"], ["mrna", "mirna"]]
df_rows = []
for missing in missing_types:
    print("\n", "="*50)
    print("Missing datatype = ", missing, "\n")
    
    mask = [x in missing for x in datatypes]
    X_test = deepcopy(X_test_truth)
    X_test.loc[:,mask] = np.nan
    truth = X_test_truth.loc[:, mask].to_numpy()
    random = np.random.rand(truth.shape[0], truth.shape[1])

    nrmse_dict = {}
    for method in methods:
        print(method)
        imp = imputers[method]
        imp.fit(X_train)
        imputed[method] = imp.transform(X_test)

        print("NRMSE = ", nrmse(truth, imputed[method][:, mask]))
        nrmse_dict[method] = nrmse(truth, imputed[method][:, mask])
        print("\n", "-"*50)
    df_rows.append(nrmse_dict)

print(df_rows)
df = pd.DataFrame(df_rows, index = ["mirna", "meth", "mrna", "mirna+meth", "meth+mrna", "mrna+mirna"])
df.to_csv("imputation_nrmse.csv")