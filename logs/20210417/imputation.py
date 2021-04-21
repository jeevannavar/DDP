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
datatypes = ["mrna"]*mrna.shape[1] + ["meth"]*meth.shape[1] + ["mirna"]*mirna.shape[1]

X_train, X_test, y_train, y_test = train_test_split(all_data, labels, test_size = 0.2, random_state = NUMBER, stratify = labels)
X_train, X_val, y_train, y_val = train_test_split(X_train, y_train, test_size = 0.25, random_state = NUMBER, stratify = y_train)

X_test_truth = deepcopy(X_test)
X_val_truth = deepcopy(X_val)

mask = [x=="mirna" for x in datatypes]
X_test.loc[:,mask] = np.nan
X_val.loc[:,mask] = np.nan

methods = ["mean", "median", "enet_0.2", "enet_0.8", "knn_iter", "rf3", "rf7", "knn1", "knn5", "knn10", "knn25", "knn50", "knn75", "knn100", "knn150", "knn200"]
imputed = {}
imputers = {
    "mean": SimpleImputer(missing_values=np.nan, strategy="mean"),
    "median": SimpleImputer(missing_values=np.nan, strategy="median"),
    "enet_0.2": IterativeImputer(estimator = ElasticNet(l1_ratio=0.2), initial_strategy = "mean", imputation_order = "random", random_state = NUMBER),
    "enet_0.8": IterativeImputer(estimator = ElasticNet(l1_ratio=0.8), initial_strategy = "mean", imputation_order = "random", random_state = NUMBER),
    "knn_iter": IterativeImputer(estimator = KNeighborsRegressor(n_neighbors=15), initial_strategy="mean", imputation_order="random", random_state = NUMBER),
    "rf3": IterativeImputer(estimator = RandomForestRegressor(max_depth=3), initial_strategy="mean", imputation_order="random", random_state = NUMBER),
    "rf7": IterativeImputer(estimator = RandomForestRegressor(max_depth=7), initial_strategy="mean", imputation_order="random", random_state = NUMBER),
    "knn1": KNNImputer(n_neighbors=1),
    "knn5": KNNImputer(n_neighbors=5),
    "knn10": KNNImputer(n_neighbors=10),
    "knn25": KNNImputer(n_neighbors=25),
    "knn50": KNNImputer(n_neighbors=50),
    "knn75": KNNImputer(n_neighbors=75),
    "knn100": KNNImputer(n_neighbors=100),
    "knn150": KNNImputer(n_neighbors=150),
    "knn200": KNNImputer(n_neighbors=200)
    }

truth = X_val_truth.loc[:, mask].to_numpy()
random = np.random.rand(truth.shape[0], truth.shape[1])

nrmse_dict = {}
for method in methods:
    print(method)
    imp = imputers[method]
    imp.fit(X_train)
    imputed[method] = imp.transform(X_val)
    
    print("MSE = ", mse(truth, imputed[method][:, mask]))
    print("RMSE = ", rmse(truth, imputed[method][:, mask]))
    print("NRMSE = ", nrmse(truth, imputed[method][:, mask]))
    nrmse_dict[method] = nrmse(truth, imputed[method][:, mask])
    
    print()
    print("-"*50)
    print()

print(nrmse_dict)