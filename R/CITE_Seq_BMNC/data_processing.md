CITE-Seq Data Processing
================

# Loading data

``` r
# Gene Expression
rna_train <- read_csv("../../../CITEseq_BMNC/Train_RNA_data.csv")
```

    ## Warning: Missing column names filled in: 'X1' [1]

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   .default = col_double(),
    ##   X1 = col_character()
    ## )
    ## ℹ Use `spec()` for the full column specifications.

``` r
rna_test <- read_csv("../../../CITEseq_BMNC/Test_RNA_data.csv")
```

    ## Warning: Missing column names filled in: 'X1' [1]

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   .default = col_double(),
    ##   X1 = col_character()
    ## )
    ## ℹ Use `spec()` for the full column specifications.

``` r
# Antibody Derived Tags (ADT)
adt_train <- read_csv("../../../CITEseq_BMNC/Train_ADT_data.csv")
```

    ## Warning: Missing column names filled in: 'X1' [1]

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   .default = col_double(),
    ##   X1 = col_character()
    ## )
    ## ℹ Use `spec()` for the full column specifications.

``` r
adt_test <- read_csv("../../../CITEseq_BMNC/Test_ADT_data.csv")
```

    ## Warning: Missing column names filled in: 'X1' [1]

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   .default = col_double(),
    ##   X1 = col_character()
    ## )
    ## ℹ Use `spec()` for the full column specifications.

``` r
# ADT predicted using XGBoost
xgboost_train <- read_csv("../../../CITEseq_BMNC/XGBoost_train_Predictions.csv")
```

    ## Warning: Missing column names filled in: 'X1' [1]

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   .default = col_double(),
    ##   X1 = col_character()
    ## )
    ## ℹ Use `spec()` for the full column specifications.

``` r
xgboost_test <- read_csv("../../../CITEseq_BMNC/XGBoost_test_Predictions.csv")
```

    ## Warning: Missing column names filled in: 'X1' [1]

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   .default = col_double(),
    ##   X1 = col_character()
    ## )
    ## ℹ Use `spec()` for the full column specifications.

``` r
# ADT predicted using cTPnet
ctpnet_train <- read_csv("../../../CITEseq_BMNC/Train_ADT_cTPnet_data.csv")
```

    ## Warning: Missing column names filled in: 'X1' [1]

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   .default = col_double(),
    ##   X1 = col_character()
    ## )
    ## ℹ Use `spec()` for the full column specifications.

``` r
ctpnet_test <- read_csv("../../../CITEseq_BMNC/Test_ADT_cTPnet_data.csv")
```

    ## Warning: Missing column names filled in: 'X1' [1]

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   .default = col_double(),
    ##   X1 = col_character()
    ## )
    ## ℹ Use `spec()` for the full column specifications.

``` r
# Labels
labels_train <- read_csv("../../../CITEseq_BMNC/Train_cell_idents.csv")
```

    ## Warning: Missing column names filled in: 'X1' [1]

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   X1 = col_character(),
    ##   x = col_character()
    ## )

``` r
labels_test <- read_csv("../../../CITEseq_BMNC/Test_cell_idents.csv")
```

    ## Warning: Missing column names filled in: 'X1' [1]

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   X1 = col_character(),
    ##   x = col_character()
    ## )

# Writing labels and train-test split to file

``` r
labels <- rbind(labels_train, labels_test) %>%
  dplyr::rename(cell_id = X1, cell_type = x)

write_csv(labels, "cell_type.csv")

fileConn<-file("trte_partition.txt")
  writeLines(c("training cell_id", paste(unlist(labels_train["X1"]), collapse = ","), "testing cell_id", paste(unlist(labels_test["X1"]), collapse = ",")), fileConn)
close(fileConn)
```

# Normalising data

``` r
# Gene Expression
mean <- apply(rna_train[2:ncol(rna_train)], 2, mean)
std <- apply(rna_train[2:ncol(rna_train)], 2, sd)

rna_train[2:ncol(rna_train)] <- scale(rna_train[2:ncol(rna_train)], center = mean, scale = std)
rna_test[2:ncol(rna_test)] <- scale(rna_test[2:ncol(rna_test)], center = mean, scale = std)

# Antibody Derived Tags (ADT)
mean <- apply(adt_train[2:ncol(adt_train)], 2, mean)
std <- apply(adt_train[2:ncol(adt_train)], 2, sd)

adt_train[2:ncol(adt_train)] <- scale(adt_train[2:ncol(adt_train)], center = mean, scale = std)
adt_test[2:ncol(adt_test)] <- scale(adt_test[2:ncol(adt_test)], center = mean, scale = std)

# ADT predicted using XGBoost
mean <- apply(xgboost_train[2:ncol(xgboost_train)], 2, mean)
std <- apply(xgboost_train[2:ncol(xgboost_train)], 2, sd)

xgboost_train[2:ncol(xgboost_train)] <- scale(xgboost_train[2:ncol(xgboost_train)], center = mean, scale = std)
xgboost_test[2:ncol(xgboost_test)] <- scale(xgboost_test[2:ncol(xgboost_test)], center = mean, scale = std)

# ADT predicted using cTPnet
mean <- apply(ctpnet_train[2:ncol(ctpnet_train)], 2, mean)
std <- apply(ctpnet_train[2:ncol(ctpnet_train)], 2, sd)

ctpnet_train[2:ncol(ctpnet_train)] <- scale(ctpnet_train[2:ncol(ctpnet_train)], center = mean, scale = std)
ctpnet_test[2:ncol(ctpnet_test)] <- scale(ctpnet_test[2:ncol(ctpnet_test)], center = mean, scale = std)
```

# Combining train and test data

``` r
# Gene Expression
rna <- rbind(rna_train, rna_test) %>%
  dplyr::rename(cell_id = X1)
rna$cell_id <- str_replace(rna$cell_id, "\\.", "\\-")
rna <- rna[match(labels$cell_id, rna$cell_id),]

# Antibody Derived Tags (ADT)
adt <- rbind(adt_train, adt_test) %>%
  dplyr::rename(cell_id = X1)
adt <- adt[match(labels$cell_id, adt$cell_id),]

# ADT predicted using XGBoost
xgboost <- rbind(xgboost_train, xgboost_test) %>%
  dplyr::rename(cell_id = X1)
xgboost$cell_id <- str_replace(xgboost$cell_id, "\\.", "\\-")
xgboost <- xgboost[match(labels$cell_id, xgboost$cell_id),]

# ADT predicted using cTPnet
ctpnet <- rbind(ctpnet_train, ctpnet_test) %>%
  dplyr::rename(cell_id = X1)
ctpnet <- ctpnet[match(labels$cell_id, ctpnet$cell_id),]
```

# Writing data to files

``` r
# Gene Expression
write_csv(rna, "../../../CITEseq_BMNC/rna.csv")

# Antibody Derived Tags (ADT)
write_csv(adt, "adt.csv")

# ADT predicted using XGBoost
write_csv(xgboost, "adt_xgboost.csv")

# ADT predicted using cTPnet
write_csv(ctpnet, "adt_ctpnet.csv")
```
