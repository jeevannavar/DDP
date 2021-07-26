Huntington’s - Pairwise Interactions
================

# Custom functions

``` r
# Custom function to transpose while preserving names
transpose_df <- function(df) {
  t_df <- data.table::transpose(df)
  colnames(t_df) <- rownames(df)
  rownames(t_df) <- colnames(df)
  t_df <- t_df %>%
    tibble::rownames_to_column(.data = .) %>%
    tibble::as_tibble(.)
  return(t_df)
}

# To perform ANOVA
run_anova <- function(train_data, train_labels, num_cpus=50) {
  cl <- makeCluster(num_cpus)
  doParallel::registerDoParallel(cl)
  anova_summary <- foreach(i = 2:ncol(train_data)) %dopar% {
    column <- names(train_data[i]) 
    avz <- broom:: tidy(aov(train_data[,i][[1]] ~ unlist(train_labels$class))) # Each ANOVA test iterating through each column
    return(c(gene=column, f_value=avz$statistic[[1]], p_value=avz$p.value[[1]]))
  }
  stopCluster(cl)
  anova_summary <- as.data.frame(do.call(rbind, anova_summary), stringsAsFactors = FALSE)
  anova_summary <- transform(anova_summary, p_value = as.numeric(p_value), f_value = as.numeric(f_value))
  return(anova_summary)
}

run_anova_ <- function(train_data, train_labels, num_cpus=50) {
  cl <- makeCluster(num_cpus)
  doParallel::registerDoParallel(cl)
  anova_summary <- foreach(i = 2:ncol(train_data)) %dopar% {
    column <- names(train_data[i]) 
    avz <- broom:: tidy(aov(train_data[,i] ~ unlist(train_labels$class))) # Each ANOVA test iterating through each column
    return(c(gene=column, f_value=avz$statistic[[1]], p_value=avz$p.value[[1]]))
  }
  stopCluster(cl)
  anova_summary <- as.data.frame(do.call(rbind, anova_summary), stringsAsFactors = FALSE)
  anova_summary <- transform(anova_summary, p_value = as.numeric(p_value), f_value = as.numeric(f_value))
  return(anova_summary)
}

# A greedy approach to selecting the top non-correlated features
greedySelect <- function(dframe, featureList, corrThreshold = 0.8, max_features = 1000){
  selected = c(featureList[1])
    for (i in c(2:length(featureList))) {
      tempList <- c(cor(dframe[featureList[i]], dframe[selected[1]])) 
      for (s in c(1:length(selected))){
        temp <- cor(dframe[featureList[i]], dframe[selected[s]])
        tempList <- append(tempList, temp)
        if (abs(temp) > corrThreshold) {break}
        }
      if (all( abs(tempList) < corrThreshold)) {
        selected <- append(selected, featureList[i])
      }
      if (length(selected) == max_features) {break}
    }
  return(selected)
}
```

# HD or Not

For each of the pairs, the data is generated, top features selected, the
data normalized, and then saved to files.

``` r
# Some preliminaries
# Selecting only training samples for doing the processing
cere <- read_csv("cerebellum_HD.csv")
visual <- read_csv("visualCortex_HD.csv")
prefrontal <- read_csv("prefrontalCortex_HD.csv")
labels <- read_csv("labels_HD.csv")

# Getting the training indices
trte <- file("trte_partition_HD.txt", open = "r")
lines <- readLines(trte)
close(trte)
train_idx <- unlist(strsplit(lines[2], ","))
test_idx <- unlist(strsplit(lines[4], ","))


cere_train <- filter(cere, cere$patient_id %in% train_idx)
visual_train <- filter(visual, visual$patient_id %in% train_idx)
prefrontal_train <- filter(prefrontal, prefrontal$patient_id %in% train_idx)
labels_train <- filter(labels, labels$patient_id %in% train_idx)

# Check that the labels are identical
identical(labels_train$patient_id, cere_train$patient_id) #This has to be true
```

    ## [1] TRUE

``` r
identical(labels_train$patient_id, visual_train$patient_id) #This has to be true
```

    ## [1] TRUE

``` r
identical(labels_train$patient_id, prefrontal_train$patient_id) #This has to be true
```

    ## [1] TRUE

``` r
# Setting row.names
cere_train <- column_to_rownames(cere_train, var = "patient_id")
visual_train <- column_to_rownames(visual_train, var = "patient_id")
prefrontal_train <- column_to_rownames(prefrontal_train, var = "patient_id")

cere <- column_to_rownames(cere, var = "patient_id")
visual <- column_to_rownames(visual, var = "patient_id")
prefrontal <- column_to_rownames(prefrontal, var = "patient_id")
```

## Inter-tissue interaction

### Cerebellum X Primary Visual Cortex

``` r
# Just train data
col_names <- apply(expand.grid(names(cere_train), names(visual_train)), 1, paste, collapse="_")
cere_visual_train <- setNames(as.data.frame(lapply(cere_train, `*`, visual_train)), col_names)
row.names(cere_visual_train) <- row.names(cere_train)
cere_visual_train <- rownames_to_column(cere_visual_train, var="patient_id")
cere_visual_train <- cere_visual_train[match(labels_train$patient_id, cere_visual_train$patient_id),]

#ANOVA
anova_summary <- run_anova_(cere_visual_train, labels_train, num_cpus = 50)
anova_sorted <- anova_summary %>% 
                     filter(p_value<=0.0001/nrow(anova_summary)) %>% 
                     arrange(desc(f_value))

#Only selecting the top 25000 or less to be supplied
M = min(25000, nrow(anova_sorted))
topM <- cere_visual_train[, which(names(cere_visual_train) %in% anova_sorted$gene[1:M])]
top_names <- greedySelect(topM, anova_sorted["gene"][[1]][1:M], corrThreshold = 0.7)

# Full data
col_names <- apply(expand.grid(names(cere), names(visual)), 1, paste, collapse="_")
cere_visual <- setNames(as.data.frame(lapply(cere, `*`, visual)), col_names)
row.names(cere_visual) <- row.names(cere)

index <- match(top_names, names(cere_visual))
df <- cere_visual[,index]
row.names(df) <- row.names(cere)
df <- rownames_to_column(df, var = "patient_id")

# Normalization
data_train <- filter(df, patient_id %in% train_idx)
data_test <- filter(df, patient_id %in% test_idx)

n = length(data_train)
mean <- apply(data_train[2:n], 2, mean)
std <- apply(data_train[2:n], 2, sd)

data_train[2:n] <- scale(data_train[2:n], center = mean, scale = std)
data_test[2:n] <- scale(data_test[2:n], center = mean, scale = std)

data <- rbind(data_train, data_test)
data <- data[match(labels$patient_id, data$patient_id), ]
write_csv(data, "cere_visual_HD.csv")

rm(cere_visual, cere_visual_train, df, topM)
```

### Primary Visual Cortex X Prefrontal Cortex

``` r
# Just train data
col_names <- apply(expand.grid(names(visual_train), names(prefrontal_train)), 1, paste, collapse="_")
visual_prefrontal_train <- setNames(as.data.frame(lapply(visual_train, `*`, prefrontal_train)), col_names)
row.names(visual_prefrontal_train) <- row.names(visual_train)
visual_prefrontal_train <- rownames_to_column(visual_prefrontal_train, var="patient_id")
visual_prefrontal_train <- visual_prefrontal_train[match(labels_train$patient_id, visual_prefrontal_train$patient_id),]

#ANOVA
anova_summary <- run_anova_(visual_prefrontal_train, labels_train, num_cpus = 40)
anova_sorted <- anova_summary %>% 
                     filter(p_value<=0.0001/nrow(anova_summary)) %>% 
                     arrange(desc(f_value))

#Only selecting the top 25000 or less features to be supplied
M = min(25000, nrow(anova_sorted))
topM <- visual_prefrontal_train[, which(names(visual_prefrontal_train) %in% anova_sorted$gene[1:M])]
top_names <- greedySelect(topM, anova_sorted["gene"][[1]][1:M], corrThreshold = 0.7)

# Full data
col_names <- apply(expand.grid(names(visual), names(prefrontal)), 1, paste, collapse="_")
visual_prefrontal <- setNames(as.data.frame(lapply(visual, `*`, prefrontal)), col_names)
row.names(visual_prefrontal) <- row.names(visual)

index <- match(top_names, names(visual_prefrontal))
df <- visual_prefrontal[,index]
row.names(df) <- row.names(visual)
df <- rownames_to_column(df, var = "patient_id")

# Normalization
data_train <- filter(df, patient_id %in% train_idx)
data_test <- filter(df, patient_id %in% test_idx)

n = length(data_train)
mean <- apply(data_train[2:n], 2, mean)
std <- apply(data_train[2:n], 2, sd)

data_train[2:n] <- scale(data_train[2:n], center = mean, scale = std)
data_test[2:n] <- scale(data_test[2:n], center = mean, scale = std)

data <- rbind(data_train, data_test)
data <- data[match(labels$patient_id, data$patient_id), ]
write_csv(data, "visual_prefrontal_HD.csv")

rm(visual_prefrontal, visual_prefrontal_train, df, topM)
```

### Prefrontal Cortex X Cerebellum

``` r
# Just train data
col_names <- apply(expand.grid(names(prefrontal_train), names(cere_train)), 1, paste, collapse="_")
prefrontal_cere_train <- setNames(as.data.frame(lapply(prefrontal_train, `*`, cere_train)), col_names)
row.names(prefrontal_cere_train) <- row.names(prefrontal_train)
prefrontal_cere_train <- rownames_to_column(prefrontal_cere_train, var="patient_id")
prefrontal_cere_train <- prefrontal_cere_train[match(labels_train$patient_id, prefrontal_cere_train$patient_id),]

#ANOVA
anova_summary <- run_anova_(prefrontal_cere_train, labels_train, num_cpus = 40)
anova_sorted <- anova_summary %>% 
                     filter(p_value<=0.0001/nrow(anova_summary)) %>% 
                     arrange(desc(f_value))

#Only selecting the top 25000 or less features to be supplied
M = min(25000, nrow(anova_sorted))
topM <- prefrontal_cere_train[, which(names(prefrontal_cere_train) %in% anova_sorted$gene[1:M])]
top_names <- greedySelect(topM, anova_sorted["gene"][[1]][1:M], corrThreshold = 0.7)

# Full data
col_names <- apply(expand.grid(names(prefrontal), names(cere)), 1, paste, collapse="_")
prefrontal_cere <- setNames(as.data.frame(lapply(prefrontal, `*`, cere)), col_names)
row.names(prefrontal_cere) <- row.names(prefrontal)

index <- match(top_names, names(prefrontal_cere))
df <- prefrontal_cere[,index]
row.names(df) <- row.names(prefrontal)
df <- rownames_to_column(df, var = "patient_id")

# Normalization
data_train <- filter(df, patient_id %in% train_idx)
data_test <- filter(df, patient_id %in% test_idx)

n = length(data_train)
mean <- apply(data_train[2:n], 2, mean)
std <- apply(data_train[2:n], 2, sd)

data_train[2:n] <- scale(data_train[2:n], center = mean, scale = std)
data_test[2:n] <- scale(data_test[2:n], center = mean, scale = std)

data <- rbind(data_train, data_test)
data <- data[match(labels$patient_id, data$patient_id), ]
write_csv(data, "prefrontal_cere_HD.csv")

rm(prefrontal_cere, prefrontal_cere_train, df, topM)
```

## Intra-tissue interaction

### Cerebellum X Cerebellum

``` r
# Just train data
col_names <- sapply(combn(colnames(cere_train), 2, simplify = FALSE), paste, collapse = "_")
data_train <- setNames(data.frame(combn(cere_train, 2, FUN = Reduce, f = `*`)), col_names)
row.names(data_train) <- row.names(cere_train)
data_train <- rownames_to_column(data_train, var="patient_id")
data_train <- data_train[match(labels_train$patient_id, data_train$patient_id),]

#ANOVA
anova_summary <- run_anova_(data_train, labels_train, num_cpus = 40)
anova_sorted <- anova_summary %>% 
                     filter(p_value<=0.0001/nrow(anova_summary)) %>% 
                     arrange(desc(f_value))

#Only selecting the top 25000 or less features to be supplied
M = min(25000, nrow(anova_sorted))
topM <- data_train[, which(names(data_train) %in% anova_sorted$gene[1:M])]
top_names <- greedySelect(topM, anova_sorted["gene"][[1]][1:M], corrThreshold = 0.7)

# Full data
col_names <- sapply(combn(colnames(cere), 2, simplify = FALSE), paste, collapse = "_")
data <- setNames(data.frame(combn(cere, 2, FUN = Reduce, f = `*`)), col_names)
row.names(data) <- row.names(cere)

index <- match(top_names, names(data))
df <- data[,index]
row.names(df) <- row.names(cere)
df <- rownames_to_column(df, var = "patient_id")

# Normalization
data_train <- filter(df, patient_id %in% train_idx)
data_test <- filter(df, patient_id %in% test_idx)

n = length(data_train)
mean <- apply(data_train[2:n], 2, mean)
std <- apply(data_train[2:n], 2, sd)

data_train[2:n] <- scale(data_train[2:n], center = mean, scale = std)
data_test[2:n] <- scale(data_test[2:n], center = mean, scale = std)

data <- rbind(data_train, data_test)
data <- data[match(labels$patient_id, data$patient_id), ]
write_csv(data, "cere_cere_HD.csv")

rm(data, data_train, df, topM)
```

### Primary Visual Cortex X Primary Visual Cortex

``` r
# Just train data
col_names <- sapply(combn(colnames(visual_train), 2, simplify = FALSE), paste, collapse = "_")
data_train <- setNames(data.frame(combn(visual_train, 2, FUN = Reduce, f = `*`)), col_names)
row.names(data_train) <- row.names(visual_train)
data_train <- rownames_to_column(data_train, var="patient_id")
data_train <- data_train[match(labels_train$patient_id, data_train$patient_id),]

#ANOVA
anova_summary <- run_anova_(data_train, labels_train, num_cpus = 40)
anova_sorted <- anova_summary %>% 
                     filter(p_value<=0.0001/nrow(anova_summary)) %>% 
                     arrange(desc(f_value))

#Only selecting the top 25000 or less features to be supplied
M = min(25000, nrow(anova_sorted))
topM <- data_train[, which(names(data_train) %in% anova_sorted$gene[1:M])]
top_names <- greedySelect(topM, anova_sorted["gene"][[1]][1:M], corrThreshold = 0.7)

# Full data
col_names <- sapply(combn(colnames(visual), 2, simplify = FALSE), paste, collapse = "_")
data <- setNames(data.frame(combn(visual, 2, FUN = Reduce, f = `*`)), col_names)
row.names(data) <- row.names(visual)

index <- match(top_names, names(data))
df <- data[,index]
row.names(df) <- row.names(visual)
df <- rownames_to_column(df, var = "patient_id")

# Normalization
data_train <- filter(df, patient_id %in% train_idx)
data_test <- filter(df, patient_id %in% test_idx)

n = length(data_train)
mean <- apply(data_train[2:n], 2, mean)
std <- apply(data_train[2:n], 2, sd)

data_train[2:n] <- scale(data_train[2:n], center = mean, scale = std)
data_test[2:n] <- scale(data_test[2:n], center = mean, scale = std)

data <- rbind(data_train, data_test)
data <- data[match(labels$patient_id, data$patient_id), ]
write_csv(data, "visual_visual_HD.csv")

rm(data, data_train, df, topM)
```

### Prefrontal Cortex X Prefrontal Cortex

``` r
# Just train data
col_names <- sapply(combn(colnames(prefrontal_train), 2, simplify = FALSE), paste, collapse = "_")
data_train <- setNames(data.frame(combn(prefrontal_train, 2, FUN = Reduce, f = `*`)), col_names)
row.names(data_train) <- row.names(prefrontal_train)
data_train <- rownames_to_column(data_train, var="patient_id")
data_train <- data_train[match(labels_train$patient_id, data_train$patient_id),]

#ANOVA
anova_summary <- run_anova_(data_train, labels_train, num_cpus = 40)
anova_sorted <- anova_summary %>% 
                     filter(p_value<=0.0001/nrow(anova_summary)) %>% 
                     arrange(desc(f_value))

#Only selecting the top 25000 or less features to be supplied
M = min(25000, nrow(anova_sorted))
topM <- data_train[, which(names(data_train) %in% anova_sorted$gene[1:M])]
top_names <- greedySelect(topM, anova_sorted["gene"][[1]][1:M], corrThreshold = 0.7)

# Full data
col_names <- sapply(combn(colnames(prefrontal), 2, simplify = FALSE), paste, collapse = "_")
data <- setNames(data.frame(combn(prefrontal, 2, FUN = Reduce, f = `*`)), col_names)
row.names(data) <- row.names(prefrontal)

index <- match(top_names, names(data))
df <- data[,index]
row.names(df) <- row.names(prefrontal)
df <- rownames_to_column(df, var = "patient_id")

# Normalization
data_train <- filter(df, patient_id %in% train_idx)
data_test <- filter(df, patient_id %in% test_idx)

n = length(data_train)
mean <- apply(data_train[2:n], 2, mean)
std <- apply(data_train[2:n], 2, sd)

data_train[2:n] <- scale(data_train[2:n], center = mean, scale = std)
data_test[2:n] <- scale(data_test[2:n], center = mean, scale = std)

data <- rbind(data_train, data_test)
data <- data[match(labels$patient_id, data$patient_id), ]
write_csv(data, "prefrontal_prefrontal_HD.csv")

rm(data, data_train, df, topM)
```

# Vonsattel Scores

For each of the pairs, the data is generated, top features selected, the
data normalized, and then saved to files.

``` r
# Some preliminaries
# Selecting only training samples for doing the processing
cere <- read_csv("cerebellum_vonsattel.csv")
visual <- read_csv("visualCortex_vonsattel.csv")
prefrontal <- read_csv("prefrontalCortex_vonsattel.csv")
labels <- read_csv("labels_vonsattel.csv")

# Getting the training indices
trte <- file("trte_partition_vonsattel.txt", open = "r")
lines <- readLines(trte)
close(trte)
train_idx <- unlist(strsplit(lines[2], ","))
test_idx <- unlist(strsplit(lines[4], ","))


cere_train <- filter(cere, cere$patient_id %in% train_idx)
visual_train <- filter(visual, visual$patient_id %in% train_idx)
prefrontal_train <- filter(prefrontal, prefrontal$patient_id %in% train_idx)
labels_train <- filter(labels, labels$patient_id %in% train_idx)

# Check that the labels are identical
identical(labels_train$patient_id, cere_train$patient_id) #This has to be true
```

    ## [1] TRUE

``` r
identical(labels_train$patient_id, visual_train$patient_id) #This has to be true
```

    ## [1] TRUE

``` r
identical(labels_train$patient_id, prefrontal_train$patient_id) #This has to be true
```

    ## [1] TRUE

``` r
# Setting row.names
cere_train <- column_to_rownames(cere_train, var = "patient_id")
visual_train <- column_to_rownames(visual_train, var = "patient_id")
prefrontal_train <- column_to_rownames(prefrontal_train, var = "patient_id")

cere <- column_to_rownames(cere, var = "patient_id")
visual <- column_to_rownames(visual, var = "patient_id")
prefrontal <- column_to_rownames(prefrontal, var = "patient_id")
```

## Inter-tissue interaction

### Cerebellum X Primary Visual Cortex

``` r
# Just train data
col_names <- apply(expand.grid(names(cere_train), names(visual_train)), 1, paste, collapse="_")
cere_visual_train <- setNames(as.data.frame(lapply(cere_train, `*`, visual_train)), col_names)
row.names(cere_visual_train) <- row.names(cere_train)
cere_visual_train <- rownames_to_column(cere_visual_train, var="patient_id")
cere_visual_train <- cere_visual_train[match(labels_train$patient_id, cere_visual_train$patient_id),]

#ANOVA
anova_summary <- run_anova_(cere_visual_train, labels_train, num_cpus = 50)
anova_sorted <- anova_summary %>% 
                     filter(p_value<=0.0001/nrow(anova_summary)) %>% 
                     arrange(desc(f_value))

#Only selecting the top 25000 or less to be supplied
M = min(25000, nrow(anova_sorted))
topM <- cere_visual_train[, which(names(cere_visual_train) %in% anova_sorted$gene[1:M])]
top_names <- greedySelect(topM, anova_sorted["gene"][[1]][1:M], corrThreshold = 0.7)

# Full data
col_names <- apply(expand.grid(names(cere), names(visual)), 1, paste, collapse="_")
cere_visual <- setNames(as.data.frame(lapply(cere, `*`, visual)), col_names)
row.names(cere_visual) <- row.names(cere)

index <- match(top_names, names(cere_visual))
df <- cere_visual[,index]
row.names(df) <- row.names(cere)
df <- rownames_to_column(df, var = "patient_id")

# Normalization
data_train <- filter(df, patient_id %in% train_idx)
data_test <- filter(df, patient_id %in% test_idx)

n = length(data_train)
mean <- apply(data_train[2:n], 2, mean)
std <- apply(data_train[2:n], 2, sd)

data_train[2:n] <- scale(data_train[2:n], center = mean, scale = std)
data_test[2:n] <- scale(data_test[2:n], center = mean, scale = std)

data <- rbind(data_train, data_test)
data <- data[match(labels$patient_id, data$patient_id), ]
write_csv(data, "cere_visual_vonsattel.csv")

rm(cere_visual, cere_visual_train, df, topM)
```

### Primary Visual Cortex X Prefrontal Cortex

``` r
# Just train data
col_names <- apply(expand.grid(names(visual_train), names(prefrontal_train)), 1, paste, collapse="_")
visual_prefrontal_train <- setNames(as.data.frame(lapply(visual_train, `*`, prefrontal_train)), col_names)
row.names(visual_prefrontal_train) <- row.names(visual_train)
visual_prefrontal_train <- rownames_to_column(visual_prefrontal_train, var="patient_id")
visual_prefrontal_train <- visual_prefrontal_train[match(labels_train$patient_id, visual_prefrontal_train$patient_id),]

#ANOVA
anova_summary <- run_anova_(visual_prefrontal_train, labels_train, num_cpus = 40)
anova_sorted <- anova_summary %>% 
                     filter(p_value<=0.0001/nrow(anova_summary)) %>% 
                     arrange(desc(f_value))

#Only selecting the top 25000 or less features to be supplied
M = min(25000, nrow(anova_sorted))
topM <- visual_prefrontal_train[, which(names(visual_prefrontal_train) %in% anova_sorted$gene[1:M])]
top_names <- greedySelect(topM, anova_sorted["gene"][[1]][1:M], corrThreshold = 0.7)

# Full data
col_names <- apply(expand.grid(names(visual), names(prefrontal)), 1, paste, collapse="_")
visual_prefrontal <- setNames(as.data.frame(lapply(visual, `*`, prefrontal)), col_names)
row.names(visual_prefrontal) <- row.names(visual)

index <- match(top_names, names(visual_prefrontal))
df <- visual_prefrontal[,index]
row.names(df) <- row.names(visual)
df <- rownames_to_column(df, var = "patient_id")

# Normalization
data_train <- filter(df, patient_id %in% train_idx)
data_test <- filter(df, patient_id %in% test_idx)

n = length(data_train)
mean <- apply(data_train[2:n], 2, mean)
std <- apply(data_train[2:n], 2, sd)

data_train[2:n] <- scale(data_train[2:n], center = mean, scale = std)
data_test[2:n] <- scale(data_test[2:n], center = mean, scale = std)

data <- rbind(data_train, data_test)
data <- data[match(labels$patient_id, data$patient_id), ]
write_csv(data, "visual_prefrontal_vonsattel.csv")

rm(visual_prefrontal, visual_prefrontal_train, df, topM)
```

### Prefrontal Cortex X Cerebellum

``` r
# Just train data
col_names <- apply(expand.grid(names(prefrontal_train), names(cere_train)), 1, paste, collapse="_")
prefrontal_cere_train <- setNames(as.data.frame(lapply(prefrontal_train, `*`, cere_train)), col_names)
row.names(prefrontal_cere_train) <- row.names(prefrontal_train)
prefrontal_cere_train <- rownames_to_column(prefrontal_cere_train, var="patient_id")
prefrontal_cere_train <- prefrontal_cere_train[match(labels_train$patient_id, prefrontal_cere_train$patient_id),]

#ANOVA
anova_summary <- run_anova_(prefrontal_cere_train, labels_train, num_cpus = 40)
anova_sorted <- anova_summary %>% 
                     filter(p_value<=0.0001/nrow(anova_summary)) %>% 
                     arrange(desc(f_value))

#Only selecting the top 25000 or less features to be supplied
M = min(25000, nrow(anova_sorted))
topM <- prefrontal_cere_train[, which(names(prefrontal_cere_train) %in% anova_sorted$gene[1:M])]
top_names <- greedySelect(topM, anova_sorted["gene"][[1]][1:M], corrThreshold = 0.7)

# Full data
col_names <- apply(expand.grid(names(prefrontal), names(cere)), 1, paste, collapse="_")
prefrontal_cere <- setNames(as.data.frame(lapply(prefrontal, `*`, cere)), col_names)
row.names(prefrontal_cere) <- row.names(prefrontal)

index <- match(top_names, names(prefrontal_cere))
df <- prefrontal_cere[,index]
row.names(df) <- row.names(prefrontal)
df <- rownames_to_column(df, var = "patient_id")

# Normalization
data_train <- filter(df, patient_id %in% train_idx)
data_test <- filter(df, patient_id %in% test_idx)

n = length(data_train)
mean <- apply(data_train[2:n], 2, mean)
std <- apply(data_train[2:n], 2, sd)

data_train[2:n] <- scale(data_train[2:n], center = mean, scale = std)
data_test[2:n] <- scale(data_test[2:n], center = mean, scale = std)

data <- rbind(data_train, data_test)
data <- data[match(labels$patient_id, data$patient_id), ]
write_csv(data, "prefrontal_cere_vonsattel.csv")

rm(prefrontal_cere, prefrontal_cere_train, df, topM)
```

## Intra-tissue interaction

### Cerebellum X Cerebellum

``` r
# Just train data
col_names <- sapply(combn(colnames(cere_train), 2, simplify = FALSE), paste, collapse = "_")
data_train <- setNames(data.frame(combn(cere_train, 2, FUN = Reduce, f = `*`)), col_names)
row.names(data_train) <- row.names(cere_train)
data_train <- rownames_to_column(data_train, var="patient_id")
data_train <- data_train[match(labels_train$patient_id, data_train$patient_id),]

#ANOVA
anova_summary <- run_anova_(data_train, labels_train, num_cpus = 40)
anova_sorted <- anova_summary %>% 
                     filter(p_value<=0.0001/nrow(anova_summary)) %>% 
                     arrange(desc(f_value))

#Only selecting the top 25000 or less features to be supplied
M = min(25000, nrow(anova_sorted))
topM <- data_train[, which(names(data_train) %in% anova_sorted$gene[1:M])]
top_names <- greedySelect(topM, anova_sorted["gene"][[1]][1:M], corrThreshold = 0.7)

# Full data
col_names <- sapply(combn(colnames(cere), 2, simplify = FALSE), paste, collapse = "_")
data <- setNames(data.frame(combn(cere, 2, FUN = Reduce, f = `*`)), col_names)
row.names(data) <- row.names(cere)

index <- match(top_names, names(data))
df <- data[,index]
row.names(df) <- row.names(cere)
df <- rownames_to_column(df, var = "patient_id")

# Normalization
data_train <- filter(df, patient_id %in% train_idx)
data_test <- filter(df, patient_id %in% test_idx)

n = length(data_train)
mean <- apply(data_train[2:n], 2, mean)
std <- apply(data_train[2:n], 2, sd)

data_train[2:n] <- scale(data_train[2:n], center = mean, scale = std)
data_test[2:n] <- scale(data_test[2:n], center = mean, scale = std)

data <- rbind(data_train, data_test)
data <- data[match(labels$patient_id, data$patient_id), ]
write_csv(data, "cere_cere_vonsattel.csv")

rm(data, data_train, df, topM)
```

### Primary Visual Cortex X Primary Visual Cortex

``` r
# Just train data
col_names <- sapply(combn(colnames(visual_train), 2, simplify = FALSE), paste, collapse = "_")
data_train <- setNames(data.frame(combn(visual_train, 2, FUN = Reduce, f = `*`)), col_names)
row.names(data_train) <- row.names(visual_train)
data_train <- rownames_to_column(data_train, var="patient_id")
data_train <- data_train[match(labels_train$patient_id, data_train$patient_id),]

#ANOVA
anova_summary <- run_anova_(data_train, labels_train, num_cpus = 40)
anova_sorted <- anova_summary %>% 
                     filter(p_value<=0.0001/nrow(anova_summary)) %>% 
                     arrange(desc(f_value))

#Only selecting the top 25000 or less features to be supplied
M = min(25000, nrow(anova_sorted))
topM <- data_train[, which(names(data_train) %in% anova_sorted$gene[1:M])]
top_names <- greedySelect(topM, anova_sorted["gene"][[1]][1:M], corrThreshold = 0.7)

# Full data
col_names <- sapply(combn(colnames(visual), 2, simplify = FALSE), paste, collapse = "_")
data <- setNames(data.frame(combn(visual, 2, FUN = Reduce, f = `*`)), col_names)
row.names(data) <- row.names(visual)

index <- match(top_names, names(data))
df <- data[,index]
row.names(df) <- row.names(visual)
df <- rownames_to_column(df, var = "patient_id")

# Normalization
data_train <- filter(df, patient_id %in% train_idx)
data_test <- filter(df, patient_id %in% test_idx)

n = length(data_train)
mean <- apply(data_train[2:n], 2, mean)
std <- apply(data_train[2:n], 2, sd)

data_train[2:n] <- scale(data_train[2:n], center = mean, scale = std)
data_test[2:n] <- scale(data_test[2:n], center = mean, scale = std)

data <- rbind(data_train, data_test)
data <- data[match(labels$patient_id, data$patient_id), ]
write_csv(data, "visual_visual_vonsattel.csv")

rm(data, data_train, df, topM)
```

### Prefrontal Cortex X Prefrontal Cortex

``` r
# Just train data
col_names <- sapply(combn(colnames(prefrontal_train), 2, simplify = FALSE), paste, collapse = "_")
data_train <- setNames(data.frame(combn(prefrontal_train, 2, FUN = Reduce, f = `*`)), col_names)
row.names(data_train) <- row.names(prefrontal_train)
data_train <- rownames_to_column(data_train, var="patient_id")
data_train <- data_train[match(labels_train$patient_id, data_train$patient_id),]

#ANOVA
anova_summary <- run_anova_(data_train, labels_train, num_cpus = 40)
anova_sorted <- anova_summary %>% 
                     filter(p_value<=0.0001/nrow(anova_summary)) %>% 
                     arrange(desc(f_value))

#Only selecting the top 25000 or less features to be supplied
M = min(25000, nrow(anova_sorted))
topM <- data_train[, which(names(data_train) %in% anova_sorted$gene[1:M])]
top_names <- greedySelect(topM, anova_sorted["gene"][[1]][1:M], corrThreshold = 0.7)

# Full data
col_names <- sapply(combn(colnames(prefrontal), 2, simplify = FALSE), paste, collapse = "_")
data <- setNames(data.frame(combn(prefrontal, 2, FUN = Reduce, f = `*`)), col_names)
row.names(data) <- row.names(prefrontal)

index <- match(top_names, names(data))
df <- data[,index]
row.names(df) <- row.names(prefrontal)
df <- rownames_to_column(df, var = "patient_id")

# Normalization
data_train <- filter(df, patient_id %in% train_idx)
data_test <- filter(df, patient_id %in% test_idx)

n = length(data_train)
mean <- apply(data_train[2:n], 2, mean)
std <- apply(data_train[2:n], 2, sd)

data_train[2:n] <- scale(data_train[2:n], center = mean, scale = std)
data_test[2:n] <- scale(data_test[2:n], center = mean, scale = std)

data <- rbind(data_train, data_test)
data <- data[match(labels$patient_id, data$patient_id), ]
write_csv(data, "prefrontal_prefrontal_vonsattel.csv")

rm(data, data_train, df, topM)
```
