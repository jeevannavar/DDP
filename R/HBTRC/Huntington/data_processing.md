Data Processing
================

This markdown document illustrates the following processing steps: \*
Pre-processing \* Ranking of features using ANOVA \* Selection of the
most informative features using a greedy method \* Scaling and
normalization \* Creation of pairwise interaction features \* Ranking
and selection of the pairwise interaction features \* Scaling and
normalization of the interaction features

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

# AD or Control

## Loading data

``` r
cerebellum <- read_delim("../../../../HBTRC_data/GN326_MeanDataAnnotated_rev081815.txt", 
                         "\t", escape_double = FALSE, trim_ws = TRUE, skip = 32)
visualCortex <- read_delim("../../../../HBTRC_data/GN327_MeanDataAnnotated_rev081815.txt",
                           "\t", escape_double = FALSE, trim_ws = TRUE, skip = 32)
prefrontalCortex <- read_delim("../../../../HBTRC_data/GN328_MeanDataAnnotated_rev081815.txt",
                               "\t", escape_double = FALSE, trim_ws = TRUE, skip = 32)

labels <- read_csv("labels_HD.csv")

# Getting the training indices
trte <- file("trte_partition_HD.txt", open = "r")
lines <- readLines(trte)
close(trte)
train_idx <- unlist(strsplit(lines[2], ","))
test_idx <- unlist(strsplit(lines[4], ","))
```

## Pre-processing

### Cerebellum

``` r
cerebellum <- cerebellum %>%
  dplyr::select("Gene Symbol" | "Gene Id" | starts_with("HB")) %>%
  dplyr::rename(Gene_Symbol = "Gene Symbol", GeneID = "Gene Id") %>%
  filter(GeneID != "None") %>%
  mutate(Gene = paste(Gene_Symbol, GeneID, sep = "|")) %>%
  dplyr::select(!GeneID & !Gene_Symbol) %>%
  group_by(Gene) %>%
  summarize(across(everything(), list(mean))) %>%
  remove_rownames() %>%
  column_to_rownames("Gene")

# Keep only finite features
cerebellum <- cerebellum[is.finite(rowSums(cerebellum)),]

# Transposing to get sample ids as row names and genes as column names
cerebellum <- transpose_df(cerebellum)

# 
pattern <- "^(\\w{2}\\_\\w{3}\\_\\w{1,2})\\_1$"
cerebellum <- cerebellum %>%  
  mutate(rowname = sub(pattern, "\\1", rowname)) %>%
  dplyr::rename(patient_id = rowname)

cerebellum <- cerebellum[match(labels$patient_id, cerebellum$patient_id),]
```

### Primary Visual Cortex

``` r
visualCortex <- visualCortex %>%
  dplyr::select("Gene Symbol" | "Gene Id" | starts_with("HB")) %>%
  dplyr::rename(Gene_Symbol = "Gene Symbol", GeneID = "Gene Id") %>%
  filter(GeneID != "None") %>%
  mutate(Gene = paste(Gene_Symbol, GeneID, sep = "|")) %>%
  dplyr::select(!GeneID & !Gene_Symbol) %>%
  group_by(Gene) %>%
  summarize(across(everything(), list(mean))) %>%
  remove_rownames() %>%
  column_to_rownames("Gene")

# Keep only finite features
visualCortex <- visualCortex[is.finite(rowSums(visualCortex)),]

# Transposing to get sample ids as row names and genes as column names
visualCortex <- transpose_df(visualCortex)

# 
pattern <- "^(\\w{2}\\_\\w{3}\\_\\w{1,2})\\_1$"
visualCortex <- visualCortex %>%  
  mutate(rowname = sub(pattern, "\\1", rowname)) %>%
  dplyr::rename(patient_id = rowname)

visualCortex <- visualCortex[match(labels$patient_id, visualCortex$patient_id),]
```

### Pre-frontal Cortex

``` r
prefrontalCortex <- prefrontalCortex %>%
  dplyr::select("Gene Symbol" | "Gene Id" | starts_with("HB")) %>%
  dplyr::rename(Gene_Symbol = "Gene Symbol", GeneID = "Gene Id") %>%
  filter(GeneID != "None") %>%
  mutate(Gene = paste(Gene_Symbol, GeneID, sep = "|")) %>%
  dplyr::select(!GeneID & !Gene_Symbol) %>%
  group_by(Gene) %>%
  summarize(across(everything(), list(mean))) %>%
  remove_rownames() %>%
  column_to_rownames("Gene")

# Keep only finite features
prefrontalCortex <- prefrontalCortex[is.finite(rowSums(prefrontalCortex)),]

# Transposing to get sample ids as row names and genes as column names
prefrontalCortex <- transpose_df(prefrontalCortex)

# 
pattern <- "^(\\w{2}\\_\\w{3}\\_\\w{1,2})\\_1$"
prefrontalCortex <- prefrontalCortex %>%  
  mutate(rowname = sub(pattern, "\\1", rowname)) %>%
  dplyr::rename(patient_id = rowname)

prefrontalCortex <- prefrontalCortex[match(labels$patient_id, prefrontalCortex$patient_id),]
```

## Ranking of fearures using ANOVA

### Cerebellum

``` r
# Selecting only training samples for doing the processing
cerebellum_train <- filter(cerebellum, cerebellum$patient_id %in% train_idx)
labels_train <- filter(labels, labels$patient_id %in% train_idx)

#Check that the labels are identical
identical(cerebellum_train$patient_id, labels_train$patient_id) 
```

    ## [1] TRUE

``` r
cerebellum_summary <- run_anova(cerebellum_train, labels_train, num_cpus=40)

cerebellum_sorted <- cerebellum_summary %>% 
                       filter(p_value<=0.01/nrow(cerebellum_summary)) %>% 
                       arrange(desc(f_value))

namelist <- c(cerebellum_sorted["gene"])
index <- match(namelist[[1]], names(cerebellum))
index <- prepend(index,1)     #To add the column with the patient_ids too to the output
cerebellum_anova <- cerebellum[,index]
cerebellum_anova_train <- cerebellum_train[,index]

#write.csv(cerebellum_anova, "cerebellum_anova.csv", row.names = FALSE)
```

### Primary Visual Cortex

``` r
# Selecting only training samples for doing the processing
visualCortex_train <- filter(visualCortex, visualCortex$patient_id %in% train_idx)
labels_train <- filter(labels, labels$patient_id %in% train_idx)

#Check that the labels are identical
identical(visualCortex_train$patient_id, labels_train$patient_id) 
```

    ## [1] TRUE

``` r
visualCortex_summary <- run_anova(visualCortex_train, labels_train, num_cpus=40)

visualCortex_sorted <- visualCortex_summary %>% 
                       filter(p_value<=0.01/nrow(visualCortex_summary)) %>% 
                       arrange(desc(f_value))

namelist <- c(visualCortex_sorted["gene"])
index <- match(namelist[[1]], names(visualCortex))
index <- prepend(index,1)     #To add the column with the patient_ids too to the output
visualCortex_anova <- visualCortex[,index]
visualCortex_anova_train <- visualCortex_train[,index]

#write.csv(visualCortex_anova, "visualCortex_anova.csv", row.names = FALSE)
```

### Pre-frontal Cortex

``` r
# Selecting only training samples for doing the processing
prefrontalCortex_train <- filter(prefrontalCortex, prefrontalCortex$patient_id %in% train_idx)
labels_train <- filter(labels, labels$patient_id %in% train_idx)

#Check that the labels are identical
identical(prefrontalCortex_train$patient_id, labels_train$patient_id) 
```

    ## [1] TRUE

``` r
prefrontalCortex_summary <- run_anova(prefrontalCortex_train, labels_train, num_cpus=40)

prefrontalCortex_sorted <- prefrontalCortex_summary %>% 
                       filter(p_value<=0.01/nrow(prefrontalCortex_summary)) %>% 
                       arrange(desc(f_value))

namelist <- c(prefrontalCortex_sorted["gene"])
index <- match(namelist[[1]], names(prefrontalCortex))
index <- prepend(index,1)     #To add the column with the patient_ids too to the output
prefrontalCortex_anova <- prefrontalCortex[,index]
prefrontalCortex_anova_train <- prefrontalCortex_train[,index]

#write.csv(prefrontalCortex_anova, "prefrontalCortex_anova.csv", row.names = FALSE)
```

## Selection of the most informative features using a greedy method

### Cerebellum

``` r
cerebellum_best <- greedySelect(cerebellum_anova_train, cerebellum_sorted["gene"][[1]], corrThreshold = 0.7)

index <- match(cerebellum_best, names(cerebellum_anova))
index <- prepend(index,1)     #To add the column with the patient_ids too to the output
cerebellum_top <- cerebellum_anova[,index]

#write.csv(cerebellum_top, "cerebellum_top.csv", row.names = FALSE)
```

### Primary Visual Cortex

``` r
visualCortex_best <- greedySelect(visualCortex_anova_train, visualCortex_sorted["gene"][[1]], corrThreshold = 0.7)

index <- match(visualCortex_best, names(visualCortex_anova))
index <- prepend(index,1)     #To add the column with the patient_ids too to the output
visualCortex_top <- visualCortex_anova[,index]

#write.csv(visualCortex_top, "visualCortex_top.csv", row.names = FALSE)
```

### Pre-frontal Cortex

``` r
prefrontalCortex_best <- greedySelect(prefrontalCortex_anova_train, prefrontalCortex_sorted["gene"][[1]], corrThreshold = 0.7)

index <- match(prefrontalCortex_best, names(prefrontalCortex_anova))
index <- prepend(index,1)     #To add the column with the patient_ids too to the output
prefrontalCortex_top <- prefrontalCortex_anova[,index]

#write.csv(prefrontalCortex_top, "prefrontalCortex_top.csv", row.names = FALSE)
```

## Scaling and Normalization (and saving to files)

### Cerebellum

``` r
data <- cerebellum_top

data_train <- filter(data, patient_id %in% train_idx)
data_test <- filter(data, patient_id %in% test_idx)

n = length(data_train)
mean <- apply(data_train[2:n], 2, mean)
std <- apply(data_train[2:n], 2, sd)

data_train[2:n] <- scale(data_train[2:n], center = mean, scale = std)
data_test[2:n] <- scale(data_test[2:n], center = mean, scale = std)

data <- rbind(data_train, data_test)
data <- data[match(labels$patient_id, data$patient_id), ]
write_csv(data, "cerebellum_HD.csv")
```

### Primary Visual Cortex

``` r
data <- visualCortex_top

data_train <- filter(data, patient_id %in% train_idx)
data_test <- filter(data, patient_id %in% test_idx)

n = length(data_train)
mean <- apply(data_train[2:n], 2, mean)
std <- apply(data_train[2:n], 2, sd)

data_train[2:n] <- scale(data_train[2:n], center = mean, scale = std)
data_test[2:n] <- scale(data_test[2:n], center = mean, scale = std)

data <- rbind(data_train, data_test)
data <- data[match(labels$patient_id, data$patient_id), ]
write_csv(data, "visualCortex_HD.csv")
```

### Pre-frontal Cortex

``` r
data <- prefrontalCortex_top

data_train <- filter(data, patient_id %in% train_idx)
data_test <- filter(data, patient_id %in% test_idx)

n = length(data_train)
mean <- apply(data_train[2:n], 2, mean)
std <- apply(data_train[2:n], 2, sd)

data_train[2:n] <- scale(data_train[2:n], center = mean, scale = std)
data_test[2:n] <- scale(data_test[2:n], center = mean, scale = std)

data <- rbind(data_train, data_test)
data <- data[match(labels$patient_id, data$patient_id), ]
write_csv(data, "prefrontalCortex_HD.csv")
```

# Braak Severity

## Loading data

``` r
labels <- read_csv("labels_vonsattel.csv")

cerebellum <- filter(cerebellum, patient_id %in% labels$patient_id)
visualCortex <- filter(visualCortex, patient_id %in% labels$patient_id)
prefrontalCortex <- filter(prefrontalCortex, patient_id %in% labels$patient_id)

# Getting the training indices
trte <- file("trte_partition_vonsattel.txt", open = "r")
lines <- readLines(trte)
close(trte)
train_idx <- unlist(strsplit(lines[2], ","))
test_idx <- unlist(strsplit(lines[4], ","))
```

## Ranking of fearures using ANOVA

### Cerebellum

``` r
# Selecting only training samples for doing the processing
cerebellum_train <- filter(cerebellum, cerebellum$patient_id %in% train_idx)
labels_train <- filter(labels, labels$patient_id %in% train_idx)

#Check that the labels are identical
identical(cerebellum_train$patient_id, labels_train$patient_id) 
```

    ## [1] TRUE

``` r
cerebellum_summary <- run_anova(cerebellum_train, labels_train, num_cpus=40)

cerebellum_sorted <- cerebellum_summary %>% 
                       filter(p_value<=0.01/nrow(cerebellum_summary)) %>% 
                       arrange(desc(f_value))

namelist <- c(cerebellum_sorted["gene"])
index <- match(namelist[[1]], names(cerebellum))
index <- prepend(index,1)     #To add the column with the patient_ids too to the output
cerebellum_anova <- cerebellum[,index]
cerebellum_anova_train <- cerebellum_train[,index]

#write.csv(cerebellum_anova, "cerebellum_anova.csv", row.names = FALSE)
```

### Primary Visual Cortex

``` r
# Selecting only training samples for doing the processing
visualCortex_train <- filter(visualCortex, visualCortex$patient_id %in% train_idx)
labels_train <- filter(labels, labels$patient_id %in% train_idx)

#Check that the labels are identical
identical(visualCortex_train$patient_id, labels_train$patient_id) 
```

    ## [1] TRUE

``` r
visualCortex_summary <- run_anova(visualCortex_train, labels_train, num_cpus=40)

visualCortex_sorted <- visualCortex_summary %>% 
                       filter(p_value<=0.01/nrow(visualCortex_summary)) %>% 
                       arrange(desc(f_value))

namelist <- c(visualCortex_sorted["gene"])
index <- match(namelist[[1]], names(visualCortex))
index <- prepend(index,1)     #To add the column with the patient_ids too to the output
visualCortex_anova <- visualCortex[,index]
visualCortex_anova_train <- visualCortex_train[,index]

#write.csv(visualCortex_anova, "visualCortex_anova.csv", row.names = FALSE)
```

### Pre-frontal Cortex

``` r
# Selecting only training samples for doing the processing
prefrontalCortex_train <- filter(prefrontalCortex, prefrontalCortex$patient_id %in% train_idx)
labels_train <- filter(labels, labels$patient_id %in% train_idx)

#Check that the labels are identical
identical(prefrontalCortex_train$patient_id, labels_train$patient_id) 
```

    ## [1] TRUE

``` r
prefrontalCortex_summary <- run_anova(prefrontalCortex_train, labels_train, num_cpus=40)

prefrontalCortex_sorted <- prefrontalCortex_summary %>% 
                       filter(p_value<=0.01/nrow(prefrontalCortex_summary)) %>% 
                       arrange(desc(f_value))

namelist <- c(prefrontalCortex_sorted["gene"])
index <- match(namelist[[1]], names(prefrontalCortex))
index <- prepend(index,1)     #To add the column with the patient_ids too to the output
prefrontalCortex_anova <- prefrontalCortex[,index]
prefrontalCortex_anova_train <- prefrontalCortex_train[,index]

#write.csv(prefrontalCortex_anova, "prefrontalCortex_anova.csv", row.names = FALSE)
```

## Selection of the most informative features using a greedy method

### Cerebellum

``` r
cerebellum_best <- greedySelect(cerebellum_anova_train, cerebellum_sorted["gene"][[1]], corrThreshold = 0.7)

index <- match(cerebellum_best, names(cerebellum_anova))
index <- prepend(index,1)     #To add the column with the patient_ids too to the output
cerebellum_top <- cerebellum_anova[,index]

#write.csv(cerebellum_top, "cerebellum_top.csv", row.names = FALSE)
```

### Primary Visual Cortex

``` r
visualCortex_best <- greedySelect(visualCortex_anova_train, visualCortex_sorted["gene"][[1]], corrThreshold = 0.7)

index <- match(visualCortex_best, names(visualCortex_anova))
index <- prepend(index,1)     #To add the column with the patient_ids too to the output
visualCortex_top <- visualCortex_anova[,index]

#write.csv(visualCortex_top, "visualCortex_top.csv", row.names = FALSE)
```

### Pre-frontal Cortex

``` r
prefrontalCortex_best <- greedySelect(prefrontalCortex_anova_train, prefrontalCortex_sorted["gene"][[1]], corrThreshold = 0.7)

index <- match(prefrontalCortex_best, names(prefrontalCortex_anova))
index <- prepend(index,1)     #To add the column with the patient_ids too to the output
prefrontalCortex_top <- prefrontalCortex_anova[,index]

#write.csv(prefrontalCortex_top, "prefrontalCortex_top.csv", row.names = FALSE)
```

## Scaling and Normalization (and saving to files)

### Cerebellum

``` r
data <- cerebellum_top

data_train <- filter(data, patient_id %in% train_idx)
data_test <- filter(data, patient_id %in% test_idx)

n = length(data_train)
mean <- apply(data_train[2:n], 2, mean)
std <- apply(data_train[2:n], 2, sd)

data_train[2:n] <- scale(data_train[2:n], center = mean, scale = std)
data_test[2:n] <- scale(data_test[2:n], center = mean, scale = std)

data <- rbind(data_train, data_test)
data <- data[match(labels$patient_id, data$patient_id), ]
write_csv(data, "cerebellum_vonsattel.csv")
```

### Primary Visual Cortex

``` r
data <- visualCortex_top

data_train <- filter(data, patient_id %in% train_idx)
data_test <- filter(data, patient_id %in% test_idx)

n = length(data_train)
mean <- apply(data_train[2:n], 2, mean)
std <- apply(data_train[2:n], 2, sd)

data_train[2:n] <- scale(data_train[2:n], center = mean, scale = std)
data_test[2:n] <- scale(data_test[2:n], center = mean, scale = std)

data <- rbind(data_train, data_test)
data <- data[match(labels$patient_id, data$patient_id), ]
write_csv(data, "visualCortex_vonsattel.csv")
```

### Pre-frontal Cortex

``` r
data <- prefrontalCortex_top

data_train <- filter(data, patient_id %in% train_idx)
data_test <- filter(data, patient_id %in% test_idx)

n = length(data_train)
mean <- apply(data_train[2:n], 2, mean)
std <- apply(data_train[2:n], 2, sd)

data_train[2:n] <- scale(data_train[2:n], center = mean, scale = std)
data_test[2:n] <- scale(data_test[2:n], center = mean, scale = std)

data <- rbind(data_train, data_test)
data <- data[match(labels$patient_id, data$patient_id), ]
write_csv(data, "prefrontalCortex_vonsattel.csv")
```

``` r
# Removing some tibbles that will not be used later
rm(cerebellum_anova, cerebellum_anova_train, cerebellum_sorted, cerebellum_summary, cerebellum_train,
   visualCortex_anova, visualCortex_anova_train, visualCortex_sorted, visualCortex_summary, visualCortex_train,
   prefrontalCortex_anova, prefrontalCortex_anova_train, prefrontalCortex_sorted, prefrontalCortex_summary, prefrontalCortex_train,
   data, data_train, data_test, mean, std)
```
