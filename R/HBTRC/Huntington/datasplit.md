Data Split and Label Recovery
================

In this notebook, I will get the labels from the data and then format
them as needed and also split the data into training and testing
splits.

## Loading data

``` r
data <- read_delim("../../../../HBTRC_data/GN326_MeanDataAnnotated_rev081815.txt", "\t", escape_double = FALSE, trim_ws = TRUE, skip = 32)
```

``` r
# Getting the Gene Symbol and description and storing them as an annotation file
annotation <- data %>%
  select("Gene Symbol", "Gene Id", Description) %>%
  rename(Gene_Symbol = "Gene Symbol", GeneID = "Gene Id") %>%
  filter(GeneID != "None") %>%
  distinct(Gene_Symbol, .keep_all = TRUE)

write_csv(annotation, "gene_description.csv")
```

## Labels

We will load the header row of each tissue type to get the sample ids
present. Then we will find the common samples among them and use only
those for the
analysis.

``` r
tissue1 <- read_delim("../../../../HBTRC_data/GN326_MeanDataAnnotated_rev081815.txt", "\t", escape_double = FALSE, trim_ws = TRUE, skip = 32, n_max = 0) %>% select(starts_with("HB")) %>% names()
tissue2 <- read_delim("../../../../HBTRC_data/GN327_MeanDataAnnotated_rev081815.txt", "\t", escape_double = FALSE, trim_ws = TRUE, skip = 32, n_max = 0) %>% select(starts_with("HB")) %>% names()
tissue3 <- read_delim("../../../../HBTRC_data/GN328_MeanDataAnnotated_rev081815.txt", "\t", escape_double = FALSE, trim_ws = TRUE, skip = 32, n_max = 0) %>% select(starts_with("HB")) %>% names()

common_samples <- intersect(intersect(tissue1, tissue2), tissue3)
full_set <- union(union(tissue1, tissue2), tissue3)
```

Below is code to get the disease labels, i.e., whether the sample
belongs to one among the class of Normal, Huntington’s, or Alzheimer’s.

``` r
#Below is the format for getting label infromation from the sample id
pattern <- "^\\w{2}\\_\\w{3}\\_(\\w*)$"

#Now to get samples' labels using that format
labels <- sub(pattern, "\\1", common_samples)
labels <- recode(labels, "N"="Normal", "AD"="Alzheimer's", "HD"="Huntington's")

# label_df <- tibble(Samples=common_samples, Labels=labels)
```

We, however, in this bit want more relevant information with respect to
just Huntington’s. So, I will load in some meta-data with restricted
access. It is still not personally recognizable, but it cannot be made
public. No attempt at personally identifying has been made.

``` r
meta <- read_csv("../../../../HBTRC_data/HBTRC_trait_data5_TID.txt")
matcher <- read_csv("../../../../HBTRC_data/patient_id_to_biosample_id.csv")

meta_data <- left_join(matcher, meta, by = c("Biosample Name")) %>%
                 filter(patient_id %in% common_samples)


table(meta_data$`HD case or control (combined)`)
```

    ## 
    ## affected  control 
    ##       85      299

``` r
table(meta_data$`VonSattel severity (from pathology reports)`)
```

    ## 
    ##   0   2   3   4 
    ## 299  18  41  26

``` r
labels_HD <- tibble(patient_id = meta_data$patient_id,
                    label = meta_data$`HD case or control (combined)`) %>%
                  filter(!is.na(label))

labels_vonsattel <- tibble(patient_id = meta_data$patient_id,
                           label = meta_data$`VonSattel severity (from pathology reports)`) %>%
                      filter(!is.na(label))
```

## Stratified Random Split

One for HD or not. One for VonSattel severity.

``` r
# HD or Control
train.index <- createDataPartition(labels_HD$label, p = .8, list = FALSE)
train_HD <- labels_HD[ train.index,]
test_HD  <- labels_HD[-train.index,]

# VonSattel Severity
train.index <- createDataPartition(labels_vonsattel$label, p = .8, list = FALSE)
train_vonsattel <- labels_vonsattel[ train.index,]
test_vonsattel  <- labels_vonsattel[-train.index,]
```

## Writing files

### The disease\_class file

``` r
# HD or not
labels_HD %>%
  rename(class = "label") %>%
  write.csv("labels_HD.csv", row.names = FALSE)

# VonSattel Severity
labels_vonsattel %>%
  rename(class = "label") %>%
  write.csv("labels_vonsattel.csv", row.names = FALSE)
```

### The train-test partition file

``` r
# HD or not
fileConn<-file("trte_partition_HD.txt")
writeLines(c("training patient_id", paste(unlist(train_HD["patient_id"]), collapse = ","), "testing patient_id", paste(unlist(test_HD["patient_id"]), collapse = ",")), fileConn)
close(fileConn)

# VonSattel Severity
fileConn<-file("trte_partition_vonsattel.txt")
writeLines(c("training patient_id", paste(unlist(train_vonsattel["patient_id"]), collapse = ","), "testing patient_id", paste(unlist(test_vonsattel["patient_id"]), collapse = ",")), fileConn)
close(fileConn)
```
