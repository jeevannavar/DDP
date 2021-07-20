library(tidyverse)
library(ggrepel)

scores <- read_csv("bmnc_spearman_testset.csv")

scores %>%
  select("ADT", "smote500_3k", "XGBOOST", "CTPNET") %>%
  ggplot(aes(x=CTPNET, y=smote500_3k, label=ADT))+
  geom_point()+
  geom_text_repel(max.overlaps = 20)+
  xlim(0,1)+
  ylim(0,1)+
  geom_abline(slope=1)+
  labs(y = "KNN (N=3000)", 
       title = "Spearman Correlation of ADTs imputed using different techniques")
