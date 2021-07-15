library(tidyverse)

logs <- tibble(data = c("Cerebellum", "Visual Cortex", "Prefrontal Cortex", "Primary tissues", "Primary tissues\n+ Inter-Tissue", "Primary tissues\n+ Intra-Tissue", "Primary\n+ All Interactions"),
               Train = c(.9548, .9428, .9687, .9732, .9754, .9776, .9754),
               Test = c(.9064, .8935, .9065, .9331, .9331, .9331, .9331))

logs %>%
  mutate(data = factor(logs$data, levels = c("Cerebellum", "Visual Cortex", "Prefrontal Cortex", "Primary tissues", "Primary tissues\n+ Inter-Tissue", "Primary tissues\n+ Intra-Tissue", "Primary\n+ All Interactions"))) %>%
  pivot_longer(c('Train', 'Test'), values_to = "Values", names_to = "Metrics") %>%
  ggplot(aes(x=data, y=Values, fill=factor(Metrics, levels = c("Train", "Test"))))+
  geom_bar(stat = "identity", position = "dodge", colour="black")+
  #geom_errorbar(aes(ymin=Values-SD, ymax=Values+SD), width=0.5, position=position_dodge(0.9))+
  labs(y="F1 Score", fill="Values", title="Brain Disease Classification of HBTRC Multi-Tissue Data")+
  theme(axis.text.x = element_text(angle = 0),
        panel.grid.minor.y = element_blank())+
  scale_y_continuous(breaks = seq(0,1,0.05))
