rm(list = ls())
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)
library(tidyverse)
library(RColorBrewer)

# Read data
data <- read.csv("eval_final.csv")
data$Name <- gsub(':','-',data$Name)
data$parameter_size <- as.numeric(gsub('B','',data$parameter_size))
data$Source <- Hmisc::capitalize(data$Source)

# Sort by score
data <- data %>%
  arrange(desc(avg_total_score)) %>%
  mutate(model_order = row_number())

data$avg_total_score=rowMeans(data[,8:17])

data2 <- data[data$Series=='Qwen',]
data2 <- data2[data2$parameter_size<50,]
colnames(data2)

data3 <- data2[,c('Name','Type',"avg_total_score")]
data3$Param <- str_extract(data3$Name, "\\d+(\\.\\d+)?B")
data3$Param <- factor(data3$Param,levels = rev(unique(data3$Param)))
data3 <- data3[order(data3$Param),]

t.test(data3$avg_total_score[data3$Type=='Coder'],data3$avg_total_score[data3$Type!='Coder'],paired = T)

ggplot(data3,aes(Param,avg_total_score,fill=Type))+
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7)+
  labs(x = 'Parameter Size of Qwen2.5 Models',y='Overall Average Score',
       title = 'Performance Comparison Between Qwen2.5 General and Coder Models',fill='Model Type')+
  scale_fill_manual(values = c("General" = "#FF6666",
                               "Coder" = "#ffc857"))+
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold",hjust = 0.5),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.text = element_text(size = 11, color = "black"),
    legend.title = element_text(size = 12, color = "black",face = 'bold'),
    legend.position = c(0,1),
    legend.background = element_blank(),
    legend.justification = c(0,1),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(),
    panel.grid.major.y = element_line()
  )

library(ggpattern)
data3$Type <- factor(data3$Type,c('General','Coder'))

ggplot(data3,aes(Param,avg_total_score))+
  ggpattern::geom_bar_pattern(
    aes(pattern=Type, fill=Type),
    pattern_density = 0.2,size=0.1,
    stat = "identity", position = position_dodge(width = 0.7), 
    width = 0.7
  )+
  # geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7)+
  labs(x = 'Parameter Size of Qwen2.5 Models',y='Overall Average Score',
       title = 'Comparison Between Qwen2.5 General and Coder Models',
       fill='Model Type',pattern = 'Model Type')+
  scale_fill_manual(values = c("General" = "#FF6666",
                               "Coder" = "#ffc857"))+
  scale_pattern_manual(values=c('circle','weave'))+
  scale_y_continuous(expand = c(0,0))+
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold",hjust = 0.5),
    axis.title.x = element_text(face = "bold"),
    axis.text.x = element_text(size = 11, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.y = element_text(face = "bold"),
    legend.text = element_text(size = 11, color = "black"),
    legend.title = element_text(size = 12, color = "black",face = 'bold'),
    legend.position = c(0,1),
    legend.background = element_blank(),
    legend.justification = c(0,1),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(),
    panel.grid.major.y = element_line()
  )

ggsave(filename = "Qwen_parameter_vs_performance.pdf", width = 7, height = 3.7)