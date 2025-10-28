rm(list = ls())
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)
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

data3 <- data[,c('Name','Type',"avg_total_score")]
data3 <- data3[grepl('Claude|DeepSeek|Gemini|Grok',data3$Name),]
data3 <- data3[!grepl('Coder|1.5|3.5',data3$Name),]
data3 <- data3[data3$Name!='Claude-Opus-3',]
data3 <- data3[data3$Name!='Gemini-2.0-Flash',]
colnames(data3)
data3 <- data3[data3$Name!='Gemini-2.5-Pro',]
data3 <- data3[data3$Name!='Grok-2',]

data5 <- data3
data5 <- data5[order(data5$Name),]
rownames(data5) <- NULL
data5$Type <- ifelse(data5$Type=='General','NoThinking','Thinking')
data5$Model <- c('Claude-Opus-4','Claude-Opus-4',"Claude-Sonnet-3.7","Claude-Sonnet-3.7",
                 "Claude-Sonnet-4","Claude-Sonnet-4",'DeepSeek','DeepSeek',
                 'Gemini-2.5-Flash','Gemini-2.5-Flash','Grok-3','Grok-3')
data5$Model <- factor(data5$Model,levels = c('Grok-3','Gemini-2.5-Flash','DeepSeek',
                                             'Claude-Sonnet-3.7','Claude-Sonnet-4','Claude-Opus-4'))

library(CellScope)
ggplot(data5,aes(Model,avg_total_score,fill=Type))+
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7)+
  labs(x = NULL,y='Overall Average Score',
       title = 'Comparison Between Non-Thinking and Thinking Modes',fill='Mode Type')+
  scale_fill_manual(values = color.scp(unique(data5$Type)))+
  scale_y_continuous(expand = c(0,0),breaks = c(0,30,60,90))+
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold",hjust = 0.5),
    axis.title.x = element_text(face = "bold"),
    axis.text.x = element_text(size = 11, color = "black",angle = 45,hjust = 1),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.y = element_text(face = "bold"),
    legend.text = element_text(size = 11, color = "black"),
    legend.title = element_text(size = 12, color = "black",face = 'bold'),
    legend.position = 'top',
    legend.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(),
    panel.grid.major.y = element_line()
  )

ggsave(filename = "Thinking_Mode_vs_performance.pdf", width = 7, height = 4.7)