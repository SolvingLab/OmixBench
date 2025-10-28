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

representative_models <- data %>%
  mutate(
    param_category = case_when(
      parameter_size < 2 ~ "Small (<2B)",
      parameter_size >= 2 & parameter_size < 10 ~ "Medium (2-10B)", 
      parameter_size >= 10 & parameter_size < 50 ~ "Large (10-50B)",
      parameter_size >= 50 & parameter_size < 100 ~ "Very Large (50-200B)",
      parameter_size >= 100 ~ "Massive (>200B)"
    )
  ) %>%
  group_by(param_category) %>%
  slice_max(avg_total_score, n = 4) %>%
  ungroup() %>%
  distinct(Name, .keep_all = TRUE)

library(ggrepel)

ggplot(data[complete.cases(data$parameter_size),], 
       aes(x = parameter_size, y = avg_total_score, color = Type)) +
  geom_point(size = 3) +
  geom_text_repel(
    data = representative_models,
    aes(label = Name),
    size = 3,
    color='black'
  ) +
  scale_x_log10(
    labels = function(x) paste0(x, "B"),
    breaks = c(0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)
  ) +
  scale_color_manual(
    values = c(
      "General" = "#FF6666",
      "Coder" = "#ffc857", 
      "Reasoning" = "#119da4"
    )
  ) +
  labs(
    title = "LLM Parameter Size vs Performance",
    x = "LLM Parameter Size", 
    y = "Overall Average Score",
    color = "Model Type"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12, color = "black", face = 'bold'),
    legend.text = element_text(size = 11, color = "black"),
    legend.title = element_text(size = 12, color = "black", face = 'bold'),
    legend.position = c(1, 0.03),
    legend.background = element_blank(),
    legend.box = element_blank(),
    legend.key = element_blank(),
    legend.justification = c(1, 0),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(),
    panel.grid.major.y = element_line()
  ) +
  annotation_logticks(sides = "b", alpha = 0.5)

ggsave(filename = "llm_parameter_vs_performance.pdf", width = 8, height = 5)

data2 <- data[complete.cases(data$parameter_size),]
data3 <- data2[data2$Type=='Coder',]
data4 <- data2[data2$Type!='Coder',]

cor(data2$parameter_size, data2$avg_total_score, method = 'sp')
cor(data3$parameter_size, data3$avg_total_score, method = 'sp')
cor(data4$parameter_size, data4$avg_total_score, method = 'sp')
