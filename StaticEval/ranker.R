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

# Omics data columns
omics_cols <- c("Epigenomics", "Genomics", "Metabolomics", "Microbiome", 
                "Multi.Omics.Integration", "Pharmacogenomics", "Proteomics", 
                "Single.cell.Omics", "Spatial.Omics", "Transcriptomics")
data2 <- data[,c(colnames(data)[2:5],omics_cols,'avg_total_score')]

omics_clean <- c("Epigenomics", "Genomics", "Metabolomics", "Microbiome", 
                 "Multi-Omics", "Pharmacogenomics", "Proteomics", 
                 "Single-cell", "Spatial Omics", "Transcriptomics")
colnames(data2)[5:14] <- omics_clean

# Load necessary packages
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(CellScope)

# =============================================================================
# Data Preprocessing
# =============================================================================

# Sort by avg_total_score in descending order
data2_sorted <- data2[order(data2$avg_total_score, decreasing = TRUE), ]

# Get omics column names (columns 5-14)
omics_cols <- colnames(data2_sorted)[5:14]

# Create dataframe to store results
polar_data <- data.frame()

# Extract top 5 for each dimension
for(i in 1:length(omics_cols)) {
  omics <- omics_cols[i]
  
  # Sort by current omics dimension, get top 5
  top5 <- data2_sorted %>%
    arrange(desc(.data[[omics]])) %>%
    slice_head(n = 5)
  
  # Calculate angle (360 degrees divided by 10)
  angle <- (i - 1) * 360 / length(omics_cols)
  
  # Create data for each rank
  for(rank in 1:5) {
    polar_data <- rbind(polar_data, data.frame(
      omics = omics,
      rank = rank,
      model = top5$Name[rank],
      score = top5[[omics]][rank],
      angle = angle,
      radius = rank  # rank 1 in innermost circle, rank 5 in outermost
    ))
  }
}

# =============================================================================
# Prepare Colors Based on Models
# =============================================================================

# Get all unique models that appear in top 5
unique_models <- unique(polar_data$model)
n_models <- length(unique_models)
model_colors <- color.scp(unique_models)

# =============================================================================
# Create Polar Coordinate Plot
# =============================================================================

p <- ggplot(polar_data, aes(x = angle, y = radius)) +
  
  # Add concentric circle grid lines
  geom_hline(yintercept = 1:5, color = "grey80", size = 0.5, linetype = "dashed") +
  
  # Add radial grid lines  
  geom_vline(xintercept = seq(0, 360, length.out = length(omics_cols) + 1)[1:length(omics_cols)], 
             color = "grey80", size = 0.3) +
  
  # Add points colored by model
  geom_point(aes(color = model), size = 3.5, alpha = 0.8) +
  
  # Convert to polar coordinates
  coord_polar(theta = "x", start = 0) +
  
  # Set colors for models
  scale_color_manual(values = model_colors, name = "Models") +
  
  # Set x-axis (angle axis)
  scale_x_continuous(breaks = seq(0, 360, length.out = length(omics_cols) + 1)[1:length(omics_cols)],
                     labels = omics_cols,
                     limits = c(0, 360)) +
  
  # Set y-axis (radial axis)  
  scale_y_continuous(breaks = 1:5,position = 'right',
                     labels = c("1st", "2nd", "3rd", "4th", "5th"),
                     limits = c(0.1, 5)) +
  # Theme settings
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.y =  element_blank(),
    axis.text.x = element_text(size = 10, color = "black", face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 10, color = "black",face = "bold"),
    legend.text = element_text(size = 9,color = "black"),
    legend.key.size = unit(0.6, "cm"),
    legend.background = element_blank(),
    plot.margin = margin(0,0.2,0,2,unit = 'cm')
  ) +
  
  # Add guides for better legend
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))

# Display plot
print(p)
ggsave(filename = 'radar-rank.pdf',width = 8,height = 7)