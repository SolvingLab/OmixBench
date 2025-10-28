rm(list = ls()) 
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)
library(RColorBrewer)

# 读取数据
data <- read.csv("eval_final.csv")
data$Name <- gsub(':','-',data$Name)
data$Name
data$parameter_size <- as.numeric(gsub('B','',data$parameter_size))
data$Source <- Hmisc::capitalize(data$Source)

data <- data %>%
  arrange(desc(avg_total_score)) %>%
  mutate(model_order = row_number())
data$avg_total_score=rowMeans(data[,8:17])

omics_cols <- c("Epigenomics", "Genomics", "Metabolomics", "Microbiome", 
                "Multi.Omics.Integration", "Pharmacogenomics", "Proteomics", 
                "Single.cell.Omics", "Spatial.Omics", "Transcriptomics")
data2 <- data[,c(colnames(data)[2:5],omics_cols,'avg_total_score')]

omics_clean <- c("Epigenomics", "Genomics", "Metabolomics", "Microbiome", 
                 "Multi-Omics", "Pharmacogenomics", "Proteomics", 
                 "Single-cell", "Spatial Omics", "Transcriptomics")
colnames(data2)[5:14] <- omics_clean
write.csv(data2,file = 'forheatmap.csv',row.names = F)


library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(CellScope)

# =============================================================================
# Data Processing
# =============================================================================

data2_sorted <- data2[order(data2$avg_total_score, decreasing = TRUE), ]

heatmap_matrix <- as.matrix(data2_sorted[, 5:14])
rownames(heatmap_matrix) <- data2_sorted$Name

left_anno_df <- data.frame(
  Series = data2_sorted$Series,
  Type = data2_sorted$Type,
  Source = data2_sorted$Source
)

# =============================================================================
# Color
# =============================================================================

# Series
series_colors <- color.scp(unique(data2_sorted$Series))

# Type
type_colors <- c('#119da4', '#FF6666', '#ffc857')
names(type_colors) <- unique(data2_sorted$Type)

# Source
source_colors <- c("Closed" = "#f07167", "Open" = "#0081a7")

# Omix Type
omics_names <- colnames(heatmap_matrix)
omics_colors <- color.scp(omics_names, palette = 'Accent')

# =============================================================================
# Annotation
# =============================================================================

left_anno <- HeatmapAnnotation(
  Series = left_anno_df$Series,
  Type = left_anno_df$Type,
  Source = left_anno_df$Source,
  col = list(
    Series = series_colors,
    Type = type_colors,
    Source = source_colors
  ),
  which = "row",
  simple_anno_size = unit(0.45, "cm"),
  gap = unit(0.2, 'cm'),
  annotation_name_side = "top"
)

right_anno <- HeatmapAnnotation(
  avg_score = anno_barplot(
    data2_sorted$avg_total_score,
    bar_width = 0.7,
    ylim = c(0, 100),
    gp = gpar(fill = "#00b4d8", col = NA)
  ),
  which = "row",
  annotation_name_side = "top",
  annotation_label = 'AvgScore',
  annotation_name_rot = 0,
  width = unit(2.5, "cm")
)

top_anno <- HeatmapAnnotation(
  Omics = omics_names,
  col = list(Omics = omics_colors),
  annotation_name_side = "left",
  show_legend = FALSE,
  show_annotation_name = FALSE,
  simple_anno_size = unit(0.4, "cm")
)

# =============================================================================
# main figure
# =============================================================================

ht <- Heatmap(
  heatmap_matrix,
  name = "Score",
  
  col = color.continuous[['dPBIRdBu']], 
  
  cluster_rows = FALSE,
  show_row_names = TRUE,
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 10),
  
  cluster_columns = FALSE,
  show_column_names = TRUE,
  column_names_side = 'top',
  column_names_rot = 45,
  top_annotation = top_anno,
  
  left_annotation = left_anno,
  right_annotation = right_anno,
  
  width = unit(6, "cm"),
  height = unit(44, "cm"),
  
  rect_gp = gpar(col = "white", lwd = 0.5),
  heatmap_legend_param = list(
    title = "Score",
    legend_height = unit(4, "cm"),
    legend_width = unit(6, "mm")
  )
)

# =============================================================================
# Output
# =============================================================================

draw(ht, 
     heatmap_legend_side = "right",
     annotation_legend_side = "top",
     merge_legend = TRUE)

dev.copy2pdf(file = 'heatmap.pdf', width = 10, height = 22)

png(filename = 'heatmap.png',res = 300, width = 2500, height = 6000)
draw(ht, 
     heatmap_legend_side = "right",
     annotation_legend_side = "top",
     merge_legend = TRUE)
dev.off()

