rm(list = ls())
library(tidyverse)
library(CellScope)
d <- readxl::read_xlsx('Eval-Cor-Human.xlsx')
d <- d[order(d$Score,decreasing = T),]
d$Model <- factor(d$Model,levels = rev(d$Model))

colors <- c("#313695","#4575B4","#74ADD1","#ABD9E9","#E0F3F8","#FFFFBF","#FEE090","#FDAE61","#F46D43","#D73027","#A50026")

ggplot(d,aes(Model,Score,fill=Score,color=Score))+
  geom_segment(aes(x = Model,xend = Model,
                   y = 0,yend = Score),
               lineend = 'round',linewidth = 2,
               alpha = 0.95,color='#e5e5e5')+
  geom_point(size=3)+
  geom_point(size=5.5,shape = 21,fill=NA,color='black')+
  labs(x = NULL,y='Correlation Coefficient',fill=NULL,
       title = 'Correlation with Human Expert Evaluation')+
  scale_fill_gradientn(colours = color.continuous[['RdYlBu']])+
  scale_color_gradientn(colours = color.continuous[['RdYlBu']])+
  scale_y_continuous(expand = c(0,0),limits = c(0,1),breaks = c(0,0.3,0.6,0.9))+
  theme_classic(base_line_size = 0.7) +
  theme(
    plot.title = element_text(size = 14, face = "bold",hjust = 0.5),
    axis.title.x = element_text(face = "bold"),
    axis.text.x = element_text(size = 12, color = "black",angle = 45,hjust = 1),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11, color = "black"),
    legend.title = element_text(size = 12, color = "black",face = 'bold'),
    legend.position = 'none',
    legend.key.height = unit(0.5,'cm'),
    legend.key.width = unit(0.3,'cm'),
    legend.justification = c(0,1),
    legend.background = element_blank(),
    plot.margin = margin(0,0,0,25)
  )
ggsave(filename = 'CorWithHuman.pdf',width = 13,height = 3.7)

