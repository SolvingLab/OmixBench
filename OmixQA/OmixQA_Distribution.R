library(readxl)
library(VersePlot)
library(tidyverse)
library(tidytext)

d <- readxl::read_xlsx('OmixQA-Task-Info.xlsx')

d2 <- as.data.frame(table(d$Omics))
d2$Var1 <- as.character(d2$Var1)
d2 <- d2[order(d2$Freq),]
d2$Var1 <- factor(d2$Var1,levels = rev(d2$Var1))


library(CellScope)
ggplot(d2,aes(Freq,reorder(Var1,Freq),fill=Var1))+
  geom_bar(stat = 'identity',width = 0.6)+
  geom_text(aes(label = Freq),hjust = -0.1,size = 4)+
  theme_classic(base_line_size = 0.5)+
  scale_fill_manual(values = color.scp(d2$Var1))+
  scale_x_continuous(expand = c(0,0),limits = c(0,85),breaks = c(0,20,40,60,80),position = 'top')+
  labs(title = NULL,x='Task Counts',y=NULL)+
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.5,colour = 'black',size = 14),
        axis.text.x = element_text(colour = 'black',size = 11),
        axis.title.x = element_text(colour = 'black',size = 12),
        axis.text.y = element_text(colour = 'black',size = 12))
ggsave(filename = 'OmixQA_Distribution.pdf',width = 3.7,height = 3.5)

library(BioEnricher)
data <- d %>%
  dplyr::group_by(Omics,Difficulty) %>%
  dplyr::summarise(count = dplyr::n()) %>%
  dplyr::mutate(
    percentage = count / sum(count),
    label_y = 1 - (cumsum(percentage) - 0.5 * percentage),
    label = paste0(round(percentage * 100, 2), "%"),
    ymax = cumsum(percentage),
    ymin = lag(ymax, default = 0),
    labelPosition = (ymax + ymin) / 2
  )

ggplot(data, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = Difficulty)) +
  geom_rect(color = "black") +
  # geom_text(x = 1.5, aes(y = labelPosition, label = label, color = k), size = 4.5, fontface = "bold") +
  scale_fill_manual(values = c("#868686","#EFC000","#0073C2")) +
  scale_color_manual(values = c("#868686","#EFC000","#0073C2")) +
  coord_polar(theta = "y") +
  xlim(c(-1, 4)) +
  theme_void()+
  facet_wrap(~Omics,nrow = 2) +
  theme(
    legend.position = 'top',
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(size = 11, colour = "black"),
    legend.key.size = unit(5,'mm'),
    legend.title = element_blank(),
    strip.text = element_text(size = 12, colour = "black"),
    strip.background = element_blank(),
    strip.clip = 'off'
  )
ggsave(filename = 'OmixQA_Distribution2.pdf',width = 7,height = 3.5)


