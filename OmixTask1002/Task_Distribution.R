library(readxl)
library(VersePlot)
library(tidyverse)
library(tidytext)

file <- "tasks.xlsx"
d <- read_excel(file)
d$`Tertiary Category` <- gsub('analysis|gene|-','',d$`Tertiary Category`)
colnames(d)
word_vec <- d %>% 
  unnest_tokens(output = "token", input = `Tertiary Category`, token = "words") %>% 
  pull(token)  

word_df <- as.data.frame(table(word_vec))
word_df$word_vec <- as.character(word_df$word_vec)
p <- WordCloudPlot(word_df,word_by = 'word_vec',score_by = 'Freq')
ggsave(filename = 'WordCloud.pdf',plot = p,width = 8.1,height = 8)


d2 <- as.data.frame(table(d$`Primary Annotation`))
d2$Var1 <- as.character(d2$Var1)
d2 <- d2[order(d2$Freq),]
d2$Var1 <- factor(d2$Var1,levels = rev(d2$Var1))


library(CellScope)
ggplot(d2,aes(Freq,reorder(Var1,Freq),fill=Var1))+
  geom_bar(stat = 'identity',width = 0.7)+
  geom_text(aes(label = Freq),hjust = -0.1,size = 4)+
  theme_classic(base_line_size = 0.8)+
  scale_fill_manual(values = color.scp(d2$Var1))+
  scale_x_continuous(expand = c(0,0),limits = c(0,200))+
  labs(title = 'Distribution of OmixTask1002',x='Counts',y=NULL)+
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.5,face = 'bold',colour = 'black',size = 14),
        axis.text.x = element_text(colour = 'black',size = 11),
        axis.title.x = element_text(colour = 'black',size = 13,face = 'bold'),
        axis.text.y = element_text(colour = 'black',size = 12))
ggsave(filename = 'OmixTask1002_Distribution.pdf',width = 5,height = 6.3)

sum(d2$Freq)


