
# Libraries
library(tidyverse)
library(hrbrthemes)
library(viridis)

aa=read.csv('box_literaute_data.csv',header = T,stringsAsFactors = F)
colnames(aa)=c("Correlation",  "Specificity_score","Rank_specificity","Null_set1","Bayesian","Markov","Null_set2")

# create a dataset
library(reshape)
data=melt(as.matrix(aa))
data=data[,-1]

colnames(data)=c("Method","value")

# make V1 an ordered factor
data$Method <- factor(data$Method, levels = unique(data$Method))

# Plot
p=data %>%
  ggplot( aes(x=Method, y=value, fill=Method)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #geom_jitter(color="black", size=1, alpha=0.9) +
  theme_ipsum() +
  theme(
    legend.position="right",
    plot.title = element_text(size=25)
  ) +
  ggtitle("Literature validation from different methods") 

p2=p +theme_classic()+ 
  theme(plot.title = element_text(size=25, face="bold", margin = margin(10, 0, 10, 0)))+
  theme(
    axis.title.x = element_text(color="black", vjust=-0.55,size = 20),
    axis.title.y = element_text(color="black" , vjust=0.55,size = 20)   
  )+
  theme(axis.text.x=element_text(color="black",size=18, vjust=0.5),
        axis.text.y=element_text(color="black",size=18, vjust=0.5))+
  labs(x='Methods',y='log(counts)')+
  theme(legend.title = element_text(size=20, face="bold"))+
  #scale_color_discrete(name="Sample Type")+
  theme(legend.position="right",legend.text=element_text(size=15))


