
library(data.table)
library(FSA)

library(threejs)
library(htmlwidgets)
library(igraph)
library(visNetwork)
library(bnlearn)
library(parallel)
library(dplyr)
options(scipen = 999)

random=data.frame(fread('no_beta_padj_scores.csv',header = T,stringsAsFactors = F,check.names = F), row.names=1)


# disease specific filtering ###########
random2=random[1:1811,]
data_2=apply(random2,1, function(x){perc(x,0.25,"lt")})
#identify pos where more than 20% of cell-types are as such
pos=which(data_2>10)
#select such columns
data_3=random2[pos,]

# pathway specific filtering ###########
random3=random[1812:3140,]
part_2=apply(random3,1, function(x){perc(x,0.25,"lt")})
#identify pos where more than 30% of cell-types are as such
pos2=which(part_2>=10)
#select such columns
part_3=random3[pos2,]
dispath=rbind(data_3,part_3)

#############@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#######################################


#pathway selection for dbt based on filtered scores matrix

#paths to chosse from
beta2=as.matrix(rownames(dispath)[252:475])

file= read.csv('comb_corr_sager.csv',header = T,row.names = 1,stringsAsFactors = F,check.names = F)

######absolute correlation


file_sub=as.matrix(file[,2])
rownames(file_sub)=rownames(file)

sub_cor=merge(beta2,file_sub,by.x="V1",by.y=0,all.x=F,all.y=F)
rownames(sub_cor)=sub_cor$V1
sub_cor$V1=NULL

sub_cor2=abs(sub_cor)


set_cor=subset(sub_cor2,sub_cor2$V1.y>0.3)

#set_cor2=acast(set_cor,set_cor$X1~set_cor$X2)

###############   pval based
adult_pval=read.table('comb_pval_sager.csv',header = T,sep = ",",check.names = F,row.names = 1,stringsAsFactors = F,)
file_sub2=as.matrix(adult_pval[,2])
rownames(file_sub2)=rownames(adult_pval)

sub_pval=merge(beta2,file_sub2,by.x="V1",by.y=0,all.x=F,all.y=F)


rownames(sub_pval)=sub_pval$V1
sub_pval$V1=NULL

sub_pval2=sub_pval

set_pval=subset(sub_pval2,sub_pval2$V1.y<0.05)


###############   pval based
adult_rank=read.table('comb_rank_sager.csv',header = T,sep = ",",check.names = F,row.names = 1,stringsAsFactors = F,)
file_sub3=as.matrix(adult_rank[,2])
rownames(file_sub3)=rownames(adult_rank)

sub_rank=merge(beta2,file_sub3,by.x="V1",by.y=0,all.x=F,all.y=F)

rownames(sub_rank)=sub_rank$V1
sub_rank$V1=NULL

sub_rank2=sub_rank

set_rank=subset(sub_rank2,sub_rank2$V1.y<0.25)
###################################################     concat    ##############################
# o<-list(set_cor,set_pval,set_rank)
# net<-Reduce(function(x, y) merge(x=x, y=y, by=0, all.x=F, all.y=F), o)

test=merge(set_cor,set_pval,by=0,all.x=F,all.y=F)
test2=merge(test,set_rank,by.x = 'Row.names',by.y=0,all.x=F,all.y=F)

write.table(test2,'nor_pathset.txt',row.names = F,sep='\t',quote=F)


##extracting scores
#random=data.frame(fread('un_a_scores.csv',header = T,stringsAsFactors = F,check.names = F), row.names=1)
#####33   use filtered score file as refernece for cross check


pathset= read.table('nor_pathset.txt',header = T,sep = "\t",check.names = F,stringsAsFactors = F,)
vec=c('Diabetes_mellitus',pathset$Row.names)

names.use <- rownames(random)[(rownames(random) %in% vec)]
adult_final <- random[ names.use,]

write.csv(t(adult_final),file='nor_pathset_score.csv',quote = F)
############################################################
############################################################


####################################    NETWORK   ###################################################
library(bnlearn)
library(parallel)
library(dplyr)
#library(bnviewer)

########################################  NORMAL    ######################################################
data<-read.csv('nor_pathset_score.csv',header = T,stringsAsFactors = F,check.names = F,row.names = 1)

temp<-as.data.frame(data)
data2<-apply(temp, 2, function(x) ntile(x,4))
data3<-apply(data2, 2, as.numeric)

cl = makeCluster(2)
bay_main2 = boot.strength(as.data.frame(data3), R = 1001,algorithm = "hc",cluster = cl)

avg_main2 = averaged.network(bay_main2, threshold = 0.3)

stopCluster(cl)


library(XMRF)
library(dplyr)
library(igraph)

data<-read.csv('nor_pathset_score.csv',header = T,stringsAsFactors = F,check.names = F,row.names = 1)

data2<-apply(data, 2, as.numeric)
data2<-apply(as.data.frame(data2), 2, function(x) ntile(x,4))
data3<-t(data2)

data2_fit <- XMRF(as.matrix(data3),method="GGM",stability="STAR",nlams=15, N=1001,beta=0.3,parallel = T,nCpus=10)
#data1_fit<-XMRF(as.matrix(data3), method="PGM",lambda.path=lambda, sth=0.3,N=101,parallel = T,nCpus=10)

plotGML(data2_fit, fn="nor_pgm.gml", weight=TRUE,vars=rownames(data3))
g2_file = read.graph(file="nor_pgm.gml", format="gml")


########################################################################################################################
############################################################
############################################################    VISUALIZATION
########################################################################################################################



#bay_main2=normal,bay_main=dbt
pos=which(bay_main2$from=="Diabetes_mellitus")
wgt2=bay_main2[pos,]


plot.network <- function(structure, ht = "100vh"){
  nodes.uniq <- unique(c(structure$arcs[,1], structure$arcs[,2]))
  nodes <- data.frame(id = nodes.uniq,
                      label = nodes.uniq,
                      color = "darkturquoise",
                      shadow = TRUE,
                      size=100)
  
  edges <- data.frame(from = structure$arcs[,1],
                      to = structure$arcs[,2],
                      arrows = "to",
                      smooth = TRUE,
                      shadow = TRUE,
                      color = "black")
  return(visNetwork(nodes, edges, height = ht, width = "100%" )%>%
           visPhysics(solver = "barnesHut", 
                      barnesHut = list(avoidOverlap = 0.1,springConstant=0)))
}


graph<-plot.network(avg_main2)%>% visNodes(font = list(size=80,color="black"))

#graph

graph %>%
  visOptions(highlightNearest = TRUE)

####    $$$$$$$$$$$$$$$$$$$$$$$$$$$     MARKOV    $$$$$$$$$$$$$$$ ########################
# 
library(visNetwork)
# 
# # Convert igraph to dataframe
ig_df <- igraph::as_data_frame(g2_file, what = "both")

#weights
wgt=ig_df$edges
pos=which(wgt$from=="Diabetes_mellitus")
wgt2=wgt[pos,]

write.table(wgt2,file="nor_mark_edges.txt",row.names = F,col.names = T,quote=F,sep="\t")
##################################
# 
# ####################
# #ig_df=structure
plot.network <- function(structure, ht = "100vh"){
  nodes.uniq <- unique(structure$vertices$name)
  nodes <- data.frame(id = nodes.uniq,
                      label = nodes.uniq,
                      color = "orange",
                      shadow = TRUE,
                      size=100)
  
  edges <- data.frame(from = structure$edges$from,
                      to = structure$edges$to,
                      smooth = TRUE,
                      shadow = TRUE,
                      color = "black")
  return(visNetwork(nodes, edges, height = ht, width = "100%" )%>%
           visPhysics(solver = "barnesHut",
                      barnesHut = list(avoidOverlap = 0.1,springConstant=0)))
}


graph<-plot.network(ig_df)%>% visNodes(font = list(size=70,color="black"))
graph %>%
  visOptions(highlightNearest = TRUE)

