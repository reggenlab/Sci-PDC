
#######################     case study 1 scores   ##############################

library(UniPath)
source("PathwayCooccurrence.R")
######################  PART I ####################


############# READ GENE SET FILE
geneSet = read.csv("file.csv",header=F,stringsAsFactors = F,row.names = 1,fill=T,flush=TRUE,sep=",")

######################################################################################

######################################################################################


########### NULL MODEL OF AVRAGE CELL TYPE EXPRESSION
human_null_data = read.table("human_null_model_5000.csv",sep=",",row.names=1,stringsAsFactors = F,header=T)

########## CONVERTING GENE EXPRESSION TO P-VALUES
Pval = binorm(human_null_data)


######### COMBINING P-VALUES IN GENE SET
combp_ref = UniPath::combine(geneSet,human_null_data,rownames(human_null_data),Pval,thr=5)


####################### PART II PANCREATIC DATA  ############################################################################

############# EXPRESSION DATA
segerstolpe <- readRDS("segerstolpe.rds")

expression_matrix = segerstolpe@assays$data$counts

Pval1 = binorm(expression_matrix)
combp = UniPath::combine(geneSet,expression_matrix,rownames(expression_matrix),Pval1,thr=5)

#' Adjusting of combine p-values using null model
#'
#' @param combp combined p-value matrix obtained using gene expression matrix
#' @param combp_ref combined p-value matrix obtained using null model
#'
#' @return A list of 3 matrices, One is absolute p-value matrix, second is raw adjusted p-value matrix and third is log transformed adjusted p-value matrix
#' @export
#'
#' @examples
#' adjust()


adjust <- function(combp,combp_ref){
  
  adjpva <- matrix(0,nrow(combp),ncol(combp),dimnames=list(rownames(combp),colnames(combp)))
  for( i in 1:nrow(combp)){
    for(j in 1:ncol(combp)){
      pos = which(combp_ref[i,] < combp[i,j]) ;
      adjpva[i,j] = nrow(as.matrix(pos))/5000 ;
    }
  }
  adjpva1 = (1-adjpva)
  adjpva2 = -log2(adjpva1 +.0001)
  scores = list(adjpva,adjpva1,adjpva2)
  names(scores) <- c("adjpva","adjpvaraw","adjpvalog")
  
  return(scores)
}

set.seed(145)
########## ADJUSTED PATHWAY SCORES
scores = adjust(combp,combp_ref)

# columns annotation data
cl = as.data.frame(segerstolpe@colData)

mt = as.matrix(paste(cl[,2],cl[,3],sep="-"))
rownames(mt) = rownames(cl)

meta = mt[colnames(scores$adjpva),]


pos = which(meta=="beta-normal")
data =  scores$adjpvaraw[,pos]
write.csv(data,file = "no_beta_padj_scores.csv",quote = F)
xc_acinar_n= cooccurrence(combp_ref,data)
corr_file=xc_acinar_n$correlation
corr_file2<-corr_file[1812:3140,1:1811]
p_file=xc_acinar_n$pval
p_file2<-p_file[1812:3140,1:1811]
write.csv(corr_file,file = "no_beta_corr.csv")
write.csv(p_file,file = "no_beta_pval.csv")



##############
#####################################################
###################################
#######################     case study 1 scores   ##############################

#######################     case study 2 scores   ##############################

########## ADJUSTED PATHWAY SCORES  HUMAN


baron <- readRDS("baron-human.rds")

expres<-baron@assays$data$counts

Pval_m_p = binorm(expres)
combp_m_p = UniPath::combine(geneSet,expres,rownames(expres),Pval_m_p,thr=5)


scores_m = adjust(combp_m_p,combp_ref)

# columns annotation data
cl_m = as.data.frame(baron@colData)

mt_m = as.matrix(cl_m[,2])
rownames(mt_m) = rownames(cl_m)

meta_m = mt_m[colnames(scores_m$adjpva),]

#################################################
pos = which(meta_m=="beta")

data_m =  scores_m$adjpvaraw[,pos]
write.table(data_m,'baron_h_beta_scores.txt',quote=F)
xc_beta_m= cooccurrence(combp_ref_m,data_m)

corr_file=xc_beta_m$correlation
corr_file2<-corr_file[1812:3140,1:1811]
p_file=xc_beta_m$pval
p_file2<-p_file[1812:3140,1:1811]
write.csv(corr_file2,file = "h_beta_corr.csv")
write.csv(p_file2,file = "h_beta_pval.csv")



########## ADJUSTED PATHWAY SCORES  MOUSE

########### NULL MODEL OF AVRAGE CELL TYPE EXPRESSION
mouse_null_data = read.table("mouse_null_model_5000_cells.csv",sep=",",row.names=1,stringsAsFactors = F,header=T)

########## CONVERTING GENE EXPRESSION TO P-VALUES
Pval_mouse = binorm(mouse_null_data)


######### COMBINING P-VALUES IN GENE SET
combp_ref_m = combine(geneSet,mouse_null_data,rownames(mouse_null_data),Pval_mouse,thr=5)

################################### MOUSE DATA ##############3
############# EXPRESSION DATA
baron <- readRDS("baron-mouse.rds")

expres<-baron@assays$data$counts

Pval_m_p = binorm(expres)
combp_m_p = UniPath::combine(geneSet,expres,rownames(expres),Pval_m_p,thr=5)

########## ADJUSTED PATHWAY SCORES
scores_m = adjust(combp_m_p,combp_ref_m)

# columns annotation data
cl_m = as.data.frame(baron@colData)

mt_m = as.matrix(cl_m[,2])
rownames(mt_m) = rownames(cl_m)

meta_m = mt_m[colnames(scores_m$adjpva),]

#################################################
pos = which(meta_m=="beta")
baron_m_beta_scores.txt
data_m =  scores_m$adjpvaraw[,pos]

write.table(data_m,'baron_m_beta_scores.txt',quote=F)

xc_beta_m= cooccurrence(combp_ref_m,data_m)

corr_file=xc_beta_m$correlation
corr_file2<-corr_file[1812:3140,1:1811]
p_file=xc_beta_m$pval
p_file2<-p_file[1812:3140,1:1811]
write.csv(corr_file2,file = "m_beta_corr.csv")
write.csv(p_file2,file = "m_beta_pval.csv")



#######################     case study 2 scores   ##############################


#######################     case study 3 scores   ##############################


expression_matrix = read.table('exp_adult_beta.txt',header=T,sep=',',row.names=1,stringsAsFactors = F)

Pval1 = binorm(expression_matrix)
combp = UniPath::combine(geneSet,expression_matrix,rownames(expression_matrix),Pval1,thr=5)

########## ADJUSTED PATHWAY SCORES
scores = adjust(combp,combp_ref)

data =  scores$adjpvaraw
write.csv(data,file = "un_a_scores.csv",quote = F)
xc_ductal_n= cooccurrence(combp_ref,data)
corr_file=xc_ductal_n$correlation
corr_file2<-corr_file[1812:3140,1:1811]
p_file=xc_ductal_n$pval
p_file2<-p_file[1812:3140,1:1811]
write.csv(corr_file2,file = "a_beta_corr.csv")
write.csv(p_file2,file = "a_beta_pval.csv")


#####################################################
expression_matrix = read.table('exp_old_beta.txt',header=T,sep=',',row.names=1,stringsAsFactors = F)

Pval1 = binorm(expression_matrix)
combp = UniPath::combine(geneSet,expression_matrix,rownames(expression_matrix),Pval1,thr=5)

########## ADJUSTED PATHWAY SCORES
scores = adjust(combp,combp_ref)



data =  scores$adjpvaraw
write.csv(data,file = "un_o_scores.csv",quote = F)
xc_ductal_n= cooccurrence(combp_ref,data)
corr_file=xc_ductal_n$correlation
corr_file2<-corr_file[1812:3140,1:1811]
p_file=xc_ductal_n$pval
p_file2<-p_file[1812:3140,1:1811]
write.csv(corr_file2,file = "o_beta_corr.csv")
write.csv(p_file2,file = "o_beta_pval.csv")


#####################################################

#####################################################
expression_matrix = read.table('exp_young_beta.txt',header=T,sep=',',row.names=1,stringsAsFactors = F)

Pval1 = binorm(expression_matrix)
combp = UniPath::combine(geneSet,expression_matrix,rownames(expression_matrix),Pval1,thr=5)

########## ADJUSTED PATHWAY SCORES
scores = adjust(combp,combp_ref)


data =  scores$adjpvaraw
write.csv(data,file = "un_y_scores.csv",quote = F)

xc_ductal_n= cooccurrence(combp_ref,data)
corr_file=xc_ductal_n$correlation
corr_file2<-corr_file[1812:3140,1:1811]
p_file=xc_ductal_n$pval
p_file2<-p_file[1812:3140,1:1811]
write.csv(corr_file2,file = "y_beta_corr.csv")
write.csv(p_file2,file = "y_beta_pval.csv")


#####################################################



#######################     case study 3 scores   ##############################
