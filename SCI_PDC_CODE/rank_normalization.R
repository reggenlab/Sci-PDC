library(parallel)
options("mc.cores" = 4L)
library("data.table")
library(stringr)
mc.cores=4
sub_file<-list.files(pattern="comb_corr*",recursive=F,full.names = FALSE)
ldf <- mclapply(sub_file, function(a) data.frame(fread(a,header = T,stringsAsFactors = F,check.names = F), row.names=1))
ldf2 <- mclapply(sub_file, function(x) {
    z <- data.frame(fread(x,header = T,stringsAsFactors = F,check.names = F), row.names=1)
    abc2=abs(z)
    disease<-colnames(abc2)
    filename<-str_sub(x,1,nchar(x)-9)
    dim(z);
    j = 0;
    rankadj <- matrix(0,nrow(abc2),ncol(abc2),dimnames=list(rownames(abc2),colnames(abc2)))
    for(x in colnames(abc2)){
        j = j + 1;
        res <- lapply(ldf, function(a) a[ , grepl( x, colnames(a) ) ])
        null<-do.call(cbind,res)
        null2<-abs(null)
        rank_ref<-apply(null2, 2, function(y) rank(y,na.last = "keep",ties.method = "min")/length(y))
        userp<-abc2
        rank_user<-apply(userp, 2, function(y) rank(y,na.last = "keep",ties.method = "min")/length(y))
        print(j)
        for( i in 1:nrow(rank_user)){
            pos = which(rank_ref[i,] < rank_user[i,j]) ;
            rankadj[i,j] =1- (length(pos)/ncol(rank_ref)) ;
        }
    }
        null_path<-paste("/storage/vibhor/Phds/Madhu/human_corr/disease_rank_adjusted/",filename,".txt", sep = "");
        write.table(rankadj,null_path,quote=F,sep="\t")
 
}) 
