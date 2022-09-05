source("./MsRDB/Algorithm.R")
source("./MsRDB/ASVwise.R")
source("./MsRDB/MultiATE.R")

library(biomformat)
library(dada2)
library(phyloseq)
library(foreach)
library(doParallel)
library(doRNG)
library(ggplot2)

library(ANCOMBC)
library(ALDEx2)
library(dacomp)
library(StructFDR)

x = read_biom("./Data/Yatsunenko_Nature_2012/62018/all.biom")
seqtab<-as(biom_data(x), "matrix")
seqtab<-t(seqtab)
seqtab<-seqtab[,apply(seqtab, 2, function(x){mean(x>0)})>0.1]
seq <- getSequences(seqtab)
seqlist <- DNAStringSet(seq)
seqdist <- DistanceMatrix(seqlist)
colnames(seqdist) <- colnames(seqtab)
rownames(seqdist) <- colnames(seqtab)
tree <- phangorn::upgma(seqdist)
taxa<-read.table("./Data/Yatsunenko_Nature_2012/62018/taxa.txt", sep="\t")
taxa<-taxa[colnames(seqtab),]
rm(x)

DataGeneration <- function(method,m=c(50,50), s=10, rank="Genus", difftaxa=NULL, effect=20, seqtab, seqdist, taxa, measureerror=TRUE, balancedepth=TRUE) {
  ChosenSample <- sample(1:nrow(seqtab), m[1]+m[2], replace = FALSE)
  Data <- seqtab[ChosenSample,]
  if (method == "NN") {
    non0<-which(apply(Data, 2, function(x){mean(x>0)})>0.03)
    ChooseASV <- sample(non0, s)
    EffectSize <- ceiling(runif(s, min=effect, max = effect*2))+1
    EffectRadius <- ceiling(runif(s, min=5, max = 10))
    ASVeffect=rep(1,ncol(seqtab))
    for (j in 1:s) {
      index <- ChooseASV[j]
      tempdist <- seqdist[index,]
      tempdist[index] <- -1
      orderdist <- order(tempdist)
      tsig<-orderdist[1:EffectRadius[j]]
      tsig<-tsig[tsig %in% non0]
      ASVeffect[tsig]=EffectSize[j]
    }
  } else if (method == "Taxa") {
    ASVeffect=rep(1,ncol(seqtab))
    for (ta in difftaxa) {
      ChooseASV=which(taxa[,rank]==ta)
      EffectSize <- ceiling(runif(1, min=effect, max = effect*2))+1
      ASVeffect[ChooseASV]=EffectSize
    }
  }
  
  
  for (i in 1:m[1])
  {
    Data[i,]=Data[i,]*ASVeffect
  }
  Z <- rep(0,m[1]+m[2])
  Z[1:m[1]] <- 1
  
  if (balancedepth==FALSE)
  {
    Depth <- ceiling(runif(m[1], min=1, max = 3))
    for (i in 1:m[1]) {
      Data[i,]=Data[i,]*Depth[i]
    }
  }
  
  if (measureerror==TRUE)
  {
    MError <- ceiling(runif(ncol(seqtab), min=1, max = 10))
    for (i in 1:ncol(seqtab)) {
      Data[,i]=Data[,i]*MError[i]
    }
  }
  
  Data <- list(Data=Data, Z=Z, Truth=ASVeffect>1)
  return(Data)
}

SingleSimulation <- function(method, seqtab, seqdist, taxa, tree, m=c(50,50), s=20, rank="Genus", difftaxa=NULL, effect=20,measureerror=FALSE, balancedepth=TRUE)
{
  Data <- DataGeneration(method, m=m, s=s, rank, difftaxa, seqtab=seqtab, seqdist=seqdist, taxa=taxa, effect=effect, measureerror=measureerror, balancedepth=balancedepth)
  
  Error=matrix(0,nrow = 6, ncol = 2)
  
  ## MsRDB
  rdbresult <- msrdb(Data$Data, Data$Z, seqdist=seqdist, knn = 10)
  Error[1,1]=(sum(rdbresult$Sig,na.rm =TRUE)-sum(rdbresult$Sig&Data$Truth))/max(sum(rdbresult$Sig,na.rm =TRUE),1)
  Error[1,2]=sum(rdbresult$Sig&Data$Truth,na.rm =TRUE)/sum(Data$Truth)
  
  ## RDB
  rdb.fdrresult <- rdb(Data$Data, Data$Z)
  Error[2,1]=(sum(rdb.fdrresult,na.rm =TRUE)-sum(rdb.fdrresult&Data$Truth))/max(sum(rdb.fdrresult,na.rm =TRUE),1)
  Error[2,2]=sum(rdb.fdrresult&Data$Truth,na.rm =TRUE)/sum(Data$Truth)
  
  ## ANCOMBC
  otu_mat = t(Data$Data+1)
  meta = data.frame(group = Data$Z, row.names = colnames(otu_mat),stringsAsFactors = FALSE)
  tax_mat = as.matrix(taxa)
  OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
  META = sample_data(meta)
  TAX = tax_table(tax_mat)
  physeq = phyloseq(OTU, META, TAX)
  Re=ancombc(phyloseq=physeq,formula="group",p_adj_method="holm",alpha = 0.1)
  Re=Re$res$diff_abn[,1]
  Error[3,1]=(sum(Re,na.rm =TRUE)-sum(Re&Data$Truth,na.rm =TRUE))/max(sum(Re,na.rm =TRUE),1)
  Error[3,2]=sum(Re&Data$Truth,na.rm =TRUE)/sum(Data$Truth,na.rm =TRUE)
  
  ## DACOMP
  ref_select = dacomp.select_references(X = Data$Data+1)
  Selected_References = ref_select$selected_references
  res_welch = dacomp.test(X = Data$Data+1, y = factor(Data$Z), ind_reference_taxa = Selected_References, test = DACOMP.TEST.NAME.WILCOXON, q=0.1)
  Padj=p.adjust(res_welch$p.values.test, method = "BH")
  Re=Padj<0.1
  Error[4,1]=(sum(Re,na.rm =TRUE)-sum(Re&Data$Truth,na.rm =TRUE))/max(sum(Re,na.rm =TRUE),1)
  Error[4,2]=sum(Re&Data$Truth,na.rm =TRUE)/sum(Data$Truth,na.rm =TRUE)
  
  ##ALDEx2
  reduceData=Data$Data
  index=apply(reduceData, 2, sum)>0
  reduceData = reduceData[,index]
  reduceTruth = Data$Truth[index]
  Xmatrix <- model.matrix(as.formula("~ group"), meta)
  aclr <- aldex.clr(t(reduceData), Xmatrix)
  aldexRe <- aldex.glm(aclr)
  Re <- aldexRe$`model.group Pr(>|t|).BH`<0.1
  Error[5,1]=(sum(Re,na.rm =TRUE)-sum(Re&reduceTruth,na.rm =TRUE))/max(sum(Re,na.rm =TRUE),1)
  Error[5,2]=sum(Re&reduceTruth,na.rm =TRUE)/sum(reduceTruth,na.rm =TRUE)
  
  ##StructFDR
  tree.fdr.obj=MicrobiomeSeqTreeFDR (t(Data$Data), tree, data.frame(group=Data$Z), grp.name='group', raw.count=TRUE)
  Re=tree.fdr.obj$p.adj<0.1
  Error[6,1]=(sum(Re,na.rm =TRUE)-sum(Re&Data$Truth,na.rm =TRUE))/max(sum(Re,na.rm =TRUE),1)
  Error[6,2]=sum(Re&Data$Truth,na.rm =TRUE)/sum(Data$Truth,na.rm =TRUE)
  
  return(c(Error[,1],Error[,2]))
}

##NN

mlist <- c(50, 100, 200)
times=100
cl <- makeCluster(6)
registerDoParallel(cl)
set.seed(123)
FinalResultsM=list()
for (i in 1:length(mlist)) {
  FinalResultsM[[i]] <- foreach(t=1:times, .combine=rbind, .packages=c('igraph','dada2','ANCOMBC','ALDEx2','dacomp','StructFDR','phyloseq')) %dorng% {
    Re=SingleSimulation(method="NN", seqtab, seqdist, taxa, tree, m=c(mlist[i],mlist[i]), s=5, effect=20, balancedepth=TRUE)
    Re
  }
}
stopCluster(cl)

Method<-c("MsRDB","RDB","ANCOMBC","DACOMP","ALDEx2","StructFDR")
ComTable=data.frame()
for (i in 1:length(mlist)) {
  Table <- FinalResultsM[[i]]
  Table <- apply(Table, 2, mean)
  Table <- data.frame(Table)
  colnames(Table) <- c("Value")
  Table <- cbind(Table,Method=c(Method,Method),ValueType=c(rep("FDR",length(Method)),rep("Power",length(Method))),SampleSize=rep(mlist[i],2*length(Method)))
  ComTable=rbind(ComTable,data.frame(Table))
}
ComTable$SampleSize <- factor(ComTable$SampleSize, levels=c("50", "100", "200"))
ComTable$Method <- factor(ComTable$Method, levels=Method)

ggplot(data = ComTable, aes(x = SampleSize, y = Value,fill = Method)) + 
  geom_bar(stat="identity", position = "dodge") +
  facet_grid(. ~ ValueType) +
  geom_hline(data = subset(ComTable,ValueType=="FDR"),aes(yintercept=0.1), linetype="dashed") +
  ylab("") +
  scale_fill_brewer(palette="Set2") +
  theme_bw() +
  coord_flip()

##Genus

mlist <- c(50, 100, 200)
times=100
difftaxa=c("Ruminococcus","Bacteroides","Faecalibacterium","Oscillibacter", "Blautia")
cl <- makeCluster(6)
registerDoParallel(cl)
set.seed(123)
FinalResultsM=list()
for (i in 1:length(mlist)) {
  FinalResultsM[[i]] <- foreach(t=1:times, .combine=rbind, .packages=c('igraph','dada2','ANCOMBC','ALDEx2','dacomp','StructFDR','phyloseq')) %dorng% {
    Re=SingleSimulation(method="Taxa", seqtab, seqdist, taxa, tree, m=c(mlist[i],mlist[i]), rank="Genus", difftaxa=difftaxa, effect=20, balancedepth=TRUE)
    Re
  }
}
stopCluster(cl)

Method<-c("MsRDB","RDB","ANCOMBC","DACOMP","ALDEx2","StructFDR")
ComTable=data.frame()
for (i in 1:length(mlist)) {
  Table <- FinalResultsM[[i]]
  Table <- apply(Table, 2, mean)
  Table <- data.frame(Table)
  colnames(Table) <- c("Value")
  Table <- cbind(Table,Method=c(Method,Method),ValueType=c(rep("FDR",length(Method)),rep("Power",length(Method))),SampleSize=rep(mlist[i],2*length(Method)))
  ComTable=rbind(ComTable,data.frame(Table))
}
ComTable$SampleSize <- factor(ComTable$SampleSize, levels=c("50", "100", "200"))
ComTable$Method <- factor(ComTable$Method, levels=Method)

ggplot(data = ComTable, aes(x = SampleSize, y = Value,fill = Method)) + 
  geom_bar(stat="identity", position = "dodge") +
  facet_grid(. ~ ValueType) +
  geom_hline(data = subset(ComTable,ValueType=="FDR"),aes(yintercept=0.1), linetype="dashed") +
  ylab("") +
  scale_fill_brewer(palette="Set2") +
  theme_bw() +
  coord_flip()

##Class

mlist <- c(50, 100, 200)
times=100
difftaxa=c("Bacteroidia","Bacilli")
cl <- makeCluster(6)
registerDoParallel(cl)
set.seed(123)
FinalResultsM=list()
for (i in 1:length(mlist)) {
  FinalResultsM[[i]] <- foreach(t=1:times, .combine=rbind, .packages=c('igraph','dada2','ANCOMBC','ALDEx2','dacomp','StructFDR','phyloseq')) %dorng% {
    Re=SingleSimulation(method="Taxa", seqtab, seqdist, taxa, tree, m=c(mlist[i],mlist[i]), rank="Class", difftaxa=difftaxa, effect=20, balancedepth=TRUE)
    Re
  }
}
stopCluster(cl)

Method<-c("MsRDB","RDB","ANCOMBC","DACOMP","ALDEx2","StructFDR")
ComTable=data.frame()
for (i in 1:length(mlist)) {
  Table <- FinalResultsM[[i]]
  Table <- apply(Table, 2, mean)
  Table <- data.frame(Table)
  colnames(Table) <- c("Value")
  Table <- cbind(Table,Method=c(Method,Method),ValueType=c(rep("FDR",length(Method)),rep("Power",length(Method))),SampleSize=rep(mlist[i],2*length(Method)))
  ComTable=rbind(ComTable,data.frame(Table))
}
ComTable$SampleSize <- factor(ComTable$SampleSize, levels=c("50", "100", "200"))
ComTable$Method <- factor(ComTable$Method, levels=Method)

ggplot(data = ComTable, aes(x = SampleSize, y = Value,fill = Method)) + 
  geom_bar(stat="identity", position = "dodge") +
  facet_grid(. ~ ValueType) +
  geom_hline(data = subset(ComTable,ValueType=="FDR"),aes(yintercept=0.1), linetype="dashed") +
  ylab("") +
  scale_fill_brewer(palette="Set2") +
  theme_bw() +
  coord_flip()


##Measure Error

mlist <- c(50, 100, 200)
times=100
cl <- makeCluster(6)
registerDoParallel(cl)
set.seed(123)
FinalResultsM=list()
for (i in 1:length(mlist)) {
  FinalResultsM[[i]] <- foreach(t=1:times, .combine=rbind, .packages=c('igraph','dada2','ANCOMBC','ALDEx2','dacomp','StructFDR','phyloseq')) %dorng% {
    Re=SingleSimulation(method="NN", seqtab, seqdist, taxa, tree, m=c(mlist[i],mlist[i]), s=5, effect=20, measureerror=TRUE, balancedepth=TRUE)
    Re
  }
}
stopCluster(cl)

Method<-c("MsRDB","RDB","ANCOMBC","DACOMP","ALDEx2","StructFDR")
ComTable=data.frame()
for (i in 1:length(mlist)) {
  Table <- FinalResultsM[[i]]
  Table <- apply(Table, 2, mean)
  Table <- data.frame(Table)
  colnames(Table) <- c("Value")
  Table <- cbind(Table,Method=c(Method,Method),ValueType=c(rep("FDR",length(Method)),rep("Power",length(Method))),SampleSize=rep(mlist[i],2*length(Method)))
  ComTable=rbind(ComTable,data.frame(Table))
}
ComTable$SampleSize <- factor(ComTable$SampleSize, levels=c("50", "100", "200"))
ComTable$Method <- factor(ComTable$Method, levels=Method)

ggplot(data = ComTable, aes(x = SampleSize, y = Value,fill = Method)) + 
  geom_bar(stat="identity", position = "dodge") +
  facet_grid(. ~ ValueType) +
  geom_hline(data = subset(ComTable,ValueType=="FDR"),aes(yintercept=0.1), linetype="dashed") +
  ylab("") +
  scale_fill_brewer(palette="Set2") +
  theme_bw() +
  coord_flip()


##Unbalanced Depth

mlist <- c(50, 100, 200)
times=100
cl <- makeCluster(6)
registerDoParallel(cl)
set.seed(123)
FinalResultsM=list()
for (i in 1:length(mlist)) {
  FinalResultsM[[i]] <- foreach(t=1:times, .combine=rbind, .packages=c('igraph','dada2','ANCOMBC','ALDEx2','dacomp','StructFDR','phyloseq')) %dorng% {
    Re=SingleSimulation(method="NN", seqtab, seqdist, taxa, tree, m=c(mlist[i],mlist[i]), s=5, effect=20, measureerror=TRUE, balancedepth=FALSE)
    Re
  }
}
stopCluster(cl)

Method<-c("MsRDB","RDB","ANCOMBC","DACOMP","ALDEx2","StructFDR")
ComTable=data.frame()
for (i in 1:length(mlist)) {
  Table <- FinalResultsM[[i]]
  Table <- apply(Table, 2, mean)
  Table <- data.frame(Table)
  colnames(Table) <- c("Value")
  Table <- cbind(Table,Method=c(Method,Method),ValueType=c(rep("FDR",length(Method)),rep("Power",length(Method))),SampleSize=rep(mlist[i],2*length(Method)))
  ComTable=rbind(ComTable,data.frame(Table))
}
ComTable$SampleSize <- factor(ComTable$SampleSize, levels=c("50", "100", "200"))
ComTable$Method <- factor(ComTable$Method, levels=Method)

ggplot(data = ComTable, aes(x = SampleSize, y = Value,fill = Method)) + 
  geom_bar(stat="identity", position = "dodge") +
  facet_grid(. ~ ValueType) +
  geom_hline(data = subset(ComTable,ValueType=="FDR"),aes(yintercept=0.1), linetype="dashed") +
  ylab("") +
  scale_fill_brewer(palette="Set2") +
  theme_bw() +
  coord_flip()
