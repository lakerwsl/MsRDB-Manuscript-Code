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
library(tidyverse)
library(microbiome)

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

SingleSimulation <- function(seqtab, seqdist, taxa, tree, m=c(50,50), s=20, rank="Genus", effect=20, measureerror=FALSE, balancedepth=TRUE)
{
  Data <- DataGeneration(method='NN', m=m, s=s, seqtab=seqtab, seqdist=seqdist, taxa=taxa, effect=effect, measureerror=measureerror, balancedepth=balancedepth)
  TrueTaxa = unique(taxa[Data$Truth,rank])
  
  taxatab <-as.data.frame(t(Data$Data))
  taxatab <-cbind(taxatab,taxa[,rank,drop=FALSE])
  taxatab <- taxatab %>% group_by_at(rank) %>% summarise_all(sum)
  taxanames <- as.vector(taxatab[,1][[1]])
  taxatab <- as.matrix(t(taxatab[,-1]))
  colnames(taxatab) <- taxanames
  
  Error=matrix(0,nrow = 6, ncol = 2)
  
  ## MsRDB
  rdbresult <- msrdb(Data$Data, Data$Z, seqdist=seqdist, knn = 10)
  ReTaxa = unique(taxa[rdbresult$Sig,rank])
  Error[1,1]=(length(ReTaxa)-sum(ReTaxa %in% TrueTaxa))/max(length(ReTaxa),1)
  Error[1,2]=sum(ReTaxa %in% TrueTaxa)/length(TrueTaxa)
  
  ## RDB-ASV
  rdb.fdrresult <- rdb(Data$Data, Data$Z)
  ReTaxa = unique(taxa[rdb.fdrresult,rank])
  Error[2,1]=(length(ReTaxa)-sum(ReTaxa %in% TrueTaxa))/max(length(ReTaxa),1)
  Error[2,2]=sum(ReTaxa %in% TrueTaxa)/length(TrueTaxa)
  
  ## RDB-Taxa
  rdb.fdrresult <- rdb(taxatab, Data$Z)
  ReTaxa = unique(taxanames[rdb.fdrresult])
  Error[3,1]=(length(ReTaxa)-sum(ReTaxa %in% TrueTaxa))/max(length(ReTaxa),1)
  Error[3,2]=sum(ReTaxa %in% TrueTaxa)/length(TrueTaxa)
  
  ## ANCOMBC-ASV
  otu_mat = t(Data$Data+1)
  meta = data.frame(group = Data$Z, row.names = colnames(otu_mat),stringsAsFactors = FALSE)
  tax_mat = as.matrix(taxa)
  OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
  META = sample_data(meta)
  TAX = tax_table(tax_mat)
  physeq = phyloseq(OTU, META, TAX)
  Re=ancombc(phyloseq=physeq,formula="group",p_adj_method="holm",alpha = 0.1)
  Re=Re$res$diff_abn[,1]
  ReTaxa = unique(taxa[Re,rank])
  Error[4,1]=(length(ReTaxa)-sum(ReTaxa %in% TrueTaxa))/max(length(ReTaxa),1)
  Error[4,2]=sum(ReTaxa %in% TrueTaxa)/length(TrueTaxa)
  
  ## ANCOMBC-Taxa
  physeq = aggregate_taxa(physeq, rank)
  Re=ancombc(phyloseq=physeq,formula="group",p_adj_method="holm",alpha = 0.1)
  ReTaxa = rownames(Re$res$diff_abn)[Re$res$diff_abn[,1]]
  if ("Unknown" %in% ReTaxa) {
    ReTaxa[ReTaxa=="Unknown"]=NA
  }
  Error[5,1]=(length(ReTaxa)-sum(ReTaxa %in% TrueTaxa))/max(length(ReTaxa),1)
  Error[5,2]=sum(ReTaxa %in% TrueTaxa)/length(TrueTaxa)
  
  ##StructFDR
  tree.fdr.obj=MicrobiomeSeqTreeFDR (t(Data$Data), tree, data.frame(group=Data$Z), grp.name='group', raw.count=TRUE)
  Re=tree.fdr.obj$p.adj<0.1
  ReTaxa = unique(taxa[Re,rank])
  Error[6,1]=(length(ReTaxa)-sum(ReTaxa %in% TrueTaxa))/max(length(ReTaxa),1)
  Error[6,2]=sum(ReTaxa %in% TrueTaxa)/length(TrueTaxa)
  
  return(c(Error[,1],Error[,2]))
}

##Genus

mlist <- c(50, 100, 200)
times=100
cl <- makeCluster(8)
registerDoParallel(cl)
set.seed(123)
FinalResultsM=list()
for (i in 1:length(mlist)) {
  FinalResultsM[[i]] <- foreach(t=1:times, .combine=rbind, .packages=c('igraph','dada2','ANCOMBC','microbiome','tidyverse','StructFDR','phyloseq')) %dorng% {
    Re=SingleSimulation(seqtab, seqdist, taxa, tree, m=c(mlist[i],mlist[i]), s=10, effect=20, balancedepth=TRUE)
    Re
  }
}
stopCluster(cl)

Method<-c("MsRDB","RDB-ASV","RDB-Taxa","ANCOMBC-ASV","ANCOMBC-Taxa","StructFDR")
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
  scale_fill_brewer(palette="Accent") +
  theme_bw() +
  coord_flip()

##Family

mlist <- c(50, 100, 200)
times=100
cl <- makeCluster(8)
registerDoParallel(cl)
set.seed(123)
FinalResultsM=list()
for (i in 1:length(mlist)) {
  FinalResultsM[[i]] <- foreach(t=1:times, .combine=rbind, .packages=c('igraph','dada2','ANCOMBC','microbiome','tidyverse','StructFDR','phyloseq')) %dorng% {
    Re=SingleSimulation(seqtab, seqdist, taxa, tree, m=c(mlist[i],mlist[i]), s=10, rank="Family", effect=20, balancedepth=TRUE)
    Re
  }
}
stopCluster(cl)

Method<-c("MsRDB","RDB-ASV","RDB-Taxa","ANCOMBC-ASV","ANCOMBC-Taxa","StructFDR")
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
  scale_fill_brewer(palette="Accent") +
  theme_bw() +
  coord_flip()



