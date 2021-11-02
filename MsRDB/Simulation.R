library(biomformat)
library(dada2)
library("ANCOMBC")
library("dacomp")
library(phyloseq)
library(foreach)
library(doParallel)
library(doRNG)
library(ggplot2)

x = read_biom("./Data/Yatsunenko_Nature_2012/62018/all.biom")
seqtab<-as(biom_data(x), "matrix")
seqtab<-t(seqtab)
seqtab<-seqtab[,apply(seqtab, 2, function(x){mean(x>0)})>0.03]
seq <- getSequences(seqtab)
seqlist <- DNAStringSet(seq)
seqdist <- DistanceMatrix(seqlist)

mean(seqtab>0)
apply(seqtab, 2, function(x){mean(x>0)})

DataGeneration <- function(m=c(50,50),s=10, effect=20,seqtab,seqdist, measureerror=TRUE,balancedepth=TRUE) {
  ChosenSample <- sample(1:nrow(seqtab), m[1]+m[2], replace = TRUE)
  Data <- seqtab[ChosenSample,]
  non0<-which(apply(Data, 2, function(x){mean(x>0)})>0.03)
  ChooseASV <- sample(non0, s)
  EffectSize <- ceiling(runif(s, min=effect, max = effect*2))+1
  EffectRadius <- ceiling(runif(s, min=10, max = 15))
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
  
  for (i in 1:m[1])
  {
    Data[i,]=Data[i,]*ASVeffect
  }
  Z <- rep(0,m[1]+m[2])
  Z[1:m[1]] <- 1
  
  if (balancedepth==FALSE)
  {
    Depth <- ceiling(runif(m[1], min=2, max = 5))
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


SingleSimulation <- function(seqtab, seqdist, m=c(50,50), s=20, effect=20,fdr=TRUE,balancedepth=TRUE)
{
  Data <- DataGeneration(m=m,s=s,seqtab=seqtab,seqdist=seqdist,effect=effect, balancedepth=balancedepth)
  
  Error=matrix(0,nrow = 6, ncol = 2)
  
  ## MsRDB Methods k=20
  rdbresult <- msrdb(Data$Data, Data$Z, seqdist=seqdist, knn = 20)
  Error[1,1]=(sum(rdbresult$Sig)-sum(rdbresult$Sig&Data$Truth))/max(sum(rdbresult$Sig),1)
  Error[1,2]=sum(rdbresult$Sig&Data$Truth)/sum(Data$Truth)
  
  ## MsRDB Methods k=10
  rdbresult <- msrdb(Data$Data, Data$Z, seqdist=seqdist, knn = 10)
  Error[2,1]=(sum(rdbresult$Sig)-sum(rdbresult$Sig&Data$Truth))/max(sum(rdbresult$Sig),1)
  Error[2,2]=sum(rdbresult$Sig&Data$Truth)/sum(Data$Truth)
  
  ## RDB.KNN K=10
  rdb.knnresult <- rdb.knn(Data$Data, Data$Z, seqdist=seqdist, knn = 10)
  Error[3,1]=(sum(rdb.knnresult)-sum(rdb.knnresult&Data$Truth))/max(sum(rdb.knnresult),1)
  Error[3,2]=sum(rdb.knnresult&Data$Truth)/sum(Data$Truth)
  
  ## RDB.KNN K=5
  rdb.knnresult <- rdb.knn(Data$Data, Data$Z, seqdist=seqdist, knn = 5)
  Error[4,1]=(sum(rdb.knnresult)-sum(rdb.knnresult&Data$Truth))/max(sum(rdb.knnresult),1)
  Error[4,2]=sum(rdb.knnresult&Data$Truth)/sum(Data$Truth)
  
  ## RDB.KNN K=3
  rdb.knnresult <- rdb.knn(Data$Data, Data$Z, seqdist=seqdist, knn = 3)
  Error[5,1]=(sum(rdb.knnresult)-sum(rdb.knnresult&Data$Truth))/max(sum(rdb.knnresult),1)
  Error[5,2]=sum(rdb.knnresult&Data$Truth)/sum(Data$Truth)
  
  ## RDB.ASV
  rdb.fdrresult <- rdb(Data$Data, Data$Z)
  Error[6,1]=(sum(rdb.fdrresult)-sum(rdb.fdrresult&Data$Truth))/max(sum(rdb.fdrresult),1)
  Error[6,2]=sum(rdb.fdrresult&Data$Truth)/sum(Data$Truth)
  
  return(c(Error[,1],Error[,2]))
}


mlist <- c(50, 100, 200, 400)
times=100
cl <- makeCluster(5)
registerDoParallel(cl)
set.seed(123)
FinalResultsM=list()
for (i in 1:length(mlist)) {
  FinalResultsM[[i]] <- foreach(t=1:times, .combine=rbind, .packages=c('igraph','dada2')) %dorng% {
    Re=SingleSimulation(seqtab, seqdist,m=c(mlist[i],mlist[i]), effect=10)
    Re
  }
}
stopCluster(cl)

Method<-c("MsRDB(K=20)","MsRDB(K=10)","WtRDB(K=10)","WtRDB(K=5)","WtRDB(K=3)","Vanilla RDB")
ComTable=data.frame()
for (i in 1:length(mlist)) {
  Table <- FinalResultsM[[i]]
  Table <- apply(Table, 2, mean)
  Table <- data.frame(Table)
  colnames(Table) <- c("Value")
  Table <- cbind(Table,Method=c(Method,Method),ValueType=c(rep("FDR",length(Method)),rep("Power",length(Method))),SampleSize=rep(mlist[i],2*length(Method)))
  ComTable=rbind(ComTable,data.frame(Table))
}
ComTable$SampleSize <- factor(ComTable$SampleSize, levels=c("50", "100", "200","400"))
ComTable$Method <- factor(ComTable$Method, levels=Method)

ggplot(data = ComTable, aes(x = SampleSize, y = Value,fill = Method)) + 
  geom_bar(stat="identity", position = "dodge") +
  facet_grid(. ~ ValueType) +
  geom_hline(data = subset(ComTable,ValueType=="FDR"),aes(yintercept=0.1), linetype="dashed") +
  ylab("") +
  scale_fill_brewer(palette="Spectral") +
  theme_bw()






effectlist <- c(2, 4, 6, 8)
times=100
cl <- makeCluster(5)
registerDoParallel(cl)
set.seed(123)
FinalResultsEffect=list()
for (i in 1:length(effectlist)) {
  FinalResultsEffect[[i]] <- foreach(t=1:times, .combine=rbind, .packages=c('igraph','dada2')) %dorng% {
    Re=SingleSimulation(seqtab, seqdist,m=c(300,300), effect=effectlist[i])
    Re
  }
}
stopCluster(cl)

Method<-c("MsRDB(K=20)","MsRDB(K=10)","WtRDB(K=10)","WtRDB(K=5)","WtRDB(K=3)","Vanilla RDB")
ComTable=data.frame()
for (i in 1:length(effectlist)) {
  Table <- FinalResultsEffect[[i]]
  Table <- apply(Table, 2, mean)
  Table <- data.frame(Table)
  colnames(Table) <- c("Value")
  Table <- cbind(Table,Method=c(Method,Method),ValueType=c(rep("FDR",length(Method)),rep("Power",length(Method))),EffectSize=rep(effectlist[i],2*length(Method)))
  ComTable=rbind(ComTable,data.frame(Table))
}
ComTable$EffectSize <- factor(ComTable$EffectSize, levels=c("2", "4", "6","8"))
ComTable$Method <- factor(ComTable$Method, levels=Method)

ggplot(data = ComTable, aes(x = EffectSize, y = Value,fill = Method)) + 
  geom_bar(stat="identity", position = "dodge") +
  facet_grid(. ~ ValueType) +
  geom_hline(data = subset(ComTable,ValueType=="FDR"),aes(yintercept=0.1), linetype="dashed") +
  ylab("") +
  scale_fill_brewer(palette="Spectral") +
  theme_bw()








sparsitylist <- c(5,10,15,20)
times=100
cl <- makeCluster(5)
registerDoParallel(cl)
set.seed(123)
FinalResultsSparsity=list()
for (i in 1:length(sparsitylist)) {
  FinalResultsSparsity[[i]] <- foreach(t=1:times, .combine=rbind, .packages=c('igraph','dada2')) %dorng% {
    Re=SingleSimulation(seqtab, seqdist,s=sparsitylist[i], m=c(200,200), effect=8)
    Re
  }
}
stopCluster(cl)

Method<-c("MsRDB(K=20)","MsRDB(K=10)","WtRDB(K=10)","WtRDB(K=5)","WtRDB(K=3)","Vanilla RDB")
ComTable=data.frame()
for (i in 1:length(effectlist)) {
  Table <- FinalResultsSparsity[[i]]
  Table <- apply(Table, 2, mean)
  Table <- data.frame(Table)
  colnames(Table) <- c("Value")
  Table <- cbind(Table,Method=c(Method,Method),ValueType=c(rep("FDR",length(Method)),rep("Power",length(Method))),Sparsity=rep(sparsitylist[i],2*length(Method)))
  ComTable=rbind(ComTable,data.frame(Table))
}
ComTable$Sparsity <- factor(ComTable$Sparsity, levels=c("5", "10", "15","20"))
ComTable$Method <- factor(ComTable$Method, levels=Method)

ggplot(data = ComTable, aes(x = Sparsity, y = Value,fill = Method)) + 
  geom_bar(stat="identity", position = "dodge") +
  facet_grid(. ~ ValueType) +
  geom_hline(data = subset(ComTable,ValueType=="FDR"),aes(yintercept=0.1), linetype="dashed") +
  ylab("") +
  scale_fill_brewer(palette="Spectral") +
  theme_bw()


