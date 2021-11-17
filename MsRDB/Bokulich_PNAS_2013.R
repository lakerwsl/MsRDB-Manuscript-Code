source("./MsRDB/Algorithm.R")
source("./MsRDB/ASVwise.R")
source("./MsRDB/MultiATE.R")

library(biomformat)
library(tidyverse)
library(phyloseq)
library(ggplot2)
library(dada2)
library(DECIPHER)
library(microbiome)

########################################
## Read Data
########################################

x = read_biom("./Data/Bokulich_PNAS_2013/59547/all.biom")
seqtab<-as(biom_data(x), "matrix")
seqtab<-t(seqtab)

#taxa <- assignTaxonomy(seqtab, "../Tax/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
taxa<-read.table("./Data/Bokulich_PNAS_2013/59547/taxa.txt", sep="\t")

seq <- getSequences(seqtab)
seqlist <- DNAStringSet(seq)
seqdist <- DistanceMatrix(seqlist)


meta = read.delim("./Data/Bokulich_PNAS_2013/2019_20180418-110417.txt")
meta <- meta[meta$sample_name %in% row.names(seqtab),]
rownames(meta) <- meta$sample_name

########################################
## Group into Different Rank
########################################

sampleinterest <- meta$sample_name[meta$cultivar %in% c("Chardonnay")]

##Genus
genustab <-as.data.frame(t(seqtab[sampleinterest,]))
genustab <-cbind(genustab,taxa[row.names(genustab),c("Genus"),drop=FALSE])
genustab <- genustab %>% group_by(Genus) %>% summarise_all(sum)
rowsum=apply(genustab[,-1],1,sum)
genustab <- genustab[rowsum!=0,]

genusunassignratio<-genustab[nrow(genustab),2:ncol(genustab)]/apply(genustab[,2:ncol(genustab)], 2, sum)
mean(as.numeric(genusunassignratio), na.rm=TRUE)

##Family

familytab <-as.data.frame(t(seqtab[sampleinterest,]))
familytab <-cbind(familytab,taxa[row.names(familytab),c("Family"),drop=FALSE])
familytab <- familytab %>% group_by(Family) %>% summarise_all(sum)
rowsum=apply(familytab[,-1],1,sum)
familytab <- familytab[rowsum!=0,]

familyunassignratio<-familytab[nrow(familytab),2:ncol(familytab)]/apply(familytab[,2:ncol(familytab)], 2, sum)
mean(as.numeric(familyunassignratio), na.rm=TRUE)

##Order

ordertab <-as.data.frame(t(seqtab[sampleinterest,]))
ordertab <-cbind(ordertab,taxa[row.names(ordertab),c("Order"),drop=FALSE])
ordertab <- ordertab %>% group_by(Order) %>% summarise_all(sum)
rowsum=apply(ordertab[,-1],1,sum)
ordertab <- ordertab[rowsum!=0,]

orderunassignratio<-ordertab[nrow(ordertab),2:ncol(ordertab)]/apply(ordertab[,2:ncol(ordertab)], 2, sum)
mean(as.numeric(orderunassignratio), na.rm=TRUE)

n=length(genusunassignratio)
df=data.frame(x=c(t(genusunassignratio),t(familyunassignratio),t(orderunassignratio)),y=rep(c("Genus","Family","Order"),rep(n,3)))
df$y_f = factor(df$y, levels=c("Genus","Family","Order"))

## Figure 6A

ggplot(df, aes(x=x)) +
  geom_histogram(position="identity", binwidth=0.05, alpha=0.5, color='#DE94C8',fill='#DE94C8') +
  theme_bw() +
  xlab('Proportion of NA') +
  ylab('Number of Samples') +
  xlim(-0.05, 1.05) +
  facet_grid(~ y_f)

########################################
## ASV-wise analysis
########################################

CompareGroup <- list(c("Sonoma","Napa"),c("Central_Coast","Napa"),c("Central_Coast","Sonoma"))
SampleList <- meta$cultivar %in% c("Chardonnay")


ASVwiseResult <- list()
for (i in 1:length(CompareGroup)) {
  GroupInterest <- CompareGroup[[i]]
  Group1<-meta$sample_name[meta$area==GroupInterest[1]&SampleList]
  Group2<-meta$sample_name[meta$area==GroupInterest[2]&SampleList]
  
  seqtabInterest=seqtab[c(Group1,Group2),]
  Z=c(rep(1,length(Group1)),rep(0,length(Group2)))
  X=NULL
  
  rdbRe <- rdb(seqtabInterest,Z,X=X)
  
  ASVwiseResult[[i]] <- rdbRe
}

########################################
## Taxon-wise analysis
########################################

##Order

OrderwiseResult <- list()
for (i in 1:length(CompareGroup)) {
  GroupInterest <- CompareGroup[[i]]
  Group1<-meta$sample_name[meta$area==GroupInterest[1]&SampleList]
  Group2<-meta$sample_name[meta$area==GroupInterest[2]&SampleList]
  
  ordertabInterest=ordertab[,c(Group1,Group2)]
  ordertabInterest=as.matrix(ordertabInterest)
  rownames(ordertabInterest)<-ordertab$Order
  ordertabInterest=t(ordertabInterest[1:(nrow(ordertab)-1),])
  Z=c(rep(1,length(Group1)),rep(0,length(Group2)))
  X=NULL
  
  rdbRe <- rdb(ordertabInterest,Z,X=X)
  
  OrderwiseResult[[i]] <- rdbRe
}

########################################
## MsRDB analysis
########################################

MsRDBResult15 <- list()
for (i in 1:length(CompareGroup)) {
  GroupInterest <- CompareGroup[[i]]
  Group1<-meta$sample_name[meta$area==GroupInterest[1]&SampleList]
  Group2<-meta$sample_name[meta$area==GroupInterest[2]&SampleList]
  seqtabInterest=seqtab[c(Group1,Group2),]
  Z=c(rep(1,length(Group1)),rep(0,length(Group2)))
  X=NULL
  MsRDBResult15[[i]] <- msrdb(seqtabInterest, Z, X=X, seqdist=seqdist, knn = 15)
}

MsRDBResult15.cluster <- list()
for (i in 1:length(CompareGroup)) {
  GroupInterest <- CompareGroup[[i]]
  Group1<-meta$sample_name[meta$area==GroupInterest[1]&SampleList]
  Group2<-meta$sample_name[meta$area==GroupInterest[2]&SampleList]
  
  seqtabInterest=seqtab[c(Group1,Group2),]
  Z=c(rep(1,length(Group1)),rep(0,length(Group2)))
  X=NULL
  MsRDBResult15.cluster[[i]] <- msrdb.cluster(MsRDBResult15[[i]],15, 10,seqtabInterest) 
}

MsRDBResult15.cluster.tax <- list()
for (i in 1:length(CompareGroup)) {
  GroupInterest <- CompareGroup[[i]]
  Group1<-meta$sample_name[meta$area==GroupInterest[1]&SampleList]
  Group2<-meta$sample_name[meta$area==GroupInterest[2]&SampleList]
  
  seqtabInterest=seqtab[c(Group1,Group2),]
  Z=c(rep(1,length(Group1)),rep(0,length(Group2)))
  X=NULL
  members <- MsRDBResult15.cluster[[i]]$membership
  numcluster <- max(members)
  taxacluster <- list()
  for (j in 1:numcluster) {
    taxacluster[[j]]<-taxa[members==j,]
    #rownames(taxacluster[[j]]) <- NULL
  }
  MsRDBResult15.cluster.tax[[i]] <- taxacluster
}

########################################
## Overlapping Between Different Methods
########################################

display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

## Figure S7

## ASV Level 
k=1
x=list(ASVwise=row.names(taxa)[ASVwiseResult[[k]]],MsRDB=row.names(taxa)[MsRDBResult15.cluster[[k]]$sigindex])
display_venn(x, fill = c("#FF9933","#FFFF33"))

########################################
## Cluster Analysis
########################################
library(GGally)
library(intergraph)
library(colorspace)
library(igraph)

## Figure 6B-D
k=2
Ut=MsRDBResult15[[k]]$Sig
seqdistK=MsRDBResult15[[k]]$seqdistK
filter=MsRDBResult15[[k]]$filter
kk=15

sigindex <- which(Ut[filter]==TRUE)
sigK <- seqdistK$indexK[1:kk,sigindex]

adj <- matrix(0,nrow = sum(Ut), ncol = sum(Ut))
for (i in 1:sum(Ut)) {
  adj[i,] <- as.numeric(sigindex %in% sigK[,i])
  adj[i,i] <- 0
}
sigg <- graph_from_adjacency_matrix(adj, mode="undirected")
sigcluster <- clusters(sigg)

sigmembership=rep(0,length(Ut))
sigmembership[Ut]=sigcluster$membership
for (i in 1:sigcluster$no) {
  index <- which(sigmembership == i)
  seqtabcmb <- seqtabInterest[,index,drop=FALSE]
  clssum <- apply(seqtabcmb, 1, sum)
  if (sum(clssum)<=10) {
    Ut[index]=FALSE
  }
}
sigindex <- which(Ut[filter]==TRUE)
sigK <- seqdistK$indexK[1:kk,sigindex]

adj <- matrix(0,nrow = sum(Ut), ncol = sum(Ut))
for (i in 1:sum(Ut)) {
  adj[i,] <- as.numeric(sigindex %in% sigK[,i])
  adj[i,i] <- 0
}
sigg <- graph_from_adjacency_matrix(adj, mode="undirected")
sigcluster <- clusters(sigg)

V(sigg)$Cluster <- sigcluster$membership
cpl=rainbow_hcl(max(sigcluster$membership))
names(cpl)=1:max(sigcluster$membership)
ggnet2(sigg, color = "Cluster", label = sigcluster$membership, label.color = "white", palette = cpl,size =7)

########################################
## Detailed Cluster 
########################################
MsRDBResult15.cluster.tax[[1]]
MsRDBResult15.cluster.tax[[2]]
MsRDBResult15.cluster.tax[[3]]

k=2
j1=8
j2=9
j3=10
GroupInterest <- CompareGroup[[k]]
Group1<-meta$sample_name[meta$area==GroupInterest[1]&SampleList]
Group2<-meta$sample_name[meta$area==GroupInterest[2]&SampleList]

seqtabInterest=seqtab[c(Group1,Group2),]
Z=c(rep(GroupInterest[1],length(Group1)),rep(GroupInterest[2],length(Group2)))
X=NULL
members <- MsRDBResult15.cluster[[k]]$membership
df=cbind(Micrococcaceae=apply(seqtabInterest[,members==j1,drop=FALSE],1,sum),Geodermatophilaceae=apply(seqtabInterest[,members==j2,drop=FALSE],1,sum), Kineosporiaceae=apply(seqtabInterest[,members==j3,drop=FALSE],1,sum))
Benchmark <- apply(seqtabInterest[,!MsRDBResult15[[k]]$Sig],1,sum)
df<-data.frame(df/Benchmark)
df=cbind(df,Z)
df$Group <- as.factor(df$Z)
dff <- df %>% pivot_longer(1:3, names_to = "Family", values_to = "Value")
ggplot(dff, aes(x=Group, y=Value)) + 
  geom_boxplot() +
  geom_jitter(aes(shape=Group,color=Group), position=position_jitter(0.2)) +
  theme_bw() + 
  ylab("Relative Abundance w.r.t. Reference") +
  facet_grid(~ Family)


########################################
## Save Cluster
########################################

## Table S1

library(openxlsx)
sav.file <- "grapecluster.xlsx"
OUT <- createWorkbook()
for (i in 1:length(CompareGroup)) {
  GroupInterest <- CompareGroup[[i]]
  Cluster <- MsRDBResult15.cluster.tax[[i]]
  ClusterTable = cbind(cluster=rep("cluster 1",nrow(Cluster[[1]])),ASV=rownames(Cluster[[1]]),Cluster[[1]])
  for (j in 2:length(Cluster)) {
    TClusterTable = cbind(cluster=rep(paste("cluster",j),nrow(Cluster[[j]])),ASV=rownames(Cluster[[j]]),Cluster[[j]])
    ClusterTable = rbind(ClusterTable, TClusterTable)
  }
  sname<-paste(GroupInterest[1],GroupInterest[2],sep = "_")
  addWorksheet(OUT, sname)
  writeData(OUT, sheet = sname, x = ClusterTable)
}
saveWorkbook(OUT,sav.file)

