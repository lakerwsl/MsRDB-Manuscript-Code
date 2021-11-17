source("./MsRDB/Algorithm.R")
source("./MsRDB/ASVwise.R")
source("./MsRDB/MultiATE.R")

library(biomformat)
library(tidyverse)
library(phyloseq)
library(ggplot2)
library(dada2)
library(DECIPHER)
library(ANCOMBC)
library(ALDEx2)
library(microbiome)

########################################
## Read Data
########################################

x = read_biom("./Data/Vangay_Cell_2018/63800/all.biom")
seqtab<-as(biom_data(x), "matrix")
seqtab<-t(seqtab)

#taxa <- assignTaxonomy(seqtab, "../Tax/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
taxa<-read.table("./Data/Vangay_Cell_2018/63800/taxa.txt", sep="\t")

seq <- getSequences(seqtab)
seqlist <- DNAStringSet(seq)
seqdist <- DistanceMatrix(seqlist)

meta = read.delim("./Data/Vangay_Cell_2018/12080_20210508-071948.txt")
meta <- meta[meta$sample_name %in% row.names(seqtab),]
rownames(meta) <- meta$sample_name

########################################
## Preprocessing
########################################

selectedasv<-taxa$Kingdom %in% c("Bacteria")
seqtab<-seqtab[,selectedasv]
taxa<-taxa[selectedasv,]
seqdist <- seqdist[selectedasv,selectedasv]

########################################
## Group into Different Rank
########################################

sampleinterest <- meta$sample_name[meta$sample_group %in% c("HmongThai","Control","Hmong2nd")]

##Genus

genustab <-as.data.frame(t(seqtab[sampleinterest,]))
genustab <-cbind(genustab,taxa[row.names(genustab),c("Genus"),drop=FALSE])
genustab <- genustab %>% group_by(Genus) %>% summarise_all(sum)
rowsum=apply(genustab[,-1],1,sum)
genustab <- genustab[rowsum!=0,]

genusunassignratio<-genustab[nrow(genustab),2:ncol(genustab)]/apply(genustab[,2:ncol(genustab)], 2, sum)
hist(as.numeric(genusunassignratio),100)
mean(as.numeric(genusunassignratio), na.rm=TRUE)

ggplot(data.frame(x=t(genusunassignratio),y=rep(1,length(genusunassignratio))), aes(x=x)) +
  geom_histogram(position="identity", binwidth=0.01, alpha=0.5, color='red',fill='red') +
  theme_minimal() +
  xlab('Proportion of NA Genus') +
  ylab('Number of Samples') +
  xlim(-0.01, 0.5) 

##Class

classtab <-as.data.frame(t(seqtab[sampleinterest,]))
classtab <-cbind(classtab,taxa[row.names(classtab),c("Class"),drop=FALSE])
classtab <- classtab %>% group_by(Class) %>% summarise_all(sum)
rowsum=apply(classtab[,-1],1,sum)
classtab <- classtab[rowsum!=0,]

classunassignratio<-classtab[nrow(classtab),2:ncol(classtab)]/apply(classtab[,2:ncol(classtab)], 2, sum)
hist(as.numeric(classunassignratio),100)
mean(as.numeric(classunassignratio), na.rm=TRUE)

ggplot(data.frame(x=t(classunassignratio),y=rep(1,length(classunassignratio))), aes(x=x)) +
  geom_histogram(position="identity", binwidth=0.01, alpha=0.5, color='red',fill='red') +
  theme_minimal() +
  xlab('Proportion of NA Class') +
  ylab('Number of Samples') +
  xlim(-0.01, 0.5) 

########################################
## Covariate Balancing
########################################

cbsample<- c("HmongThai","Control","Hmong2nd")
cbmeta<-meta[meta$sample_group %in% cbsample, ]

## Figure S3

## age distribution
ggplot(cbmeta, aes(x=sample_group, y=age)) + 
  geom_boxplot() + ylab("Age") + xlab("Group") +
  geom_jitter(aes(shape=sample_group,color=sample_group), position=position_jitter(0.2)) +
  theme_bw() 

## bmi distribution
ggplot(cbmeta, aes(x=sample_group, y=bmi)) + 
  geom_boxplot() + ylab("BMI") + xlab("Group") +
  geom_jitter(aes(shape=sample_group,color=sample_group), position=position_jitter(0.2)) +
  theme_bw() 

########################################
## ASV-wise analysis
########################################

CompareGroup <- list(c("HmongThai","Control"),c("Hmong2nd","Control"),c("Hmong2nd","HmongThai"))
CovariateInterest <- c("age","bmi")

ASVwiseResult <- list()
for (i in 1:length(CompareGroup)) {
  GroupInterest <- CompareGroup[[i]]
  Group1<-meta$sample_name[meta$sample_group==GroupInterest[1]]
  Group2<-meta$sample_name[meta$sample_group==GroupInterest[2]]
  
  seqtabInterest=seqtab[c(Group1,Group2),]
  Z=c(rep(1,length(Group1)),rep(0,length(Group2)))
  X=meta[c(Group1,Group2),CovariateInterest]
  
  rdbRe <- rdb(seqtabInterest,Z,X=X)
  
  otu_mat = t(seqtabInterest)
  meta_mat = cbind(group=Z,X)
  tax_mat = as.matrix(taxa)
  OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
  META = sample_data(meta_mat)
  TAX = tax_table(tax_mat)
  physeq = phyloseq(OTU, META, TAX)
  ancombcRe=ancombc(phyloseq=physeq,formula="group+age+bmi",p_adj_method="bonferroni",alpha = 0.1, zero_cut=1)

  Xmatrix <- model.matrix(as.formula("~ group+age+bmi"), meta_mat)
  aclr <- aldex.clr(t(seqtabInterest), Xmatrix)
  aldexRe <- aldex.glm(aclr)
  
  ASVwiseResult[[i]] <- list(rdbRe=rdbRe,ancombcRe=ancombcRe,aldexRe=aldexRe)
}

for (i in 1:length(CompareGroup)) {
  TempRe <- ASVwiseResult[[i]]
  print(sum(TempRe$rdbRe))
  print(sum(TempRe$ancombcRe$res$diff_abn$group))
  print(sum(TempRe$aldexRe$`model.group Pr(>|t|).BH`<0.1))
}




########################################
## Taxon-wise analysis
########################################

##Genus

GenuswiseResult <- list()
for (i in 1:length(CompareGroup)) {
  GroupInterest <- CompareGroup[[i]]
  Group1<-meta$sample_name[meta$sample_group==GroupInterest[1]]
  Group2<-meta$sample_name[meta$sample_group==GroupInterest[2]]
  
  genustabInterest=genustab[,c(Group1,Group2)]
  genustabInterest=as.matrix(genustabInterest)
  rownames(genustabInterest)<-genustab$Genus
  genustabInterest=t(genustabInterest[1:(nrow(genustab)-1),])
  Z=c(rep(1,length(Group1)),rep(0,length(Group2)))
  X=meta[c(Group1,Group2),CovariateInterest]
  
  rdbRe <- rdb(genustabInterest,Z,X=X)
  
  seqtabInterest=seqtab[c(Group1,Group2),]
  otu_mat = t(seqtabInterest)
  meta_mat = cbind(group=Z,X)
  tax_mat = as.matrix(taxa)
  OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
  META = sample_data(meta_mat)
  TAX = tax_table(tax_mat)
  physeq = phyloseq(OTU, META, TAX)
  physeq = aggregate_taxa(physeq, "Genus")
  ancombcRe=ancombc(phyloseq=physeq,formula="group+age+bmi",p_adj_method="bonferroni",alpha = 0.1, zero_cut=1)

  Xmatrix <- model.matrix(as.formula("~ group+age+bmi"), meta_mat)
  aclr <- aldex.clr(t(genustabInterest), Xmatrix)
  aldexRe <- aldex.glm(aclr)
  
  GenuswiseResult[[i]] <- list(rdbRe=rdbRe,ancombcRe=ancombcRe,aldexRe=aldexRe)
}

for (i in 1:length(CompareGroup)) {
  TempRe <- GenuswiseResult[[i]]
  print(sum(TempRe$rdbRe))
  print(sum(TempRe$ancombcRe$res$diff_abn$group))
  print(sum(TempRe$aldexRe$`model.group Pr(>|t|).BH`<0.1))
}


##Class

ClasswiseResult <- list()
for (i in 1:length(CompareGroup)) {
  GroupInterest <- CompareGroup[[i]]
  Group1<-meta$sample_name[meta$sample_group==GroupInterest[1]]
  Group2<-meta$sample_name[meta$sample_group==GroupInterest[2]]
  
  classtabInterest=classtab[,c(Group1,Group2)]
  classtabInterest=as.matrix(classtabInterest)
  rownames(classtabInterest)<-classtab$Class
  classtabInterest=t(classtabInterest[1:(nrow(classtab)-1),])
  Z=c(rep(1,length(Group1)),rep(0,length(Group2)))
  X=meta[c(Group1,Group2),CovariateInterest]
  
  rdbRe <- rdb(classtabInterest,Z,X=X)
  
  seqtabInterest=seqtab[c(Group1,Group2),]
  otu_mat = t(seqtabInterest)
  meta_mat = cbind(group=Z,X)
  tax_mat = as.matrix(taxa)
  OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
  META = sample_data(meta_mat)
  TAX = tax_table(tax_mat)
  physeq = phyloseq(OTU, META, TAX)
  physeq = aggregate_taxa(physeq, "Class")
  ancombcRe=ancombc(phyloseq=physeq,formula="group+age+bmi",p_adj_method="bonferroni",alpha = 0.1, zero_cut=1)
  
  Xmatrix <- model.matrix(as.formula("~ group+age+bmi"), meta_mat)
  aclr <- aldex.clr(t(classtabInterest), Xmatrix)
  aldexRe <- aldex.glm(aclr)
  
  ClasswiseResult[[i]] <- list(rdbRe=rdbRe,ancombcRe=ancombcRe,aldexRe=aldexRe)
}

for (i in 1:length(CompareGroup)) {
  TempRe <- ClasswiseResult[[i]]
  print(sum(TempRe$rdbRe))
  print(sum(TempRe$ancombcRe$res$diff_abn$group))
  print(sum(TempRe$aldexRe$`model.group Pr(>|t|).BH`<0.1))
}



########################################
## MsRDB analysis
########################################

MsRDBResult20 <- list()
for (i in 1:length(CompareGroup)) {
  GroupInterest <- CompareGroup[[i]]
  Group1<-meta$sample_name[meta$sample_group==GroupInterest[1]]
  Group2<-meta$sample_name[meta$sample_group==GroupInterest[2]]
  
  seqtabInterest=seqtab[c(Group1,Group2),]
  Z=c(rep(1,length(Group1)),rep(0,length(Group2)))
  X=meta[c(Group1,Group2),CovariateInterest]
  
  MsRDBResult20[[i]] <- msrdb(seqtabInterest, Z, X=X, seqdist=seqdist, knn = 20)
}

for (i in 1:length(CompareGroup)) {
  TempRe <- MsRDBResult20[[i]]
  print(sum(TempRe$Sig))
}

########################################
## Overlapping Between Different Methods
########################################
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

## ASV Level 
comparisons=c("","_CH2","_HH2")
k=3

filename=paste("ASVVenn",comparisons[k],".png",sep = "")
png(filename, width = 400, height = 400,pointsize = 16)
x=list(ASVwise=row.names(taxa)[ASVwiseResult[[k]]$rdbRe],MsRDB=row.names(taxa)[MsRDBResult20[[k]]$Sig])
display_venn(x, fill = c("#FF9933","#FFFF33"))
dev.off()

filename=paste("ASVVenn_ancombc",comparisons[k],".png",sep = "")
png(filename, width = 400, height = 400,pointsize = 16)
ancombcASV=row.names(ASVwiseResult[[k]]$ancombcRe$res$diff_abn)[ASVwiseResult[[k]]$ancombcRe$res$diff_abn$group]
x=list(ASVwise=ancombcASV,MsRDB=row.names(taxa)[MsRDBResult20[[k]]$Sig])
display_venn(x, fill = c("#FF9933","#FFFF33"))
dev.off()

filename=paste("ASVVenn_aldex",comparisons[k],".png",sep = "")
png(filename, width = 400, height = 400,pointsize = 16)
aldexASV=row.names(ASVwiseResult[[k]]$aldexRe)[ASVwiseResult[[k]]$aldexRe$`model.group Pr(>|t|).BH`<0.1]
x=list(ASVwise=aldexASV,MsRDB=row.names(taxa)[MsRDBResult20[[k]]$Sig])
display_venn(x, fill = c("#FF9933","#FFFF33"))
dev.off()

## Genus Level
filename=paste("GenusVenn",comparisons[k],".png",sep = "")
png(filename, width = 400, height = 400,pointsize = 16)
x=list(ASVwise=unique(taxa$Genus[ASVwiseResult[[k]]$rdbRe]), Genuswise=genustab$Genus[GenuswiseResult[[k]]$rdbRe], MsRDB=unique(taxa$Genus[MsRDBResult20[[k]]$Sig]))
x$MsRDB[is.na(x$MsRDB)]="NA"
x$ASVwise[is.na(x$ASVwise)]="NA"
display_venn(x, fill = c("#FF3333","#FF9933", "#FFFF33"))
dev.off()

filename=paste("GenusVenn_ancombc",comparisons[k],".png",sep = "")
png(filename, width = 400, height = 400,pointsize = 16)
ancombcGenus=row.names(GenuswiseResult[[k]]$ancombcRe$res$diff_abn)[GenuswiseResult[[k]]$ancombcRe$res$diff_abn$group]
x=list(ASVwise=unique(taxa$Genus[row.names(taxa) %in% ancombcASV]), Genuswise=ancombcGenus, MsRDB=unique(taxa$Genus[MsRDBResult20[[k]]$Sig]))
x$MsRDB[is.na(x$MsRDB)]="NA"
x$ASVwise[is.na(x$ASVwise)]="NA"
display_venn(x, fill = c("#FF3333","#FF9933", "#FFFF33"))
dev.off()

filename=paste("GenusVenn_aldex",comparisons[k],".png",sep = "")
png(filename, width = 400, height = 400,pointsize = 16)
aldexGenus=row.names(GenuswiseResult[[k]]$aldexRe)[GenuswiseResult[[k]]$aldexRe$`model.group Pr(>|t|).BH`<0.1]
x=list(ASVwise=unique(taxa$Genus[row.names(taxa) %in% aldexASV]), Genuswise=aldexGenus, MsRDB=unique(taxa$Genus[MsRDBResult20[[k]]$Sig]))
x$MsRDB[is.na(x$MsRDB)]="NA"
x$ASVwise[is.na(x$ASVwise)]="NA"
display_venn(x, fill = c("#FF3333","#FF9933", "#FFFF33"))
dev.off()

## Class Level
filename=paste("ClassVenn",comparisons[k],".png",sep = "")
png(filename, width = 400, height = 400,pointsize = 16)
x=list(ASVwise=unique(taxa$Class[ASVwiseResult[[k]]$rdbRe]), Genuswise=unique(taxa$Class[taxa$Genus %in% genustab$Genus[GenuswiseResult[[k]]$rdbRe]]), Classwise=classtab$Class[ClasswiseResult[[k]]$rdbRe], MsRDB=unique(taxa$Class[MsRDBResult20[[k]]$Sig]))
x$MsRDB[is.na(x$MsRDB)]="NA"
x$ASVwise[is.na(x$ASVwise)]="NA"
x$Genuswise[is.na(x$Genuswise)]="NA"
x$Classwise[is.na(x$Classwise)]="NA"
display_venn(x, fill = c("#FF3333","#FF9933", "#FFFF33","#99FF33"))
dev.off()


filename=paste("ClassVenn_ancombc",comparisons[k],".png",sep = "")
png(filename, width = 400, height = 400,pointsize = 16)
ancombcClass=row.names(ClasswiseResult[[k]]$ancombcRe$res$diff_abn)[ClasswiseResult[[k]]$ancombcRe$res$diff_abn$group]
x=list(ASVwise=unique(taxa$Class[row.names(taxa) %in% ancombcASV]), Genuswise=unique(taxa$Class[taxa$Genus %in% ancombcGenus]), Classwise=ancombcClass, MsRDB=unique(taxa$Class[MsRDBResult20[[k]]$Sig]))
x$MsRDB[is.na(x$MsRDB)]="NA"
x$ASVwise[is.na(x$ASVwise)]="NA"
x$Genuswise[is.na(x$Genuswise)]="NA"
x$Classwise[is.na(x$Classwise)]="NA"
display_venn(x, fill = c("#FF3333","#FF9933", "#FFFF33","#99FF33"))
dev.off()


filename=paste("ClassVenn_aldex",comparisons[k],".png",sep = "")
png(filename, width = 400, height = 400,pointsize = 16)
aldexClass=row.names(ClasswiseResult[[k]]$aldexRe)[ClasswiseResult[[k]]$aldexRe$`model.group Pr(>|t|).BH`<0.1]
x=list(ASVwise=unique(taxa$Class[row.names(taxa) %in% aldexASV]), Genuswise=unique(taxa$Class[taxa$Genus %in% aldexGenus]), Classwise=aldexClass, MsRDB=unique(taxa$Class[MsRDBResult20[[k]]$Sig]))
x$MsRDB[is.na(x$MsRDB)]="NA"
x$ASVwise[is.na(x$ASVwise)]="NA"
x$Genuswise[is.na(x$Genuswise)]="NA"
x$Classwise[is.na(x$Classwise)]="NA"
display_venn(x, fill = c("#FF3333","#FF9933", "#FFFF33","#99FF33"))
dev.off()

########################################
## Detailed Difference Between Different Methods
########################################

k=1
GroupInterest <- CompareGroup[[k]]
Group1<-meta$sample_name[meta$sample_group==GroupInterest[1]]
Group2<-meta$sample_name[meta$sample_group==GroupInterest[2]]

classtabInterest=classtab[,c(Group1,Group2)]
classtabInterest=as.matrix(classtabInterest)
rownames(classtabInterest)<-classtab$Class
classtabInterest=t(classtabInterest[1:(nrow(classtab)-1),])

genustabInterest=genustab[,c(Group1,Group2)]
genustabInterest=as.matrix(genustabInterest)
rownames(genustabInterest)<-genustab$Genus
genustabInterest=t(genustabInterest[1:(nrow(genustab)-1),])

seqtabInterest=seqtab[c(Group1,Group2),]
Z=c(rep(GroupInterest[1],length(Group1)),rep(GroupInterest[2],length(Group2)))

classx=list(ASVwise=unique(taxa$Class[ASVwiseResult[[k]]$rdbRe]), Genuswise=unique(taxa$Class[taxa$Genus %in% genustab$Genus[GenuswiseResult[[k]]$rdbRe]]), Classwise=classtab$Class[ClasswiseResult[[k]]$rdbRe], MsRDB=unique(taxa$Class[MsRDBResult20[[k]]$Sig]))
genusx=list(ASVwise=unique(taxa$Genus[ASVwiseResult[[k]]$rdbRe]), Genuswise=genustab$Genus[GenuswiseResult[[k]]$rdbRe], MsRDB=unique(taxa$Genus[MsRDBResult20[[k]]$Sig]))
asvx=list(ASVwise=row.names(taxa)[ASVwiseResult[[k]]$rdbRe],MsRDB=row.names(taxa)[MsRDBResult20[[k]]$Sig])

ClassDiff1 <- classx$Genuswise[!(classx$Genuswise %in% classx$Classwise)]
ClassDiff1 <- ClassDiff1[2]
GenusClassDiff1 <- unique(taxa[taxa$Class %in% ClassDiff1,'Genus'])
GenusClassDiff1 <- genusx$Genuswise[genusx$Genuswise %in% GenusClassDiff1]
df=cbind(classtabInterest[,ClassDiff1,drop=FALSE],genustabInterest[,GenusClassDiff1])
Benchmark <- apply(seqtabInterest[,!MsRDBResult20[[k]]$Sig],1,sum)
df<-data.frame(df/Benchmark)
df=cbind(df,Z)
df$Group <- as.factor(df$Z)


## Figure 5C

dff <- df %>% pivot_longer(2:6, names_to = "Genus", values_to = "Relative Abundance w.r.t. Reference")
ggplot(dff, aes(x=Group, y=`Relative Abundance w.r.t. Reference`)) + 
  geom_boxplot() +
  geom_jitter(aes(shape=Group,color=Group), position=position_jitter(0.2)) +
  theme_bw() + 
  facet_grid(~ Genus)



## Figure 5A

df$SignificantGenus=apply(df[,2:6],1,sum)
df$Rest=df$Bacilli-df$SignificantGenus
dff <- df %>% pivot_longer(c(1,9,10), names_to = "Type", values_to = "Relative Abundance w.r.t. Reference")
ggplot(dff, aes(x=Group, y=`Relative Abundance w.r.t. Reference`)) + 
  geom_boxplot() +
  geom_jitter(aes(shape=Group,color=Group), position=position_jitter(0.2)) +
  theme_bw() + 
  facet_grid(~ Type)








GenusDiff <- genusx$MsRDB[!(genusx$MsRDB %in% genusx$Genuswise)]
GenusDiff <- GenusDiff[3]
ASVGenusDiff <- row.names(taxa)[taxa$Genus %in% GenusDiff]
ASVGenusDiff <- asvx$MsRDB[asvx$MsRDB %in% ASVGenusDiff]
df=cbind(genustabInterest[,GenusDiff,drop=FALSE],seqtabInterest[,ASVGenusDiff])
Benchmark <- apply(seqtabInterest[,!MsRDBResult20[[k]]$Sig],1,sum)
df<-data.frame(df/Benchmark)
df=cbind(df,Z)
df$Group <- as.factor(df$Z)
df$Combine=apply(df[,2:3],1,sum)

## Figure 5B

ggplot(df, aes(x=Group, y=Combine)) + 
  geom_boxplot() +
  geom_jitter(aes(shape=Group,color=Group), position=position_jitter(0.2)) +
  theme_bw() + 
  ylab("Relative Abundance w.r.t. Reference")
mean(df$Combine[df$Z==GroupInterest[1]]>0)
mean(df$Combine[df$Z==GroupInterest[2]]>0)



GenusDiff <- genusx$MsRDB[!(genusx$MsRDB %in% genusx$Genuswise)]
GenusDiff <- GenusDiff[4]
ASVGenusDiff <- row.names(taxa)[taxa$Genus %in% GenusDiff]
ASVGenusDiff <- asvx$MsRDB[asvx$MsRDB %in% ASVGenusDiff]
df=cbind(genustabInterest[,GenusDiff,drop=FALSE],seqtabInterest[,ASVGenusDiff])
Benchmark <- apply(seqtabInterest[,!MsRDBResult20[[k]]$Sig],1,sum)
df<-data.frame(df/Benchmark)
df=cbind(df,Z)
df$Group <- as.factor(df$Z)
df$Combine=apply(df[,2:5],1,sum)
df$Diff=df$X.Ruminococcus..torques.group-df$Combine

mean(df$Combine[df$Z==GroupInterest[1]]>0)
mean(df$Combine[df$Z==GroupInterest[2]]>0)


dff <- df %>% pivot_longer(c(1,8,9), names_to = "Type", values_to = "Relative Abundance w.r.t. Reference")
dff$Type[dff$Type=="X.Ruminococcus..torques.group"]="Ruminococcus torques"
dff$Type[dff$Type=="Combine"]="SignificantASV"
dff$Type[dff$Type=="Diff"]="Rest"
dff$Type_f = factor(dff$Type, levels=c("Ruminococcus torques","Rest","SignificantASV"))

## Figure S4A

ggplot(dff, aes(x=Group, y=`Relative Abundance w.r.t. Reference`)) + 
  geom_boxplot() +
  geom_jitter(aes(shape=Group,color=Group), position=position_jitter(0.2)) +
  theme_bw() + 
  facet_grid(~ Type_f)







NAASV <- asvx$MsRDB[is.na(taxa[asvx$MsRDB,'Genus'])]
df=seqtabInterest[,NAASV]
Benchmark <- apply(seqtabInterest[,!MsRDBResult20[[k]]$Sig],1,sum)
df<-data.frame(df/Benchmark)
df=cbind(df,Z)
df$Group <- as.factor(df$Z)
df$Combine=apply(df[,1:(ncol(df)-2)],1,sum)

## Figure S4B

ggplot(df, aes(x=Group, y=Combine)) + 
  geom_boxplot() +
  geom_jitter(aes(shape=Group,color=Group), position=position_jitter(0.2)) +
  theme_bw() + 
  ylab("Relative Abundance w.r.t. Reference")

NAFamily<-taxa[NAASV,'Family']
NAFamily[is.na(NAFamily)]='NA'
NAOrder<-taxa[NAASV,'Order']
NAOrder[is.na(NAOrder)]='NA'
NAClass<-taxa[NAASV,'Class']
NAClass[is.na(NAClass)]='NA'
NAPhylum<-taxa[NAASV,'Phylum']
NAPhylum[is.na(NAPhylum)]='NA'
df <- cbind(NAFamily,NAOrder,NAClass,NAPhylum)
df<-data.frame(df)

## Figure 5D

ggplot(df) +
  geom_bar(aes(y = NAFamily), color='#FFFF33', fill='#FF9933') +
  theme_bw() + 
  ylab("Family") +
  xlab("ASV Count")
ggplot(df) +
  geom_bar(aes(y = NAOrder), color='#FFFF33', fill='#FF9933') +
  theme_bw() + 
  ylab("Order") +
  xlab("ASV Count")






























  