library(dada2)
library(kmer)
library(DECIPHER)
library(igraph)

msrdb <- function(seqtab, Z, X=NULL, seqdist=NULL, disttype="align",knn=NULL, alpha=0.1, fdr=TRUE) {
  filter <- apply(seqtab, 2, sum)>0
  seqtab.filter<-seqtab[,filter]
  
  n=nrow(seqtab.filter)
  d=ncol(seqtab.filter)
  
  if (is.null(seqdist)) {
    seq <- getSequences(seqtab.filter)
    if (disttype=="align") {
      seqlist <- DNAStringSet(seq)
      seqdist <- DistanceMatrix(seqlist)
    } else if (disttype=="kmer") {
      seqlist <- lapply(as.list(seq), function(x){unlist(strsplit(x, split = ""))})
      seqdist<-kdistance(seqlist,k=5)
      seqdist<-as.matrix(seqdist)
    }
  } else {
    seqdist=seqdist[filter,filter]
  }
  
  if(is.vector(X)){
    X<- matrix(X, ncol = 1, nrow = length(X))
  }
  
  if (is.null(knn)) {
    knn=max(50,ceiling(d/200))
  }
  seqdistK <- knnfinding(seqdist,knn)
  
  seqtabP=seqtab.filter
  for (i in 1:n) {
    seqtabP[i,]=seqtab.filter[i,]/sum(seqtab.filter[i,])
  }
  
  wt=wtfinding(seqtabP, Z, X, seqdistK)
  Ut=wtrdb(seqtabP, Z, X, wt, alpha, fdr)
  ComUt<-rep(FALSE,ncol(seqtab))
  ComUt[filter]=Ut
  
  rdbresult <- list(Sig=ComUt, seqdistK=seqdistK, filter=filter)
  
  return(rdbresult)
}

knnfinding <- function(seqdist, K) {
  d=ncol(seqdist)
  distK <- matrix(0, nrow = K, ncol = d)
  indexK <- matrix(0, nrow = K, ncol = d)
  
  for (j in 1:d) {
    tempdist <- seqdist[j,]
    tempdist[j] <- -1
    orderdist <- order(tempdist)
    indexK[,j] <- orderdist[1:K]
    distK[,j] <- tempdist[indexK[,j]]
    distK[1,j] <- 0
  }
  seqdistK <- list(indexK=indexK,distK=distK,K=K)
  return(seqdistK)
}

wtfinding <- function(seqtab, Z, X, seqdistK) {
  treat=Z==1
  mtreat=sum(treat)
  control=Z==0
  mcontrol=sum(control)
  
  d=ncol(seqtab)
  h=10*sqrt(log(d))
  H=2*sqrt(log(d))
  K=seqdistK$K
  indexK=seqdistK$indexK
  wtK=matrix(0, nrow = K, ncol = d)
  wtseqtab=seqtab
  wtmeanvar=matrix(0,nrow = 5,ncol = d)
  if (is.null(X)) {
    Ptreat=seqtab[treat,]
    Pcontrol=seqtab[control,]
    wtmeanvar[1,]=apply(Ptreat, 2, mean)
    wtmeanvar[2,]=apply(Pcontrol, 2, mean)
    wtmeanvar[3,]=apply(Ptreat, 2, var)/mtreat
    wtmeanvar[5,]=apply(Pcontrol, 2, var)/mcontrol
  } else {
    wtmeanvar=MultipleATE(seqtab, Z, X)
  }
  DA=(wtmeanvar[1,]-wtmeanvar[2,])/(wtmeanvar[1,]+wtmeanvar[2,])
  DABenchmark=DA
  DAsdBenchmark=DA
  DAFinal=DA
  
  numiter <- 10
  halfnumiter = ceiling(numiter/2)
  incK <- K^{1/numiter}
  curK<-incK
  KK=min(ceiling(curK),K)
  stopstatus<-rep(FALSE,d)
  for (t in 1:numiter)
  {
    DAratio <- wtmeanvar[1,]/(wtmeanvar[1,]+wtmeanvar[2,])
    DAsd=sqrt(wtmeanvar[3,]*(1-DAratio)*(1-DAratio)-2*DAratio*(1-DAratio)*wtmeanvar[4,]+wtmeanvar[5,]*DAratio*DAratio)
    DAsd=DAsd/sqrt(wtmeanvar[1,]+wtmeanvar[2,])
    for (j in 1:d) {
      tempindex=indexK[1:KK,j]
      DAdiff=DA[j]-DA[tempindex]
      if (DAsd[j]==0) {
        wtK[1:KK,j]=rep(0,KK)
        wtK[DAdiff==0,j]=1
      } else {
        DAdiff=DAdiff/DAsd[j]
        wtK[1:KK,j]=kernel(DAdiff/h)
      }
      wtseqtab[,j] <- seqtab[,tempindex] %*% wtK[1:KK,j]
    }
    if (is.null(X)) {
      wtPtreat=wtseqtab[treat,]
      wtPcontrol=wtseqtab[control,]
      wtmeanvar[1,]=apply(wtPtreat, 2, mean)
      wtmeanvar[2,]=apply(wtPcontrol, 2, mean)
      wtmeanvar[3,]=apply(wtPtreat, 2, var)/mtreat
      wtmeanvar[5,]=apply(wtPcontrol, 2, var)/mcontrol
    } else {
      wtmeanvar=MultipleATE(wtseqtab, Z, X)
    }
    DA=(wtmeanvar[1,]-wtmeanvar[2,])/(wtmeanvar[1,]+wtmeanvar[2,])
    
    if (t==halfnumiter) {
      DABenchmark=DA
      #DavarBenchmark=(wtmeanvar[3,]-2*wtmeanvar[4,]+wtmeanvar[5,])/(wtmeanvar[1,]+wtmeanvar[2,])^2
      DAratio <- wtmeanvar[1,]/(wtmeanvar[1,]+wtmeanvar[2,])
      DAsd=sqrt(wtmeanvar[3,]*(1-DAratio)*(1-DAratio)-2*DAratio*(1-DAratio)*wtmeanvar[4,]+wtmeanvar[5,]*DAratio*DAratio)
      DAsdBenchmark=DAsd/sqrt(wtmeanvar[1,j]+wtmeanvar[2,j])
      DAFinal=DA
      stopstatus[DAsdBenchmark==0]=TRUE
    } else if (t>halfnumiter) {
      DAdiff=abs(DA-DABenchmark)/DAsdBenchmark
      DAdiff[DAsdBenchmark==0]=2*H
      index=DAdiff>H
      stopstatus=stopstatus|index
      DAFinal[!stopstatus]=DA[!stopstatus]
      DA[stopstatus]=DAFinal[stopstatus]
    }
    curK<-curK*incK
    KK=min(ceiling(curK),K)
  }
  wt <- list(wtK=wtK,indexK=indexK,K=K)
  return(wt)
}
  
wtrdb <- function(seqtab, Z, X, wt, alpha, fdr) {
  treat=Z==1
  mtreat=sum(treat)
  control=Z==0
  mcontrol=sum(control)
  d=ncol(seqtab)
  M=sqrt(2*log(d)/d)
  D=sqrt(2*log(d)-2*log(alpha))
  Dpm=D+0.2*M
  K=wt$K
  indexK=wt$indexK
  wtK=wt$wtK
  
  wtseqtab=seqtab
  wtmeanvar=matrix(0,nrow = 5,ncol = d)
  for (j in 1:d) {
    tempindex=indexK[1:K,j]
    wtseqtab[,j] <- seqtab[,tempindex] %*% wtK[,j]
  }
  if (is.null(X)) {
    wtPtreat=wtseqtab[treat,]
    wtPcontrol=wtseqtab[control,]
    wtmeanvar[1,]=apply(wtPtreat, 2, mean)
    wtmeanvar[2,]=apply(wtPcontrol, 2, mean)
    wtmeanvar[3,]=apply(wtPtreat, 2, var)/mtreat
    wtmeanvar[5,]=apply(wtPcontrol, 2, var)/mcontrol
  } else {
    wtmeanvar=MultipleATE(wtseqtab, Z, X)
  }
  
  Pmean=matrix(0,nrow = 2,ncol = d)
  Ptreat=seqtab[treat,]
  Pcontrol=seqtab[control,]
  Pmean[1,]=apply(Ptreat, 2, mean)
  Pmean[2,]=apply(Pcontrol, 2, mean)
  
  if (fdr == FALSE)
  {
    Vt <- wtrdb.core(wtmeanvar, M, D, Dpm, Pmean)
  } else {
    Dlist<-seq(0,D,length.out=100)
    pFDR<-rep(0,length(Dlist))
    for (j in 1:length(Dlist))
    {
      DD = Dlist[j]
      tVt <- wtrdb.core(wtmeanvar, M, DD, DD+0.2*M, Pmean)
      R=max(1,sum(!tVt))
      R0=2*d*pnorm(DD,lower.tail =FALSE)
      #R0=d*pchisq(DD*DD,df=2,lower.tail =FALSE)
      pFDR[j]=R0/R
    }
    if (pFDR[100]<alpha)
    {
      FThre <- Dlist[min(which(pFDR<=alpha))]
    } else {
      FThre <- Dlist[100]
    }
    Vt <- wtrdb.core(wtmeanvar, M, FThre, FThre+0.2*M, Pmean)
  }
  return(!Vt)
}

resultsort <- function(cluster, seq) {
  numcluster <- max(cluster$membership)
  result <- list()
  for (i in 1:numcluster) {
    indexcluster <- cluster$sigindex[which(cluster$membership==i)]
    seqcluster <- seq[indexcluster]
    result[[i]] <- seqcluster
  }
  return(result)
}

wtrdb.core <- function(meanvar, M, D, Dpm, Pmean) {
  d=ncol(meanvar)
  Vt=rep(TRUE,d)
  while (sum(Vt)>0) {
    Rtreat=sum(Pmean[1,Vt])
    Rcontrol=sum(Pmean[2,Vt])
    if (Rtreat<=0 && Rcontrol>0) {
      tstat=meanvar[2,]/sqrt(meanvar[5,])
    } else if (Rtreat>0 && Rcontrol<=0) {
      tstat=meanvar[1,]/sqrt(meanvar[3,])
    } else if (Rtreat>0 && Rcontrol>0) {
      tstat=(meanvar[1,]/Rtreat-meanvar[2,]/Rcontrol)/sqrt(meanvar[3,]/Rtreat/Rtreat-2*meanvar[4,]/Rtreat/Rcontrol+meanvar[5,]/Rcontrol/Rcontrol)
    } else {
      break
    }
    Mtstat=median(tstat[Vt])
    if (Mtstat>M) {
      Wt=Vt&(tstat<(-D))
    } else if (Mtstat<(-M)) {
      Wt=Vt&(tstat>D)
    } else {
      Wt=Vt&(abs(tstat)>Dpm)
    }
    if (sum(Wt)==0){
      break
    }
    Vt[Wt]=FALSE
  }
  return(Vt)
}

kernel <- function(x) {
  x=x*x
  y=exp(-2*x)
  return(y)
} 

msrdb.cluster <- function(rdbresult,kk, c_cut, seqtab) {
  Ut=rdbresult$Sig
  seqdistK=rdbresult$seqdistK
  filter=rdbresult$filter
  
  if (kk>seqdistK$K) kk=seqdistK$K
  
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
    seqtabcmb <- seqtab[,index,drop=FALSE]
    clssum <- apply(seqtabcmb, 1, sum)
    if (sum(clssum)<=c_cut) {
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
  
  sigmembership=rep(0,length(Ut))
  sigmembership[Ut]=sigcluster$membership

  recluster <- list(sigindex=Ut, membership=sigmembership, csize=sigcluster$csize)
  
  return(recluster)
}

