rdb <- function(P, Z, X=NULL, alpha=0.1, fdr=TRUE)
{
  if (nrow(P)!=length(Z))
    stop("Please make sure the number of rows in P equal to the length of Z")
  if (!all(Z %in% c(0, 1)))
    stop("Please make sure Z only contain 0 and 1")
  for (i in 1:nrow(P))
  {
    P[i,]=P[i,]/sum(P[i,])
  }
  if (!is.null(X))
  {
    if (is.vector(X)) X<-matrix(X,ncol=1)
    if (nrow(P)!=nrow(X))
      stop("Please make sure the number of rows in P equal to number of rows in X")
  }
  
  treat=Z==1
  mtreat=sum(treat)
  control=Z==0
  mcontrol=sum(control)
  
  filter=apply(P, 2, sum)>0
  P.filter=P[,filter]
  
  d=ncol(P.filter)
  M=sqrt(2*log(d)/d)
  D=sqrt(2*log(d)-2*log(alpha))
  Dpm=D+0.2*M
  
  meanvar=matrix(0,nrow = 5,ncol = d)
  if (is.null(X)) {
    Ptreat=P.filter[treat,]
    Pcontrol=P.filter[control,]
    meanvar[1,]=apply(Ptreat, 2, mean)
    meanvar[2,]=apply(Pcontrol, 2, mean)
    meanvar[3,]=apply(Ptreat, 2, var)/mtreat
    meanvar[5,]=apply(Pcontrol, 2, var)/mcontrol
  } else {
    meanvar=MultipleATE(P.filter, Z, X)
  }
  
  if (fdr == FALSE)
  {
    Vt <- rdb.core(meanvar, M, D, Dpm)
  } else {
    Dlist<-seq(0,D,length.out=100)
    pFDR<-rep(0,length(Dlist))
    for (j in 1:length(Dlist))
    {
      DD = Dlist[j]
      tVt <- rdb.core(meanvar, M, DD, DD+0.2*M)
      R=max(1,sum(!tVt))
      R0=2*d*pnorm(DD,lower.tail =FALSE)
      pFDR[j]=R0/R
    }
    if (pFDR[100]<alpha)
    {
      FThre <- Dlist[min(which(pFDR<=alpha))]
    } else {
      FThre <- Dlist[100]
    }
    Vt <- rdb.core(meanvar, M, FThre, FThre+0.2*M)
  }
  
  ComUt<-rep(FALSE,ncol(P))
  ComUt[filter]=!Vt
  
  return(ComUt)
}

rdb.core <- function(meanvar, M, D, Dpm)
{
  d=ncol(meanvar)
  Vt=rep(TRUE,d)
  while (sum(Vt)>0) {
    Rtreat=sum(meanvar[1,Vt])
    Rcontrol=sum(meanvar[2,Vt])
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

rdb.knn <- function(seqtab, Z, X=NULL, seqdist=NULL, disttype="align",knn=NULL, alpha=0.1, fdr=TRUE) {
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
  
  if (is.null(knn)) {
    knn=max(50,ceiling(d/200))
  }
  seqdistK <- knnfinding(seqdist,knn)
  
  seqtabP=seqtab.filter
  for (i in 1:n) {
    seqtabP[i,]=seqtab.filter[i,]/sum(seqtab.filter[i,])
  }
  
  wtK=matrix(1, nrow = knn, ncol = d)
  wt <- list(wtK=wtK,indexK=seqdistK$indexK,K=knn)
  
  Ut=wtrdb(seqtabP, Z, X, wt, alpha, fdr)
  ComUt<-rep(FALSE,ncol(seqtab))
  ComUt[filter]=Ut

  return(ComUt)
}

