## The MultipleATE is a generalized version of ATE in ATE package, which is tailed for our purpose. 

MultipleATE <- function(seqtab, Z, X) {
  
  weights <- sampleweightestimation(Z, X)
  means <- getmean(seqtab, Z, weights)
  variance <- getvariance(seqtab, Z, weights, means)
  
  return(rbind(means, variance))
}

getvariance <- function(seqtab, Z, weights, means) {
  d=ncol(seqtab)
  N<- length(Z)
  n1<- sum(Z)
  lam<- as.vector(weights$lam.p)
  be<- as.vector(weights$lam.q)
  FUNu<-  function(x) c(1,x)
  umat <- t(apply(X, 1, FUNu))
  K<- length(lam)
  theta=0
  rho2<-function(x){dd.cr.rho(x,theta=theta)}
  
  variance=matrix(0, nrow = 3, ncol = d)
  for (j in 1:d) {
    A<- matrix(0,ncol = 2*K, nrow = 2*K)
    C<- matrix(0, ncol = 2*K, nrow = 2)
    meat<- matrix(0, ncol = 2*K+2, nrow = 2*K+2)
    for(i in 1:N){
      A[(1:K),(1:K)]<- A[(1:K),(1:K)] + 
        Z[i]*rho2(crossprod(lam, umat[i,]))*tcrossprod(umat[i,])
      A[((K+1):(2*K)),((K+1):(2*K))]<- A[((K+1):(2*K)),((K+1):(2*K))] + 
        (1-Z[i])*rho2(crossprod(be, umat[i,]))*tcrossprod(umat[i,])
      
      C[1,(1:K)]<- C[1,(1:K)] + Z[i]*rho2(crossprod(lam, umat[i,]))*seqtab[i,j]*umat[i,]
      C[2,((K+1):(2*K))]<- C[2,((K+1):(2*K))] +
        (1-Z[i])*rho2(crossprod(be, umat[i,]))*seqtab[i,j]*umat[i,]
      
      meat<- meat +  tcrossprod(get.gk(X[i,], seqtab[i,j],Z[i], lam,be, means[,j]))
    }
    A<- A/N
    C<- C/N
    meat<- meat/N
    
    bread<-  matrix(0, nrow = 2*K+2, ncol = 2*K+2)
    bread[1:(2*K),1:(2*K)]<- A
    bread[2*K+1:2, ]<- cbind(C,diag(c(-1,-1)))
    bread<- solve(bread)
    cov.mat<-(bread%*%meat%*%t(bread))/N
    cov.mat<- cov.mat[-(1:(2*K) ),-(1:(2*K))]
    variance[1,j]=cov.mat[1,1]
    variance[2,j]=cov.mat[1,2]
    variance[3,j]=cov.mat[2,2]
  }
  return(variance)
}

getmean <- function(seqtab, Z, weights) {
  d=ncol(seqtab)
  means=matrix(0, nrow = 2, ncol = d)
  for (j in 1:d) {
    means[1,j]=sum((weights$weights.p*seqtab[,j])[Z==1])
    means[2,j]=sum((weights$weights.q*seqtab[,j])[Z==0])
  }
  return(means)
}

sampleweightestimation <- function(Z, X) {
  K<- ncol(X)+1
  J<- 2
  
  p.hat<- newton.r(Z, X)
  q.hat<- newton.r(1-Z, X)
  
  w.p<- p.hat$weights
  w.q<- q.hat$weights
  
  lam.p<- tail(p.hat$res,1)
  lam.q<- tail(q.hat$res,1)
  
  res.l<- list("lam.p" = lam.p, "lam.q" = lam.q, "weights.p" = w.p,
               "weights.q" = w.q)
  return(res.l)
  
}

newton.r<-function (Z, X) {
  max.iter = 100
  theta=0
  tol = 0.001
  backtrack = TRUE
  bt.a = 0.3 
  bt.b = 0.5
  
  rho<-function(x){cr.rho(x,theta=theta)}
  rho1<-function(x){d.cr.rho(x,theta=theta)}
  rho2<-function(x){dd.cr.rho(x,theta=theta)}
  FUNu<-  function(x) c(1,x)
  ini<- rep(0, ncol(X)+1)
  N <- length(Z)
  
  res <- matrix(ini, nrow = 1)
  umat <- t(apply(X, 1, FUNu))
  ubar <- apply(umat, 2, mean)
  dervs<- c(0)
  
  for (i in 2:max.iter) {
    x.val <- res[i-1, ]
    objectiveValue<- obj(res[i-1,],umat,ubar,Z,rho)
    
    if(objectiveValue <= -1e30){
      warning("The objective function is unbounded, a different rho() function 
              might be needed.")
      lam.hat <- as.vector(tail(res, 1))
      l.t.u <- apply(umat, 1, crossprod, lam.hat)
      nom <- rep(0, N)
      nom[Z == 1] <- rho1(l.t.u)[Z == 1]/N
      weights <- nom
      return(list("res" = res, "weights" = weights, "conv" = FALSE))
    }
    del.f <- derv.obj(x.val, Z, rho1, u = umat, ubar = ubar)
    H <- derv2.obj(x.val, Z, rho2, u = umat)
    del.x <- -solve(H, del.f)
    nabla.f <- del.f
    dervs<- c(dervs,sum(nabla.f^2))
    
    step <- 1
    if (backtrack) {
      step <- backtrack(bt.a, bt.b, x.val, del.x, nabla.f, 
                        obj, umat, ubar, Z, rho)
    }
    
    res <- rbind(res, x.val + step * del.x)
    
    if (sum(nabla.f^2) < tol) {
      lam.hat <- as.vector(tail(res, 1))
      l.t.u <- apply(umat, 1, crossprod, lam.hat)
      nom <- rep(0, N)
      nom[Z == 1] <- rho1(l.t.u)[Z == 1]/N
      weights <- nom
      return(list("res" = res, "weights" = weights, "conv" = TRUE))
    }
  }
  
  warning("The algorithm did not converge")
  lam.hat <- as.vector(tail(res, 1))
  l.t.u <- apply(umat, 1, crossprod, lam.hat)
  nom <- rep(0, N)
  nom[Z == 1] <- rho1(l.t.u)[Z == 1]/N
  weights <- nom
  return(list("res" = res, "weights" = weights, "conv" = FALSE))
}

obj<- function(lam,u,ubar,Z,rho){
  lam.t.u<- apply(u,1,crossprod,lam)
  #plot(lam.t.u)
  lam.t.ubar<- crossprod(lam,ubar)
  
  -mean(Z*rho(lam.t.u),na.rm = TRUE) + lam.t.ubar
}

cr.rho<-function (v, theta) 
{
  v<-as.vector(v)
  if (theta == 0) {
    ans <- -exp(-v)
  }
  else if (theta == -1){
    ans <- suppressWarnings(log(1+v))
    ind <- is.na(ans)
    ans[ind] <- -Inf
  }
  else  {
    a <- -(1 - theta * v)^(1 + 1/theta)
    ans <- a/(theta + 1)
    ind <- is.na(ans)
    ans[ind] <- -Inf
  }
  return(ans)
}


#The first and second derivatives of the CR family functions

d.cr.rho<-  function (v, theta) 
{
  v<-as.vector(v)
  if (theta == 0) {
    ans <- exp(-v)
  }
  else if (theta == -1){
    ans <- 1/(1+v)
    ind <- is.na(ans)
    ans[ind] <- -Inf
  }
  else {
    ans <- (1 - theta * v)^(1/theta)
    ind <- is.na(ans)
    ans[ind] <- -Inf
  }
  return(ans)
}

dd.cr.rho<-   function (v, theta) 
{
  v<-as.vector(v)
  if (theta == 0) {
    a <- -exp(-v)
  }
  else if (theta == -1){
    a <- -1/(1+v)^2
    ind <- is.na(a)
    a[ind] <- -Inf
  }
  else {
    a <- -(1 - theta * v)^(1/theta - 1)
    ind <- is.na(a)
    a[ind] <- -Inf
  }
  a
}

derv.obj<- function(lam,Z,rho1,u,ubar){
  lam.t.u <- apply(u,1,crossprod,lam)
  lam.t.ubar<- crossprod(lam,ubar)
  
  #get rho1(lam^Tu)*R/pi
  v<- numeric(length(Z))
  v[Z==1]<- rho1(lam.t.u)[Z==1]
  #get final answer
  -apply(apply(u,2,"*",v),2,mean,na.rm = TRUE)+ubar
}

derv2.obj<- function(lam,Z,rho2,u){
  lam.t.u<- apply(u,1,crossprod,lam)
  
  #get rho(lam^Tu)
  v<- numeric(length(Z))
  v[Z==1]<- rho2(lam.t.u)[Z==1]
  
  #Get matrices for hessian
  mats<- sapply(1:nrow(u), function(i) tcrossprod(u[i,]),simplify = "array")
  mats2<- apply(mats,c(1,2),"*",v)
  #get final answer
  -apply(mats2,c(2,3),mean,na.rm = TRUE)
}

backtrack<- function(alpha,beta,x.val,del.x,nabla.f,obj,u, ubar,Z,rho){
  u<- u
  step<- 1
  
  l.t.u<- apply(u,1,crossprod,x.val)
  f.x<- obj(x.val,u,ubar, Z, rho)
  
  df.t.dx<- crossprod(nabla.f, del.x)
  #print(f.x)
  #print(df.t.dx)
  
  while(obj(x.val+step*del.x, u,ubar,Z,rho) > f.x + alpha*step*df.t.dx ){
    step<- beta*step
  }
  return(step)
}

get.gk<- function(X.vec,Y.scaler,Z.int,lam,beta,means.ind){
  theta=0
  rho1<-function(x){d.cr.rho(x,theta=theta)}
  FUNu<-  function(x) c(1,x)
  
  lam<- as.vector(lam)
  be<- as.vector(beta)
  uk <- FUNu(as.numeric(X.vec))
  
  gk1<- Z.int*rho1(crossprod(lam, uk))*uk - uk
  gk2<- (1-Z.int)*rho1(crossprod(be, uk))*uk - uk
  gk3<- Z.int*rho1(crossprod(lam,uk))*Y.scaler- means.ind[1]
  gk4<- (1-Z.int)*rho1(crossprod(be, uk))*Y.scaler - means.ind[2]
  
  gk<- c(gk1,gk2,gk3,gk4) 
  gk
}





