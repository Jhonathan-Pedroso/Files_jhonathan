setwd("/home/jhonathan/Results/Results_M")
load("M_500g_1.RData")
m1=M.f; rm(M.f)
load("M_500g_2.RData")
m2=M.f; rm(M.f)
load("M_500g_3.RData")
m3=M.f; rm(M.f)
load("M_500g_4.RData")
m4=M.f; rm(M.f)
load("M_500g_5.RData")
m5=M.f; rm(M.f)
load("M_500g_6.RData")
m6=M.f; rm(M.f)
load("M_500g_7.RData")
m7=M.f; rm(M.f)
load("M_500g_8.RData")
m8=M.f; rm(M.f)
load("M_500g_9.RData")
m9=M.f; rm(M.f)
load("M_500g_10.RData")
m10=M.f; rm(M.f)
M.f=cbind(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10); rm(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10)

test=0
for (i in 1:ncol(M.f)) {
  
  if (all(M.f[,i]==4) | all(M.f[,i]==3) | all(M.f[,i]==2) | all(M.f[,i]==1) | all(M.f[,i]==0)) {
    if (i==1) {
      test=as.matrix(i) 
    }
    if (i!=1) {
      test=cbind(test,i)
    }
  }
}
test=test[,-c(1)]
M.f=M.f[,-test]

A.meu <- function (K) {
  ploi=max(M)
  p=colSums(K)/(ploi*nrow(K))
  Kp=(matrix(0,nrow(K),ncol(K)))
  
  for (i in 1:nrow(K)) {
    for (j in 1:ncol(K)){
      Kp[i,j]=(K[i,j]-ploi*(p[j]))/sqrt(n*p[j]*(1-p[j]))
    }
  }
  A<-tcrossprod(Kp)/ncol(K)
  A<-A+diag(nrow(A))*1e-04
}

A=A.meu(M)

setwd("/home/jhonathan/Data/")
load("geno_sbdata.RData")

octaploide=t(octaploide)
M=octaploide[-c(1:2),]
test=0
index <- is.na(M)

n1=nrow(M); n2=ncol(M)
M=as.factor(M); M=as.numeric(M) 
M=matrix(M,n1,n2); M=M-1

p.m=colMeans(M, na.rm=TRUE)
hist(p.m)

for (i in 1:ncol(M)) {
  M[is.na(M[,i]),i] <-round(p.m[i],0)
}

ploi=8
p=colSums(M)/(ploi*nrow(M))
hist(p)

for (i in 1:ncol(M)) {
  if (i==1) {
    n=rep(0,ncol(M))
  }
  count=0
  if (any(M[,i]==8)) {
    count=1
  } 
  if (any(M[,i]==7)) {
    count=count+1
  } 
  if (any(M[,i]==6)) {
    count=count+1
  } 
  if (any(M[,i]==5)) {
    count=count+1
  } 
  if (any(M[,i]==4)) {
    count=count+1
  } 
  if (any(M[,i]==3)) {
    count=count+1
  } 
  if (any(M[,i]==2)) {
    count=count+1
  } 
  if (any(M[,i]==1)) {
    count=count+1
  } 
  if (any(M[,i]==0)) {
    count=count+1
  } 
  n[i]=count  
}
hist(n)

A=A.meu(M)

A.vit= function (K) {
  ploi=max(K)
  freq=colSums(K)/(ploi*nrow(K))
  Fa=t(matrix(freq,ncol(K),nrow(K)))
  
  W=K-ploi*Fa
  cor=ploi*sum(freq*(1-freq))
  Av=W%*%t(W)/cor
  return (Av)
}

A=A.vit(M)

h2 <- 0.9
ngen=200
gen=rnorm(ngen,mean=0,sd=c(runif(ngen,0,100)))
loc=sample(1:ncol(M), ngen, replace=F)
g=matrix(0,ncol(M),1)
g[loc,] <-gen; u=M%*%g
y.all <- 10000 + u + rnorm(nrow(M),mean=0,sd=sqrt((1-h2)/h2*var(u)))
cov(y.all,u)/var(y.all); sqrt((1-h2)/h2*var(u)); var(u)

A=A.vit(M)
library(rrBLUP)
res.v=mixed.solve(y.all,K=A)
cor(res.v$u,u); cov(y.all,res.v$u)/var(y.all)
res.v$Vu
var(u)

A=A.meu(M)
res.m=mixed.solve(y.all,K=A)
cor(res.m$u,u); cov(y.all,res.m$u)/var(y.all)
res.v$Vu
var(u)
