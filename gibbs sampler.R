# rm(list=ls(all=TRUE))
set.seed(2)
library('mvtnorm')
library('Rcpp')

setwd('U:\\GIT_models\\git_cluster_rcurve')
source('gibbs functions.R')
sourceCpp('aux1.cpp')
y=data.matrix(read.csv('fake data.csv',as.is=T))
xmat=data.matrix(read.csv('fake data xmat.csv',as.is=T))
t.xmat=t(xmat)
xtx=t.xmat%*%xmat

#useful stuff
nloc=nrow(y)
nspp=ncol(y)
nparam=ncol(xmat)
ngroups=10
gamma=0.01

#initial parameter values
betas=matrix(0,nparam,ngroups)
alpha=rep(0,nspp)
omega=log(y)
omega[y==0]=0
sig2=10
cs=sample(1:ngroups,size=nspp,replace=T)
theta=rep(1/ngroups,ngroups)

#MH stuff
jump1=list(omega=matrix(0.1,nloc,nspp))
accept1=list(omega=matrix(0,nloc,nspp))
accept.output=50; nadapt=1000

#gibbs stuff
ngibbs=10000
vec.betas=matrix(NA,ngibbs,nparam*ngroups)
vec.alpha=matrix(NA,ngibbs,nspp)
vec.sig2=matrix(NA,ngibbs,1)
vec.theta=matrix(NA,ngibbs,ngroups)
vec.logl=matrix(NA,ngibbs,1)

for (i in 1:ngibbs){
  print(i)
  tmp=sample.omega(y=y,omega=omega,nspp=nspp,nloc=nloc,jump=jump1$omega,
                   xmat=xmat,betas=betas,cs=cs,alpha=alpha,sig2=sig2)
  omega=tmp$omega
  accept1$omega=accept1$omega+tmp$accept
  # omega=omega.true
    
  alpha=sample.alpha(nloc=nloc,sig2=sig2,xmat=xmat,betas=betas,omega=omega,cs=cs,nspp=nspp)
  # alpha=alpha.true
  
  betas=sample.betas(ngroups=ngroups,cs=cs,nparam=nparam,xtx=xtx,t.xmat=t.xmat,sig2=sig2,alpha=alpha)
  # betas=betas.true
  
  sig2=sample.sig2(nloc=nloc,nspp=nspp,omega=omega,alpha=alpha,xmat=xmat,betas=betas,cs=cs)
  
  cs=sample.cs(ngroups=ngroups,omega=omega,xmat=xmat,alpha=alpha,betas=betas,sig2=sig2,theta=theta)
  # cs=cs.true
  
  theta=sample.theta(cs=cs,ngroups=ngroups,gamma=gamma)
  # theta=rep(1/ngroups,ngroups)
  
  logl=get.logl(y=y,omega=omega,nspp=nspp,nloc=nloc,xmat=xmat,betas=betas,cs=cs,alpha=alpha,sig2=sig2)
  
  #adapt MH
  if (i%%accept.output==0 & i<nadapt){
    k=print.adapt(accept1z=accept1,jump1z=jump1,accept.output=accept.output)
    accept1=k$accept1
    jump1=k$jump1
  }
  
  #store results
  vec.betas[i,]=betas
  vec.alpha[i,]=alpha
  vec.sig2[i]=sig2
  vec.theta[i,]=theta
  vec.logl[i]=logl
}
plot(vec.logl[1:i],type='l')