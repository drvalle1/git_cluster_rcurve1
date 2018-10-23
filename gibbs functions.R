acceptMH <- function(p0,p1,x0,x1,BLOCK){   #accept for M, M-H
  # if BLOCK, then accept as a block,
  # otherwise, accept individually
  
  nz           <- length(x0)  #no. to accept
  if(BLOCK) nz <- 1
  
  a    <- exp(p1 - p0)       #acceptance PR
  z    <- runif(nz,0,1)
  keep <- which(z < a)
  
  if(BLOCK & length(keep) > 0) x0 <- x1
  if(!BLOCK)                   x0[keep] <- x1[keep]           
  accept <- length(keep)        
  
  list(x = x0, accept = accept)
}
#-------------------------------
sample.omega=function(y,omega,nspp,nloc,jump,xmat,betas,cs,alpha,sig2){
  omega.old=omega
  omega.new=matrix(rnorm(nspp*nloc,mean=omega.old,sd=jump),nloc,nspp)
  
  #likelihood
  p1.old=dpois(y,exp(omega.old),log=T)
  p1.new=dpois(y,exp(omega.new),log=T)
  
  #prior
  media=xmat%*%betas
  alpha1=matrix(alpha,nloc,nspp,byrow=T)
  media1=alpha1+media[,cs]
  p2.old=dnorm(omega.old,mean=media1,sd=sqrt(sig2),log=T)
  p2.new=dnorm(omega.new,mean=media1,sd=sqrt(sig2),log=T)
  
  #MH algorithm
  k=acceptMH(p1.old+p2.old,p1.new+p2.new,omega.old,omega.new,F)  
  list(omega=k$x,accept=k$x!=omega.old)
}
#-----------------------------------------------
sample.alpha=function(nloc,sig2,xmat,betas,omega,cs,nspp){
  prec=(nloc/sig2)+(1/10)
  var1=1/prec
  media=xmat%*%betas
  media1=media[,cs]
  soma.err=colSums(omega-media1)
  pmedia=(1/sig2)*soma.err
  rnorm(nspp,mean=var1*pmedia,sd=sqrt(var1))
}
#-----------------------------------------------
sample.betas=function(ngroups,cs,nparam,xtx,t.xmat,sig2,alpha){
  i1=diag(1,nparam)
  alpha1=matrix(alpha,nloc,nspp,byrow=T)
  err=omega-alpha1
  betas=matrix(NA,nparam,ngroups)
  for (i in 1:ngroups){
    cond=cs==i
    nc=sum(cond)
    if (nc==0) betas[,i]=rmvnorm(1,mean=rep(0,nparam),sigma=i1)
    if (nc> 0) {
      prec=(nc/sig2)*xtx+i1
      var1=solve(prec)
      if (nc==1) soma.err=err[,cond]
      if (nc!=1) soma.err=rowSums(err[,cond])
      pmedia=(1/sig2)*t.xmat%*%soma.err
      betas[,i]=rmvnorm(1,mean=var1%*%pmedia,sigma=var1)
    }      

  }
  betas
}
#-----------------------------------------------
sample.sig2=function(nloc,nspp,omega,alpha,xmat,betas,cs){
  a1=((nloc*nspp)-1)/2
  alpha1=matrix(alpha,nloc,nspp,byrow=T)
  tmp=xmat%*%betas
  media=alpha1+tmp[,cs]
  err=omega-media
  b1=sum(err^2)/2
  1/rgamma(1,a1,b1)
}
#-----------------------------------------------
sample.cs=function(ngroups,omega,xmat,alpha,betas,sig2,theta){
  #calculate probabilities for each species for each group
  prob=matrix(NA,ngroups,nspp)
  alpha1=matrix(alpha,nloc,nspp,byrow=T)
  for (i in 1:ngroups){
    media=xmat%*%betas[,i]
    media1=alpha1+matrix(media,nloc,nspp)
    prob[i,]=colSums(dnorm(omega,mean=media1,sd=sqrt(sig2),log=T))+log(theta[i])
  }
  
  #make max=0
  max1=apply(prob,2,max)
  max2=matrix(max1,ngroups,nspp,byrow=T)
  tmp=prob-max2
  
  #normalize probabilities
  tmp1=exp(tmp)
  soma=matrix(colSums(tmp1),ngroups,nspp,byrow=T)
  prob=tmp1/soma
  
  #sample cs
  rmultinom1(prob=t(prob),randu=runif(nspp))+1
}
#-----------------------------------------------
sample.theta=function(cs,ngroups,gamma){
  n=rep(0,ngroups)
  tmp=table(cs)
  n[as.numeric(names(tmp))]=tmp
  
  v=theta=rep(NA,ngroups)
  prod=1
  for (i in 1:(ngroups-1)){
    n.greater.k=n[-(1:i)]
    v[i]=rbeta(1,n[i]+1,sum(n.greater.k)+gamma)
    theta[i]=v[i]*prod
    prod=prod*(1-v[i])
  }
  theta[ngroups]=prod
  theta
}
#----------------------------
print.adapt = function(accept1z,jump1z,accept.output){
  accept1=accept1z; jump1=jump1z; 
  
  for (k in 1:length(accept1)){
    z=accept1[[k]]/accept.output
    print(names(accept1)[k])
    print(mean(z)); print(mean(jump1[[k]]))
  }
  
  for (k in 1:length(jump1)){
    cond=(accept1[[k]]/accept.output)>0.4 & jump1[[k]]<10000
    jump1[[k]][cond] = jump1[[k]][cond]*2       
    cond=(accept1[[k]]/accept.output)<0.2 & jump1[[k]]>0.001
    jump1[[k]][cond] = jump1[[k]][cond]*0.5
    accept1[[k]][]=0
  }
  
  return(list(jump1=jump1,accept1=accept1))
}
#----------------------------
get.logl=function(y,omega,nspp,nloc,xmat,betas,cs,alpha,sig2){
  #likelihood
  p1=dpois(y,exp(omega),log=T)

  #prior
  media=xmat%*%betas
  alpha1=matrix(alpha,nloc,nspp,byrow=T)
  media1=alpha1+media[,cs]
  p2=dnorm(omega,mean=media1,sd=sqrt(sig2),log=T)

  #MH algorithm
  sum(p1)+sum(p2)
  
}