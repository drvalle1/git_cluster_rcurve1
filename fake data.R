rm(list=ls(all=TRUE))
set.seed(1)

nloc=1000
nspp=200
ngroup=5
nparam=6
xmat=matrix(rnorm(nparam*nloc),nloc,nparam)

alpha.true=alpha=2+rnorm(nspp)
betas=matrix(0,nparam,ngroup)
ind=sample(1:(nparam*ngroup),size=nparam*ngroup/3)
betas[ind]=sample(c(-0.5,0.5),size=length(ind),replace=T); betas
betas.true=betas
sig2.true=sig2=0.1
cs.true=cs=sample(1:ngroup,size=nspp,replace=T)

omega=matrix(NA,nloc,nspp)
for (i in 1:nspp){
  media=alpha[i]+xmat%*%betas[,cs[i]]
  omega[,i]=rnorm(nloc,mean=media,sd=sqrt(sig2))
}
omega.true=omega
y=matrix(rpois(nloc*nspp,exp(omega)),nloc,nspp)
range(colMeans(y==0))

setwd('U:\\GIT_models\\git_cluster_rcurve')
write.csv(y,'fake data.csv',row.names=F)
write.csv(xmat,'fake data xmat.csv',row.names=F)
