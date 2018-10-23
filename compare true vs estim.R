compare=function(true,estim){
  rango=range(c(true,estim))
  plot(true,estim,xlim=rango,ylim=rango)
  lines(rango,rango,col='red')
}

compare(omega,omega.true)

compare(alpha,alpha.true)

plot(vec.sig2[1:i],type='l')
abline(h=sig2.true,col='red')

k=data.frame(estim.cs=cs,true.cs=cs.true)
table(k)

ind=c(10,9,5,4,2)

compare(betas[,ind],betas.true)

plot(theta,type='h')
