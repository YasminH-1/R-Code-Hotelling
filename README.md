# R-Code-Hotelling
x=read.csv("Yelmy.csv")
plot(x)
x=as.matrix(x)
x
f=x[1==x[,1],]
f
f=f[,-1]
f
m=x[1!=x[,1],]
m
m=m[,-1]
m
#hotelling with outliers assuming non-equal variance - copied the commands from my class's lecture slides
#m is matrix of males, f is matrix of females, a is number of males, b is number of females, p is number of variables
#applying hotelling on dataset with outliers and equal variance
ht2=function(x,y){
  n=dim(x)[1]; m=dim(y)[1]; p=dim(x)[2]
  xcov=cov(x); ycov=cov(y)
  SP=(n-1)*xcov+(m-1)*ycov; SP=SP/(n+m-2)
  xcentre=colMeans(x); ycentre=colMeans(y)
  d=xcentre-ycentre
  T2=t(d)%*%solve(SP)%*%d
  T2=T2*n*m/(n+m)
  F=T2*(n+m-p-1)/(p*(n+m-2))
  pv=1-pf(F,p,n+m-p-1)
  list(MaleCentre=xcentre, FemaleCentre=ycentre, MaleCovariance=xcov,
       FemaleCovariance=ycov, PooledVariance=SP, Tsquared=T2, F=F, df=c(p,n+m-p-1),PValue=pv)
  
}
ht2(m,f)
#applying hotelling on dataset with outliers and unequal variance
ht=function(x,y){
  n=dim(x)[1]; m=dim(y)[1]; p=dim(x)[2]
  xcov=cov(x); ycov=cov(y)
  S=(xcov/n)+(ycov/m);
  xcentre=colMeans(x); ycentre=colMeans(y)
  d=xcentre-ycentre
  T2=t(d)%*%solve(S)%*%d
  F=T2
  pv=1-pf(F,p,n+m-p-1)
  list(MaleCentre=xcentre, FemaleCentre=ycentre, MaleCovariance=xcov,FemaleCovariance=ycov, Variance=S, Tsquared=T2, F=F, df=c(p,n+m-p-1),PValue=pv)
  
}
ht(m,f)
library(robustX)
require(robustX)
n=dim(x)[1]
p=dim(x)[2]
