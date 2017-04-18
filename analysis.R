rm(list=ls())

################
### Analysis
################

# Import data
galaxy.data <- read.table("galaxy.txt")

str(galaxy.data)
colnames(galaxy.data)=c("galaxy")

str(galaxy.data)

attach(galaxy.data)

hist(galaxy,nclass=50,prob=T)
lines(density(galaxy))

# Metropolis-Hastings algorithm

source("jeffmixt_fun.R")
library(gtools)
library(coda)
library(truncnorm)
library(mvtnorm)
library(QRM)

nsim <- 10^5
burn <- 0.5*nsim


ysam <- galaxy

hist(ysam,nclass=50,prob=T)
lines(density(ysam))

obj <- post_Jeffreys_mixture( ysam=ysam , nsim=nsim , burn=burn , Kw=10)

param.post <- obj$param.post

par(mfrow=c(2,2))
plot(param.post[,1],type="l",ylim=c(0,1))
for(i in 2:Kw)
{
	lines(param.post[,i],col=i)
}
plot(param.post[,(Kw+1)],type="l",ylim=c(0,150))
for(i in (Kw+2):(2*Kw))
{
	lines(param.post[,i],col=(i-Kw))
}
plot(param.post[,(2*Kw+1)],type="l",ylim=c(0,10))
for(i in (2*Kw+2):(3*Kw))
{
	lines(param.post[,i],col=(i-(2*Kw)))
}
plot(param.post[,31],type="l")

param.post.mean <- apply(param.post[1:nsim,],2,mean,na.rm=T) 
param.post.sd <- apply(param.post[1:nsim,],2,sd,na.rm=T) 

ic <- apply(param.post[burn:nsim,],2,quantile,prob=c(0.025,0.975))

y.seq <- seq( min(ysam), max(ysam), 0.01 )

fy.estim <- c()
for(i in 1:length(y.seq))
{
	fy.estim[i] <- targKcomp( y.seq[i], ww=param.post.mean[1:Kw] , 
						mmean=param.post.mean[(Kw+1):(2*Kw)] ,
						ssd=param.post.mean[(2*Kw+1):(3*Kw)]
						 ) 
						 
	print(i)
}

fy.estim.025 <- c()
for(i in 1:length(y.seq))
{
	fy.estim.025[i] <- targKcomp( y.seq[i], ww=ic[1,1:Kw] , 
						mmean=ic[1,(Kw+1):(2*Kw)] ,
						ssd=ic[1,(2*Kw+1):(3*Kw)]
						 ) 
						 
	print(i)
}

fy.estim.975 <- c()
for(i in 1:length(y.seq))
{
	fy.estim.975[i] <- targKcomp( y.seq[i], ww=ic[2,1:Kw] , 
						mmean=ic[2,(Kw+1):(2*Kw)] ,
						ssd=ic[2,(2*Kw+1):(3*Kw)]
						 ) 
						 
	print(i)
}


color.transparent <- adjustcolor("lightblue",alpha.f=1)
par(mfrow=c(1,1))
hist(ysam,nclass=25,prob=T,ylim=c(0,0.25),col="grey",
		main="Galaxy",xlab="velocity of galaxy (1000km/s)")
		
for(i in 1:length(y.seq))
{
	lines( c(y.seq[i],y.seq[i]), c(fy.estim.025[i],fy.estim.975[i]) , 
			col = color.transparent )
}
lines( y.seq, fy.estim, col="red" , lwd=2 )

round(cbind(ic[1,],param.post.mean,ic[2,]), digits=3)


