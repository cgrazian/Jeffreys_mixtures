# List of functions

# Target distribution: Gaussian mixture model

targKcomp <- function(x, ww, mmean, ssd){
	sum(ww*dnorm(x, mmean, ssd))
}

# Monte Carlo integration

intMC <- function(ysam, elem){
	 ysam <- as.matrix(ysam)
	temp <- apply(ysam, 1, elem)
	return( mean( temp ) ) 
}

# Fisher information matrix for the weights

FishInfoKcompMC <- function(ww, mmean, ssd, MCsam){

	Kw <- length(ww)

	mm <- matrix(0,nrow=(Kw-1),ncol=(Kw-1))
	for(i in 1:(Kw-1)){
		for(j in i:(Kw-1)){

			elem <- function(x){
					(dnorm(x, mmean[i], ssd[i]) - dnorm(x, mmean[Kw], ssd[Kw]))*
						(dnorm(x, mmean[j], ssd[j]) - dnorm(x, mmean[Kw], ssd[Kw]))/
							(targKcomp(x, ww, mmean, ssd))^2
			}
			mm[i,j] <- intMC(MCsam, elem)
 
		}
	}

	mm1 <- mm + t( mm ) 
	diag( mm1 ) <- diag( mm )

	return(mm1)

}

# Jeffreys prior for the mixture weights

jeffKcompMC <- function(ww, mmean, ssd, MCsam){
	sqrt( det( FishInfoKcompMC( ww=ww, mmean=mmean, ssd=ssd, MCsam=MCsam ) ) )
}

# Functions for the Pareto distribution

dpareto <- function(x, a=0.5, b=1) a*b^a/x^(a+1)
ppareto <- function(x, a=0.5, b=1) (x > b)*(1-(b/x)^a)
qpareto <- function(u, a=0.5, b=1) b/(1-u)^(1/a)
rpareto <- function(n, a=0.5, b=1) qpareto(runif(n),a,b)

# Likelihood for a Gaussian mixture model

llik.mix.aug <- function(x, z, ww, mu, sigma)
{
	n <- length(x)
	llik.par <- c()
	for(i in 1:n)
	{
		llik.par[i] <- sum( log(ww[z[i,]==1])+
							
				log(dnorm(x[i],mu,sigma)[z[i,]==1] ) ) 
	}
	return(sum(llik.par))
}

# Mixture prior for sigma

prior.sigma <- function(sigma,a=1,zeta0)
{
	0.5*dunif(sigma,0,zeta0) + 0.5*dpareto(sigma,a=1,b=zeta0)
}

# Posterior distribution for a Gaussian mixture model

lpost.mix.augMC <- function(x,z,ww,mu,sigma,mu0,zeta0)
{
	Kw <- length(ww)
	
	# Likelihood
	lik <- llik.mix.aug(x,z,ww,mu,sigma)
	
	# First step prior for sigma and mu
	lprior.mu <- sum(log(dnorm(mu,mu0,zeta0)))
	lprior.sigma <- sum(log(prior.sigma(sigma,a=1,zeta0)))
	
	# First step prior for ww
	nsamMC <- 500
	ncompMC <- round (max(ww[1]*nsamMC , 2),digits=0)
	MCysam <- rnorm(ncompMC,mu[1],sigma[1])
	for(u in 2:length(ww)){
		ncompMC <- round (max(ww[u]*nsamMC , 2),digits=0)
		MCysam <- c(MCysam,rnorm(ncompMC,mu[u],sigma[u]))
	}
	lprior.w <- log(jeffKcompMC(ww,mu,sigma,MCysam))
	
	# Secondo step prior
	lprior.mu0zeta0 <- -log(zeta0)
	
	return(lik+lprior.mu+lprior.sigma+lprior.w+lprior.mu0zeta0)
}

# MCMC simulation from the posterior distribution of a Gaussian mixture model with the 
#	Jeffreys prior for the mixture weights

post_Jeffreys_mixture <- function( ysam , nsim=10^5 , burn=0.5*nsim , Kw)
{

	# ysam = sample
	# nsim = number of Monte Carlo simulations
	# burn = burnin
	# Kw = number of components in the mixtures

	nsam <- length(ysam)

	# Initialization
	post.den <- -Inf
	count<-0
	while(post.den==-Inf)
	{
		count<- count+1
		ww.inits <- rdirichlet(1,rep(1,Kw))
		xi.inits <- rnorm(1,0,10)
		K.inits <- rgamma(1,0.1,0.1)			

		bug <- c()	
		z.mat <- matrix(NA,nrow=nsam,ncol=Kw)
		for(i in 1:nsam)
		{
			temp <- 1
			if(sum(temp!=0)==0){
				bug[i] <- 0 	
			} else {
				bug[i] <- 1
				z.prob <- ww.inits
				z.mat[i,] <- rmultinom(1,1,z.prob)
			}
		}

		z.freq.inits <- apply(z.mat,2,mean)

		mu.mean <- c()
		sd.mean <- c()
		for(temp in 1:Kw)
		{
			mu.mean[temp] <- ifelse( z.freq.inits[temp]!=0 , 
							mean(ysam[z.mat[,temp]==1]) , 0 )
			sd.mean[temp] <- ifelse( sum(z.mat[,temp])>1 , 
							sd(ysam[z.mat[,temp]==1]) , 1 )
		}

		m.inits <- rnorm(Kw,mu.mean,sd.mean)
		sd.inits <- rtruncnorm(1,a=0,b=Inf,mean=sd.mean,sd=sd.mean)

		if(sum(bug)==nsam){
			post.den <- lpost.mix.augMC(ysam,z.mat,ww.inits,
								m.inits,sd.inits,xi.inits,K.inits)

			den.lik <- llik.mix.aug(ysam,z.mat,ww.inits,m.inits,sd.inits)
			z.old <- z.mat	
		} else {
			post.den <- -Inf
		}
		print(paste(“inits”,count,sep=“-“)
	}

	lik.val <- c(den.lik)
	post.val <- c(post.den)

	param.post <- matrix(NA,nrow=nsim,ncol=(3*Kw+2))
	param.post[1,1:Kw] <- ww.inits
	param.post[1,(Kw+1):(2*Kw)] <- m.inits
	param.post[1,(2*Kw+1):(3*Kw)] <- sd.inits
	param.post[1,(3*Kw+1)] <- xi.inits
	param.post[1,(3*Kw+2)] <- K.inits

	z.freq <- matrix(NA,nrow=nsim,ncol=Kw)
	z.freq[1,] <- apply(z.old,2,mean)

	acc.ww <- 0
	acc.m <- 0
	acc.sd <- 0
	acc.xi <- 0
	acc.K <- 0

	### Metropolis within Gibbs
	###		accept/reject at each step
	###		with the full conditional of a standard
	###		Gibbs Sampling as proposals
	###		
	###

	for(count in 2:nsim)
	{

		if(round(count/10^4)==(count/10^4)){
			print(paste(“sim”,count,sep=“-“))
		}
 
		### Proposals
	
		ww.old <- param.post[count-1,1:Kw]
		m.old <- param.post[count-1,(Kw+1):(2*Kw)]
		sd.old <- param.post[count-1,(2*Kw+1):(3*Kw)]
		xi.old <- param.post[count-1,(3*Kw+1)]
		K.old <- param.post[count-1,(3*Kw+2)]
		bug <- "old.values"
	
		### Update of z
		q.z.num <- 0
		q.z.den <- 0
		for(i in 1:nsam)
		{
			z.prob <- ww.old * dnorm( ysam[i] , m.old , sd.old)
			z.mat[i,] <- rmultinom(1,1,z.prob)
			q.z.num <- q.z.num + dmultinom(z.old[i,],1,z.prob)
			q.z.den <- q.z.den + dmultinom(z.mat[i,],1,z.prob)				
		} 
		
		post.num <- lpost.mix.augMC(ysam,z.mat,ww.old,
								m.old,sd.old,xi.old,K.old)	
							
		num <- post.num + q.z.num
		den <- post.den + q.z.den
	
		u <- runif(1)
	
		if(u<exp(num-den)){
			z.new <- z.mat
			post.den <- post.num
		} else {
			z.new <- z.old
			post.den <- post.den
		}
	
		n <- apply(z.new,2,sum)
		sx <- c()
		for(j in 1:Kw)
		{
			sx[j]<-sum(ysam[z.new[,j]==1])
		}

		### Update of the weights
	
		ww.prop <- rdirichlet(1,1+n)
		
		q.ww.num <- log( ddirichlet(ww.old, 1+n ) )
		q.ww.den <- log( ddirichlet(ww.prop, 1+n ) )
		
		post.num <- lpost.mix.augMC(ysam,z.new,ww.prop,
								m.old,sd.old,xi.old,K.old)	
	
		num <- post.num + q.ww.num
		den <- post.den + q.ww.den

		u <- runif(1)
		if(u<exp(num-den)){
			ww.new <- ww.prop
			post.den <- post.num
			acc.ww <- acc.ww + 1
		} else {
			ww.new <- ww.old
			post.den <- post.den
		}

		### Update of the means
	
		mean.q <- ifelse( n>0 , sx/n , 0 )
		sd.q <- ifelse( n>0 , sd.old/sqrt(n) , 1)
	
		m.prop <- rnorm(Kw,mean.q,sd.q)
		
		q.m.num <- sum( log( dnorm( m.prop , mean.q , sd.q ) ) )
		q.m.den <- sum ( log( dnorm( m.old , mean.q , sd.q ) ) )
	
		post.num <- lpost.mix.augMC(ysam,z.new,ww.new,
								m.prop,sd.old,xi.old,K.old)	
		num <- post.num + q.m.num
		den <- post.den + q.m.den

		u <- runif(1)
		if(u<exp(num-den)){
			m.new <- m.prop
			post.den <- post.num
			acc.m <- acc.m + 1
		} else {
			m.new <- m.old
			post.den <- post.den
		}
	
		### Update of the standard deviations
		
		sv <- c()
		sd2.prop <-c()	
		q.sd.num <- 0
		q.sd.den <- 0
	
		for(j in 1:Kw)
		{
			sv[j]<-sum( (ysam[z.mat[,j]==1] - m.new[j])^2 )
			alpha.sigma <- 0.01 + (n[j]+1)/2
			beta.sigma <- 0.01 + 0.5 * 0.01 * m.new[j]^2 + 0.5*sv[j]
			sd2.prop[j] <- rinvgamma( 1 , alpha.sigma, beta.sigma)	
			q.sd.num <- q.sd.num + dinvgamma( sd.old[j]^2 , alpha.sigma , beta.sigma )
			q.sd.den <- q.sd.den + dinvgamma( sd2.prop[j] , alpha.sigma , beta.sigma )

		}
	
		sd.prop <- sqrt(sd2.prop)	
		
		post.num <- lpost.mix.augMC(ysam,z.new,ww.new,
								m.new,sd.prop,xi.old,K.old)	
	
		num <- post.num + q.sd.num
		den <- post.den + q.sd.den

		u <- runif(1)
		if(u<exp(num-den)){
			sd.new <- sd.prop
			post.den <- post.num
			acc.sd <- acc.sd + 1
		} else {
			sd.new <- sd.old
			post.den <- post.den
		}

		### Update of xi
			
		xi.prop <- rnorm(1,xi.old,0.1)
				
		post.num <- lpost.mix.augMC(ysam,z.new,ww.new,
									m.new,sd.new,xi.prop,K.old)		
		num <- post.num
		den <- post.den
		
		u <- runif(1)
		if(u<exp(num-den)){
			xi.new <- xi.prop
			post.den <- post.num
			acc.xi <- acc.xi + 1
		} else {
			xi.new <- xi.old
			post.den <- post.den
		}

		### Update of K
		
		K.prop <- rtruncnorm(1,a=0,b=Inf,mean=K.old,sd=0.01)

		q.K.num <- dtruncnorm(K.old,a=0,b=Inf,mean=K.prop,sd=0.01)
		q.K.den <- dtruncnorm(K.prop,a=0,b=Inf,mean=K.old,sd=0.01)

		post.num <- lpost.mix.augMC(ysam,z.new,ww.new,
									m.new,sd.new,xi.new,K.prop)			
		num <- post.num + q.K.num
		den <- post.den + q.K.den
		
		u <- runif(1)
		if(u<exp(num-den)){
			K.new <- K.prop
			post.den <- post.num
			acc.K <- acc.K + 1
		} else {
			K.new <- K.old
			post.den <- post.den
		}
					
		param.post[count,1:Kw] <- ww.new
		param.post[count,(Kw+1):(2*Kw)] <- m.new
		param.post[count,(2*Kw+1):(3*Kw)] <- sd.new
		param.post[count,(3*Kw+1)] <- xi.new
		param.post[count,(3*Kw+2)] <- K.new	
		
		z.old <- z.new

		lik.val <- c(lik.val, llik.mix.aug(ysam,z.new,ww.new,m.new,sd.new))
		post.val <- c(post.val,post.den)	
		
		z.freq[count,] <- apply(z.old,2,mean)

		if(round(count/10^4) ==(count/10^4)){
			save.image(“jeffmix.RData")		
		
		}
  
	}

	return(list(param.post=param.post , lik.seq=lik.val , post.seq = post.val , 
			accept=c(	acc.ww , acc.m , acc.sd , acc.xi , acc.K ) ) )
}		