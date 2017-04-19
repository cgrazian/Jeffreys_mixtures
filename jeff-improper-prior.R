# Target distribution

targKcomp=function(x,ww,mmean,ssd){sum(ww*dnorm(x,mmean,ssd))}

intMC=function(ysam,elem)
{
	ysam=as.matrix(ysam)
	temp=apply(ysam,1,elem)
	return(mean(temp)) 
}

# Likelihood

llik.mm=function(x,ww,mmean,ssd)
{
	vec=c()
	for(i in 1:length(x))
	{
		vec[i]=log(targKcomp(x[i],ww,mmean,ssd))
	}
	return(sum(vec))
}

### Numerical integration

FishInfoKcomp.wms=function(ww,mmean,ssd){

	Kw=length(ww)
	limsup.targ=mmean[Kw]+5*ssd[Kw]
	liminf.targ=mmean[1]-5*ssd[1]

	mm=matrix(0,nrow=(Kw*Kw-1),ncol=(Kw*Kw-1))
 
	### Sub-matrix for the weights

	mm.ww=matrix(0,nrow=(Kw-1),ncol=(Kw-1))
	for(i in 1:(Kw-1)){
		for(j in i:(Kw-1)){

			elem=function(x){
				(dnorm(x,mmean[i],ssd[i])-dnorm(x,mmean[Kw],ssd[Kw]))*
					(dnorm(x,mmean[j],ssd[j])-dnorm(x,mmean[Kw],ssd[Kw]))/
						targKcomp(x,ww,mmean,ssd)
			}
			mm.ww[i,j]=integrate(Vectorize(elem,vectorize.args='x'),upper=limsup.targ,
				lower=liminf.targ)$value
 
		}
	}
	mm[1:(Kw-1),1:(Kw-1)]=mm.ww

	### Sub-matrix for the means

	mm.mm=matrix(0,nrow=Kw,ncol=Kw)
	for(i in 1:Kw){
		elem.diag=function(x){
			(ww[i]^2/ssd[i]^4)*((x-mmean[i])*dnorm(x,mmean[i],ssd[i]))^2/
				(targKcomp(x,ww,mmean,ssd))
		}
		mm.mm[i,i]=integrate(Vectorize(elem.diag,vectorize.args='x'),upper=limsup.targ,
			lower=liminf.targ)$value
	}

	for(i in 1:(Kw-1)){
		for(j in (i+1):Kw){
			elem.offd=function(x){
				ww[i]*((x-mmean[i])/ssd[i]^2)*dnorm(x,mmean[i],ssd[i])*
					ww[j]*((x-mmean[j])/ssd[j]^2)*dnorm(x,mmean[j],ssd[j])/				
						(targKcomp(x,ww,mmean,ssd))
			}
			mm.mm[i,j]=integrate(Vectorize(elem.offd,vectorize.args='x'),upper=limsup.targ,
				lower=liminf.targ)$value
		}
	}
	mm[Kw:(2*Kw-1),Kw:(2*Kw-1)]=mm.mm

	### Rectangular sub-matrix for weights and means

	mm.wm=matrix(0,nrow=(Kw-1),ncol=Kw)

	for(i in 1:(Kw-1)){
		for(j in 1:(Kw-1)){

			elem.diag=function(x){
				( (ww[i]*(x-mmean[i])/ssd[i]^2)*dnorm(x,mmean[i],ssd[i]) )*
				( dnorm(x,mmean[i],ssd[i])-dnorm(x,mmean[Kw],ssd[Kw]) ) / 
					targKcomp(x,ww,mmean,ssd) 
			}

			elem.offd=function(x){
				( (ww[j]*(x-mmean[j])/ssd[j]^2)*dnorm(x,mmean[j],ssd[j]) )*
				( dnorm(x,mmean[i],ssd[i])-dnorm(x,mmean[Kw],ssd[Kw]) ) / 
					targKcomp(x,ww,mmean,ssd) 
			}

			if(i==j){
				mm.wm[i,j]=integrate(Vectorize(elem.diag,vectorize.args='x'),upper=limsup.targ,
						lower=liminf.targ)$value		
			} else {
				mm.wm[i,j]=integrate(Vectorize(elem.offd,vectorize.args='x'),upper=limsup.targ,
					lower=liminf.targ)$value}		
		
		
		}

		elem.last=function(x){
			( (ww[Kw]*(x-mmean[Kw])/ssd[Kw]^2)*dnorm(x,mmean[Kw],ssd[Kw]) )*
			( dnorm(x,mmean[i],ssd[i])-dnorm(x,mmean[Kw],ssd[Kw]) ) / 
				targKcomp(x,ww,mmean,ssd) 

		}
		mm.wm[i,Kw]=integrate(Vectorize(elem.last,vectorize.args='x'),upper=limsup.targ,
				lower=liminf.targ)$value		
	}
	mm[1:(Kw-1),Kw:(2*Kw-1)]=mm.wm

	### Sub-matrix for the standard deviations

	mm.ss=matrix(0,nrow=Kw,ncol=Kw)

 	for(i in 1:Kw){
		elem.diag=function(x){
			(ww[i]*dnorm(x,mmean[i],ssd[i])*(((x-mmean[i])/ssd[i])^2-1)/
				ssd[i])^2/
					(targKcomp(x,ww,mmean,ssd))
		}    
		mm.ss[i,i]=integrate(Vectorize(elem.diag,vectorize.args='x'),upper=limsup.targ,
				lower=liminf.targ)$value
	}

	for(i in 1:(Kw-1)){
		for(j in (i+1):Kw){
			elem.offd=function(x){
				(ww[i]*dnorm(x,mmean[i],ssd[i])*(((x-mmean[i])/ssd[i])^2-1)/ssd[i])*
				(ww[j]*dnorm(x,mmean[j],ssd[j])*(((x-mmean[j])/ssd[j])^2-1)/ssd[j])/		
					(targKcomp(x,ww,mmean,ssd))
			}
			mm.ss[i,j]=integrate(Vectorize(elem.offd,vectorize.args='x'),upper=limsup.targ,
				lower=liminf.targ)$value
			}
	}
	mm[(2*Kw):(3*Kw-1),(2*Kw):(3*Kw-1)]=mm.ss

	### Rectangular sub-matrix for the standard deviations and the weights

	mm.ws=matrix(0,nrow=(Kw-1),ncol=(Kw))

	for(i in 1:(Kw-1)){
		for(j in 1:(Kw-1)){

			elem.diag=function(x){
			( ww[i] * ( (x-mmean[i])^2 / ssd[i]^2 - 1 ) / ssd[i] * 
				(dnorm(x,mmean[i],ssd[i]) - dnorm(x,mmean[Kw],ssd[Kw])) * 
					dnorm(x,mmean[i],ssd[i]) ) / 
						targKcomp(x,ww,mmean,ssd)
			}
			elem.offd=function(x){
			( ww[j] * ( (x-mmean[j])^2 / ssd[j]^2 - 1 ) / ssd[j] * 
				(dnorm(x,mmean[i],ssd[i]) - dnorm(x,mmean[Kw],ssd[Kw])) * 
					dnorm(x,mmean[j],ssd[j]) ) / 
						targKcomp(x,ww,mmean,ssd)
			}

			if(i==j){
				mm.ws[i,j]=integrate(Vectorize(elem.diag,vectorize.args='x'),upper=limsup.targ,
						lower=liminf.targ)$value		
			} else {
				mm.ws[i,j]=integrate(Vectorize(elem.offd,vectorize.args='x'),upper=limsup.targ,
					lower=liminf.targ)$value}		
		}

		elem.last=function(x){
			( ww[Kw] * ( (x-mmean[Kw])^2 / ssd[Kw]^2 - 1 ) / ssd[Kw] * 
				(dnorm(x,mmean[i],ssd[i]) - dnorm(x,mmean[Kw],ssd[Kw])) * 
					dnorm(x,mmean[Kw],ssd[Kw]) ) / 
						targKcomp(x,ww,mmean,ssd)
		}
		mm.ws[i,Kw]=integrate(Vectorize(elem.last,vectorize.args='x'),upper=limsup.targ,
			lower=liminf.targ,stop.on.error=F)$value
	}
	mm[1:(Kw-1),(2*Kw):(3*Kw-1)]=mm.ws

	### Rectangular sub-matrix for the standard deviations and the means

	mm.ms=matrix(0,nrow=Kw,ncol=Kw)

	for(i in 1:Kw){
		for(j in 1:Kw){

			elem.diag=function(x){
				( ww[i]*dnorm(x,mmean[i],ssd[i]) )^2 *
				( (x-mmean[i])/ssd[i]^3 ) * ( (x-mmean[i])^2/ssd[i]^2 - 1 ) / 
					targKcomp(x,ww,mmean,ssd)
			}

			elem.offd=function(x){
				ww[i]*dnorm(x,mmean[i],ssd[i]) *
					( (x-mmean[i])/ssd[i]^2 ) * 
				ww[j]*dnorm(x,mmean[j],ssd[j]) *
					( (x-mmean[j])^2/ssd[j]^2 - 1 ) /
					targKcomp(x,ww,mmean,ssd)
			}


			if(i==j){
				mm.ms[i,j]=integrate(Vectorize(elem.diag,vectorize.args='x'),upper=limsup.targ,
						lower=liminf.targ)$value		
			} else {
				mm.ms[i,j]=integrate(Vectorize(elem.offd,vectorize.args='x'),upper=limsup.targ,
					lower=liminf.targ)$value}		
		}
	}
	mm[Kw:(2*Kw-1),(2*Kw):(3*Kw-1)]=mm.ms

	mm1=mm+t(mm)
	diag(mm1)=diag(mm)

	return(mm1)

}

### Monte Carlo integration

FishInfoKcomp.wmsMC=function(ww,mmean,ssd,nMCMC){

	Kw=length(ww)

	nsamMC=nMCMC
	ncompMC=ww*nsamMC
	MCysam=rnorm(ncompMC[1],mmean[1],ssd[1])
	for(u in 2:length(ncompMC)){	
		MCysam=c(MCysam,rnorm(ncompMC[u],mmean[u],ssd[u]))}

	mm=matrix(0,nrow=(Kw*Kw-1),ncol=(Kw*Kw-1))
 
	### Sub-matrix for the weights

	mm.ww=matrix(0,nrow=(Kw-1),ncol=(Kw-1))
	for(i in 1:(Kw-1)){
		for(j in i:(Kw-1)){

			elem=function(x){
				(dnorm(x,mmean[i],ssd[i])-dnorm(x,mmean[Kw],ssd[Kw]))*
					(dnorm(x,mmean[j],ssd[j])-dnorm(x,mmean[Kw],ssd[Kw]))/
						(targKcomp(x,ww,mmean,ssd)^2)
			}
			mm.ww[i,j]=intMC(MCysam,elem)

 
		}
	}
	mm[1:(Kw-1),1:(Kw-1)]=mm.ww

	### Sub-matrix for the means

	mm.mm=matrix(0,nrow=Kw,ncol=Kw)
	for(i in 1:Kw){
		elem.diag=function(x){
			(ww[i]^2/ssd[i]^4)*((x-mmean[i])*dnorm(x,mmean[i],ssd[i]))^2/
				(targKcomp(x,ww,mmean,ssd)^2)
		}
		mm.mm[i,i]=intMC(MCysam,elem.diag)

	}

	for(i in 1:(Kw-1)){
		for(j in (i+1):Kw){
			elem.offd=function(x){
				ww[i]*((x-mmean[i])/ssd[i]^2)*dnorm(x,mmean[i],ssd[i])*
					ww[j]*((x-mmean[j])/ssd[j]^2)*dnorm(x,mmean[j],ssd[j])/				
						(targKcomp(x,ww,mmean,ssd)^2)
			}
			mm.mm[i,j]=intMC(MCysam,elem.offd)

		}
	}
	mm[Kw:(2*Kw-1),Kw:(2*Kw-1)]=mm.mm

	### Rectangular sub-matrix for weights and means

	mm.wm=matrix(0,nrow=(Kw-1),ncol=Kw)

	for(i in 1:(Kw-1)){
		for(j in 1:(Kw-1)){

			elem.diag=function(x){
				( ww[i] * (x-mmean[i]) * dnorm(x,mmean[i],ssd[i]) / 
					ssd[i]^2 ) *
				( dnorm(x,mmean[i],ssd[i])-dnorm(x,mmean[Kw],ssd[Kw]) ) / 
					(targKcomp(x,ww,mmean,ssd)^2)
			}

			elem.offd=function(x){
				( ww[j] * (x-mmean[j]) * dnorm(x,mmean[j],ssd[j]) / 
					ssd[j]^2 ) *
				( dnorm(x,mmean[i],ssd[i])-dnorm(x,mmean[Kw],ssd[Kw]) ) / 
					(targKcomp(x,ww,mmean,ssd)^2)
			}

			if(i==j){
				mm.wm[i,j]=intMC(MCysam,elem.diag)
		
			} else {
				mm.wm[i,j]=intMC(MCysam,elem.offd)
			}
		
		}

		elem.last=function(x){
			( (ww[Kw]*(x-mmean[Kw])/ssd[Kw]^2)*dnorm(x,mmean[Kw],ssd[Kw]) )*
			( dnorm(x,mmean[i],ssd[i])-dnorm(x,mmean[Kw],ssd[Kw]) ) / 
				(targKcomp(x,ww,mmean,ssd)^2) 
		}
		mm.wm[i,Kw]=intMC(MCysam,elem.last)
		
	}
	mm[1:(Kw-1),Kw:(2*Kw-1)]=mm.wm

	### Sub-matrix for the standard deviations

	mm.ss=matrix(0,nrow=Kw,ncol=Kw)

 	for(i in 1:Kw){
		elem.diag=function(x){
			(ww[i]*dnorm(x,mmean[i],ssd[i])*(((x-mmean[i])/ssd[i])^2-1)/
				ssd[i])^2/(targKcomp(x,ww,mmean,ssd)^2)
		}    
		mm.ss[i,i]=intMC(MCysam,elem.diag)

	}

	for(i in 1:(Kw-1)){
		for(j in (i+1):Kw){
			elem.offd=function(x){
				(ww[i]*dnorm(x,mmean[i],ssd[i])*(((x-mmean[i])/ssd[i])^2-1)/ssd[i])*
				(ww[j]*dnorm(x,mmean[j],ssd[j])*(((x-mmean[j])/ssd[j])^2-1)/ssd[j])/		
					(targKcomp(x,ww,mmean,ssd)^2)
			}
			mm.ss[i,j]=intMC(MCysam,elem.offd)

			}
	}
	mm[(2*Kw):(3*Kw-1),(2*Kw):(3*Kw-1)]=mm.ss

	### Rectangular sub-matrix for the standard deviations and the weights

	mm.ws=matrix(0,nrow=(Kw-1),ncol=(Kw))

	for(i in 1:(Kw-1)){
		for(j in 1:(Kw-1)){

			elem.diag=function(x){
			( ww[i] * ( (x-mmean[i])^2 / ssd[i]^2 - 1 ) / ssd[i] * 
				(dnorm(x,mmean[i],ssd[i]) - dnorm(x,mmean[Kw],ssd[Kw])) * 
					dnorm(x,mmean[i],ssd[i]) ) / 
						(targKcomp(x,ww,mmean,ssd)^2)
			}
			elem.offd=function(x){
			( ww[j] * ( (x-mmean[j])^2 / ssd[j]^2 - 1 ) / ssd[j] * 
				(dnorm(x,mmean[i],ssd[i]) - dnorm(x,mmean[Kw],ssd[Kw])) * 
					dnorm(x,mmean[j],ssd[j]) ) / 
						(targKcomp(x,ww,mmean,ssd)^2)
			}

			if(i==j){
				mm.ws[i,j]=intMC(MCysam,elem.diag)
		
			} else {
				mm.ws[i,j]=intMC(MCysam,elem.offd)

			}		
		}
		elem.last=function(x){
			( ww[Kw] * ( (x-mmean[Kw])^2 / ssd[Kw]^2 - 1 ) / ssd[Kw] * 
				(dnorm(x,mmean[i],ssd[i]) - dnorm(x,mmean[Kw],ssd[Kw])) * 
					dnorm(x,mmean[Kw],ssd[Kw]) ) / 
						(targKcomp(x,ww,mmean,ssd)^2)
		}
		mm.ws[i,Kw]=intMC(MCysam,elem.last)

	}
	mm[1:(Kw-1),(2*Kw):(3*Kw-1)]=mm.ws

	### Rectangular sub-matrix for the standard deviations and the means

	mm.ms=matrix(0,nrow=Kw,ncol=Kw)

	for(i in 1:Kw){
		for(j in 1:Kw){

			elem.diag=function(x){
				( ww[i]*dnorm(x,mmean[i],ssd[i]) )^2 *
				( (x-mmean[i])/ssd[i]^3 ) * ( (x-mmean[i])^2/ssd[i]^2 - 1 ) / 
					(targKcomp(x,ww,mmean,ssd)^2)
			}

			elem.offd=function(x){
				ww[i]*dnorm(x,mmean[i],ssd[i]) *
					( (x-mmean[i])/ssd[i]^2 ) * 
				ww[j]*dnorm(x,mmean[j],ssd[j]) *
					( (x-mmean[j])^2/ssd[j]^2 - 1 ) /
					(targKcomp(x,ww,mmean,ssd)^2)
			}

			if(i==j){
				mm.ms[i,j]=intMC(MCysam,elem.diag)
		
			} else {
				mm.ms[i,j]=intMC(MCysam,elem.offd)
			}		
		}
	}
	mm[Kw:(2*Kw-1),(2*Kw):(3*Kw-1)]=mm.ms

	mm1=mm+t(mm)
	diag(mm1)=diag(mm)

	return(mm1)

}

####################
### Jeffreys prior
####################

jeffKcomp.wms=function(ww,mmean,ssd,nMCMC){
	if(sum(ssd<0.2)==0)
	{
		det.num=det(FishInfoKcomp.wms(ww,mmean,ssd))
		if(det.num<0)
		{
			cont=0
			det.num=-1
			while(det.num<=0&cont<5)
			{
				det.num=det(FishInfoKcomp.wmsMC(ww,mmean,ssd,nMCMC))
				cont=cont+1
			}
			det.num=ifelse(det.num<0,0,det.num)
		}
	} else {
		cont=0
		det.num=-1
		while(det.num<=0&cont<5)
		{
			det.num=det(FishInfoKcomp.wmsMC(ww,mmean,ssd,nMCMC))
			cont=cont+1
		}
		det.num=ifelse(det.num<0,0,det.num)
	}
	return(sqrt(det.num))
}

