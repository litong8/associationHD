###try a list of sigma_err sequence while keep all other parameters the same as our main simulation setting
###for each value in the sequence, generate diseased group data with large N (for example 50000), 
## for simplicity, we can have one time point as t=0.
###Compute SNR for all Q I using the formula to calculate SNR=diag(Var(Istar))/diag((Var(I)-Var(Istar))) 
## only among those t=0.
###Compute average SNR by mSNR=mean(SNR)
###Find the three sigma_err values such that mSNR is close to value we want (i.e., 1,2,4,10)
###Then just copy those selected 3 sigma_err into the simulation code and loop through it like 
## what you do for different N.

library(MASS)
library(Matrix)

####Define variance process model here
Sigma_AR<-function(tseq,theta){
	Sigma=theta[1]*theta[2]^abs(outer(tseq,tseq,"-"))
	return(as.matrix(Sigma))
}
Sigma_CS<-function(tseq,theta){
	Sigma=theta[1]*theta[2]^(abs(outer(tseq,tseq,"-"))>0)
	return(as.matrix(Sigma))
}

####Data generation process
datagen<-function(N,K,Q,p,tseq,beta0,beta1,beta2,b0,b1,b2,b3,Sigma_a,Sigma_X,Sigma,Sigma_e,Sigma_err,rho, rho_e,rho_err,gamma,seed,type="AR"){
	###N: Sample size
	###K: Number of clinical domains
	###Q: Number of imaging domains
	###p: Number of covariates
	###beta0: imaging intercept
	###beta1: imaging slope
	###beta2: imaging coefficient of X
	###tseq: a vector of measurement times for balanced design (might extend to a matrix of different measurement times for each individual)
	###Sigma_a: a 2Q*2Q matrix for the variance-covariance matrix of random effects
	###Sigma_X: a p*p matrix for the variance-covariance matrix of covariates
	###Sigma:variance-covariance structure at a specific time for imaging domains
	###Sigma_e:variance-covariance structure at a specific time for clinical domains
	###Sigma_err:variance-covariance structure at a specific time for measurement errors
	###rho:parameter for variance-covariance structure over time for imaging domains
	###rho_e:parameter for variance-covariance structure over time for clinical domains
	###rho_err:parameter for variance-covariance structure over time for measurement errors
	###seed: simulation seed
	set.seed(seed)
	###generate random effect a_i
	a_i=mvrnorm(N,mu=rep(0,2*Q),Sigma=Sigma_a)
	###generate covariate
	X=mvrnorm(N,mu=rep(0,p),Sigma=Sigma_X)
	###variance-covariance structure over time for imaging domains
	if (type=="AR"){
		Sigma_t=Sigma_AR(tseq,c(1,rho))
		###variance-covariance structure over time for clinical domains
		Sigma_e_t=Sigma_AR(tseq,c(1,rho_e))
		###variance-covariance structure over time for measurement errors
		Sigma_err_t=Sigma_AR(tseq,c(1,rho_err))
	}
	if (type=="CS"){
		Sigma_t=Sigma_CS(tseq,c(1,rho))
		###variance-covariance structure over time for clinical domains
		Sigma_e_t=Sigma_CS(tseq,c(1,rho_e))
		###variance-covariance structure over time for measurement errors
		Sigma_err_t=Sigma_CS(tseq,c(1,rho_err))
	}
	data=NULL
	###generate all residuals
	###assume kronecker structure for variance-covariance over time at different domains. Might need more complicated form for the transformed model
	SSigma_t=kronecker(Sigma_t, Sigma)
	epi_t_all=mvrnorm(N,mu=rep(0,nrow(SSigma_t)),Sigma=SSigma_t)
	SSigma_e_t=kronecker(Sigma_e_t, Sigma_e)
	e_t_all=mvrnorm(N,mu=rep(0,nrow(SSigma_e_t)),Sigma=SSigma_e_t)
	SSigma_err_t=kronecker(Sigma_err_t, Sigma_err)
	err_t_all=mvrnorm(N,mu=rep(0,nrow(SSigma_err_t)),Sigma=SSigma_err_t)

	###generate data at each time
	for (i in 1:length(tseq)){
		t=tseq[i]
		###Generate imaging domains
		epi_t=epi_t_all[,(i-1)*Q+c(1:Q)]
		for (q in 1:Q){
			mu_I=beta0[q]+beta1[q]*t+X%*%beta2[q,]+a_i[,2*(q-1)+c(1:2)]%*%c(1,t)
			epi_t[,q]=mu_I+epi_t[,q]
		}
		###Generate clinical domains (using I(t)), add learning effect gamma*t
		e_t=e_t_all[,(i-1)*K+c(1:K)]+rep(1,N)%o%(gamma*t)
		for (k in 1:K){
			mu_S=b0[k]+b1[k]*t+X%*%b2[k,]+epi_t%*%b3[k,]
			e_t[,k]=mu_S+e_t[,k]
		}
		####add measurement error
		err_t=err_t_all[,(i-1)*(Q+K)+c(1:(Q+K))]
		tmp=cbind(1:N,t,X,cbind(epi_t,e_t)+err_t,epi_t,e_t)
		data=rbind(data,tmp)
	}
	###rename variables
	###use star to denote the truth here
	xname=paste("X",as.character(1:p),sep="")
	iname=paste("I",as.character(1:Q),sep="")
	sname=paste("S",as.character(1:K),sep="")
	iname1=paste("Istar",as.character(1:Q),sep="")
	sname1=paste("Sstar",as.character(1:K),sep="")
	data=data.frame(data)
	names(data)=c("ID","t",xname,iname,sname,iname1,sname1)
      return(list(data=data))	
}	

################################################
## set parameters for data generation

seed=1111
####run simulation example
set.seed(seed)
####example setting
K=4
## length(gamma)=K;    gamma=c(0.2,0.3,0.1,0.4,0.5,0.9)
gamma=c(0.1,0.4,0.5,0.9)

Q=3
N=50000
tseq=0
seed=1111
sigma_a=0.2
###genearte a random positive definite matrix from Wishart distribution
#Sigma_a=sigma_a*rWishart(1,2*Q,diag(2*Q))[,,1]
Sigma_a=Sigma_CS(1:(2*Q),c(1,-0.1))
p=2
rho_X=0.1
Sigma_X=Sigma_CS(1:p,c(1,rho_X))

###Define correlation over time
rho=0.3
rho_e=0.2
rho_err=0

####generate Sigma randomly from Wishart distribution with a scalar sigma
#sigma=5
#Sigma=diag(diag(sigma*rWishart(1,Q,diag(Q))[,,1]))
Sigma=Sigma_AR(1:Q,c(3,0.2))

####generate Sigma_e randomly from Wishart distribution with a scalar sigma_e
#sigma_e=0.3
#Sigma_e=sigma_e*rWishart(1,K,diag(K))[,,1]
Sigma_e=Sigma_CS(1:K,c(2.5,0.15))
Sigma_e.true=Sigma_e
Sigma_e_t.true=Sigma_AR(tseq,c(1,rho_e))

####generate random mean parameters
#beta0=rnorm(Q)
#beta1=rnorm(Q)
#beta2=matrix(data=rnorm(Q*p),nrow=Q,ncol=p)
#b0=rnorm(K)
#b1=rnorm(K)
#b2=matrix(data=rnorm(K*p),nrow=K,ncol=p)
#b3=matrix(data=rnorm(K*Q),nrow=K,ncol=Q)

beta0=rep(1,Q)
beta1=1:Q
beta2=(1:K)%o%rep(0.5,p)
b0=rep(1,K)
b1=1:K
b2=(1:K)%o%rep(0.5,p)
b3=(K:1)%o%rep(0.8,Q)

beta_true=cbind(b0,b1,b2,b3)

iname=paste("I",as.character(1:Q),sep="")
iname1=paste("Istar",as.character(1:Q),sep="")

## Signal-noise ratio
mSNR=NULL
####generate Sigma_err (include Sigma_I,Sigma_S,Sigma_IS) randomly from Wishart distribution with a scalar sigma_err
sigma_err=seq(0,1.4,0.01)
for (i in sigma_err){
Sigma_err=i*rWishart(1,K+Q,diag(K+Q))[,,1]
#Sigma_err=diag(c(rep(1,K),rep(1.2,Q)))

alldata=datagen(N,K,Q,p,tseq,beta0,beta1,beta2,b0,b1,b2,b3,Sigma_a,Sigma_X,Sigma,Sigma_e,Sigma_err,rho, rho_e,rho_err,gamma,seed,type="AR")
data=alldata$data
I=data[,iname]
Istar=data[,iname1]
SNR=diag(var(Istar))/diag(var(I)-var(Istar))
mSNR=c(mSNR,mean(SNR))
}
write.csv(data.frame(t(rbind(sigma_err,mSNR))),paste("sigma_err vs. SNR",".csv")) 



