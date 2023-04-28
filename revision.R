## revision 1: set n0=n/3
library(MASS)
library(Matrix)

## Simseed These lines cannot be directly run in R, it need to be run using the submit file. Because the code
## R CMD BATCH "--args $SLURM_ARRAY_TASK_ID" Simulation_cor2_est_parallel.R is running with the args value assigned from the name 
## of SLURM_ARRAY_TASK_ID which vary among different tasks.

args=commandArgs(trailingOnly=TRUE)
simseed=8888+5*as.numeric(args[1])

####Simulation for step 1 models
####Data generation process

####Define variance process model here
Sigma_AR<-function(tseq,theta){
	Sigma=theta[1]*theta[2]^abs(outer(tseq,tseq,"-"))
	return(as.matrix(Sigma))
}
Sigma_CS<-function(tseq,theta){
	Sigma=theta[1]*theta[2]^(abs(outer(tseq,tseq,"-"))>0)
	return(as.matrix(Sigma))
}


###Data generation code
datagen<-function(N,K,Q,p,tseq,beta0,beta1,beta2,b0,b1,b2,b3,Sigma_a,Sigma_X,Sigma,Sigma_e,Sigma_err,rho, rho_e,rho_err,gamma,seed,type){
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
	if (type=="AR"){
            ###variance-covariance structure over time for imaging domains
		Sigma_t=Sigma_AR(tseq,c(1,rho))
		###variance-covariance structure over time for clinical domains
		Sigma_e_t=Sigma_AR(tseq,c(1,rho_e))
		###variance-covariance structure over time for measurement errors
		Sigma_err_t=Sigma_AR(tseq,c(1,rho_err))
	}
	if (type=="CS"){
            ###variance-covariance structure over time for imaging domains
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
    ## the learning effect is a linear function of how many times you get the test 
    ## (i.e., t). If you use constant, then you are assuming the skills you gained in 
    ## 2 tests are the same as the skills you gained in 3 tests.

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
	
####now generate normal population using half of the sample
	N=floor(N/3)
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
	data0=NULL
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
			mu_I=beta0[q]+beta1[q]*0+X%*%beta2[q,]+a_i[,2*(q-1)+c(1:2)]%*%c(1,0)
			epi_t[,q]=mu_I
		}
		###Generate clinical domains (using I(t))
		e_t=e_t_all[,(i-1)*K+c(1:K)]
		for (k in 1:K){
            ## assume I(t) is constant over t for each individual
			mu_S=b0[k]+b1[k]*0+X%*%b2[k,]
			e_t[,k]=mu_S
		}
		e_t=e_t+rep(1,N)%o%(gamma*t)
		####add measurement error
		err_t=err_t_all[,(i-1)*(Q+K)+c(1:(Q+K))]
		tmp=cbind(1:N,t,X,cbind(epi_t,e_t)+err_t,epi_t,e_t)
		data0=rbind(data0,tmp)
	}
	###rename variables
	###use star to denote the truth here
	data0=data.frame(data0)
	names(data0)=c("ID","t",xname,iname,sname,iname1,sname1)	
	return(list(data=data,data0=data0))	
}  # end datagen

################# some useful functions

### l2 distance between 2 vectors
mydist<-function(u,v){
	sqrt(sum((u-v)^2))
}

### Box kernel
myK<-function(x){
	as.numeric(abs(x)<=1)/2
}

###trace
tr<-function(mat){
	return(sum(diag(mat)))
}

###mat to vec   lower triangle of matrix    
mat2vec<-function(mat){
	d=nrow(mat)
	v=NULL
	for (i in 1:d){
		for (j in 1:i){
			v=c(v,mat[i,j])
		}
	}
	v
}

###vec to mat
vec2mat<-function(v){
	d=floor(sqrt(2*length(v)))
	count=0
	mat=matrix(data=NA,nrow=d,ncol=d)
	for (i in 1:d){
		for (j in 1:i){
			count=count+1
			mat[i,j]=mat[j,i]=v[count]
		}
	}
	mat
}

########################################

###Model fitting code 
###OLS fitting code  (observed data) Xmat=t,X,I   Ymat=S,
## (independent working correlation over time), no measurement error correction 
myfit_ols<-function(data,K,Q,p){
	###OLS fitting
	Xmat=as.matrix(cbind(1,data[,c(2:(2+p+Q))]))
	Ymat=as.matrix(data[,2+p+Q+c(1:K)])
	YY=t(Xmat)%*%Ymat
	XX=t(Xmat)%*%Xmat
	beta_ols=t(solve(XX,YY))
	beta_ols
}
###Oracle fitting code, assuming error-free variables are observed
## independence working correlation
myfit_oracle<-function(data,K,Q,p){
	###oracle fitting  (true data Xmat=(t,X,I*)  Ymat=S*	
	Xmatstar=as.matrix(cbind(1,data[,c(2:(2+p),2+p+Q+K+c(1:Q))]))
	Ymatstar=as.matrix(data[,2+p+Q+K+Q+c(1:K)])
	YYstar=t(Xmatstar)%*%Ymatstar
	XXstar=t(Xmatstar)%*%Xmatstar
	beta_star=t(solve(XXstar,YYstar))
	beta_star
}

###Fitting code assuming independence working correlation, 
## includes measurement error adjustment
myfit_ind<-function(data,K,Q,p,Sigma_err){
	N=nrow(data)/length(tseq)
	###initial estimation 
	Sigma_I=Sigma_err[1:Q,1:Q]
	Sigma_S=Sigma_err[Q+c(1:K),Q+c(1:K)]
	Sigma_IS=Sigma_err[1:Q,Q+c(1:K)]
	###Measurement error correction
	Xmat=as.matrix(cbind(1,data[,c(2:(2+p+Q))]))
	Ymat=as.matrix(data[,2+p+Q+c(1:K)])
	YY=t(Xmat)%*%Ymat
	XX=t(Xmat)%*%Xmat
	####make sure positivity of XX, adjustment for measurement error, S has N*n_i elements
	## although XX is positive definite, after you subtract another positive definite 
      ## matrix lambda* nrow(Xmat)*Sigma_I from it, it can become negative, so the minimum 
      ## eigen can be negative. Here using the absolute value is for the purpose of solving 
      ## the lambda such that the minimum eigen value equals 0.
      
      YY[2+p+c(1:Q),1:K]=YY[2+p+c(1:Q),1:K]-nrow(Xmat)*Sigma_IS
	if (min(eigen(XX[2+p+c(1:Q),2+p+c(1:Q)]-nrow(Xmat)/(nrow(Xmat)-1)*nrow(Xmat)*Sigma_I)$values)>0){
		XX[2+p+c(1:Q),2+p+c(1:Q)]=XX[2+p+c(1:Q),2+p+c(1:Q)]-nrow(Xmat)*Sigma_I
	}else{
		myf<-function(lambda){
			abs(min(eigen(XX[2+p+c(1:Q),2+p+c(1:Q)]-lambda*nrow(Xmat)*Sigma_I)$values))
		}
		lambda=optimize(myf,c(0,nrow(Xmat)/(nrow(Xmat)-1)))$minimum
		XX[2+p+c(1:Q),2+p+c(1:Q)]=XX[2+p+c(1:Q),2+p+c(1:Q)]-(lambda-6/(nrow(Xmat)-1))*nrow(Xmat)*Sigma_I
	}	
	beta_ind=t(solve(XX,YY))
	beta_ind
}

###Oracle fitting code, assuming error-free variables are observed;
###Fitting code assuming specific working correlation matrix 
## either true or one-step working correlation Sigma_e, Sigma_e_t
myfit_oracle_cor<-function(data,K,Q,p,tseq,Sigma_e,Sigma_e_t){
	N=nrow(data)/length(tseq)
	##oracle fitting  (true data Xmatstar=(t,X,I*)  Ymatstar=S*	
	Xmatstar=as.matrix(cbind(1,data[,c(2:(2+p),2+p+Q+K+c(1:Q))]))
	Ymatstar=as.matrix(data[,2+p+Q+K+Q+c(1:K)])
	###estimate parameter
	WW=solve(Sigma_e_t)
	Hhat=matrix(data=0,nrow=2+p+Q,ncol=2+p+Q)
	for (j in 1:length(tseq)){
		for (k in 1:length(tseq)){
			Hhat=Hhat+c(WW[j,k])*t(Xmatstar[1:N+(j-1)*N,])%*%Xmatstar[1:N+(k-1)*N,]
		}
	}
	HHhat=kronecker(solve(Hhat),Sigma_e)
	###update U
	SSigma_e=kronecker(Sigma_e_t, Sigma_e)
	WWtilde=solve(SSigma_e)
	Uhat=matrix(data=0,nrow=(2+p+Q)*K,ncol=1)
	for (i in 1:N){
		for (j in 1:(length(tseq)*K)){
			for (l in 1:(length(tseq)*K)){
				jj=(j-1)%/%K+1
				ll=(l-1)%/%K+1
				jjj=(j-1)%%K+1
				lll=(l-1)%%K+1
				Uhat=Uhat+c(WWtilde[j,l])*kronecker(Xmatstar[i+((1:length(tseq))-1)*N,],diag(rep(1,K)))[j,]*Ymatstar[i+(ll-1)*N,lll]
			}
		}
	}
	beta_star_cor=matrix(data=HHhat%*%Uhat,nrow=K,ncol=2+p+Q,byrow=FALSE)
	beta_star_cor
}

###Fitting code assuming specific working correlation matrix, 
## includes measurement error adjustment
## either true or one-step working correlation Sigma_e, Sigma_e_t
myfit_cor<-function(data,K,Q,p,tseq,Sigma_err,Sigma_e,Sigma_e_t){
	N=nrow(data)/length(tseq)
	###initial estimation 
	Sigma_I=Sigma_err[1:Q,1:Q]
	Sigma_S=Sigma_err[Q+c(1:K),Q+c(1:K)]
	Sigma_IS=Sigma_err[1:Q,Q+c(1:K)]
	###Measurement error correction
	Xmat=as.matrix(cbind(1,data[,c(2:(2+p+Q))]))
	Ymat=as.matrix(data[,2+p+Q+c(1:K)])
	YY=t(Xmat)%*%Ymat
	XX=t(Xmat)%*%Xmat
	###estimate parameter
	WW=solve(Sigma_e_t)
	Hhat1=matrix(data=0,nrow=2+p+Q,ncol=2+p+Q)
	###expanded sigma_I
	Omega_I=matrix(data=0,nrow=2+p+Q,ncol=2+p+Q)
	Omega_I[2+p+(1:Q),2+p+(1:Q)]=Sigma_I
	Hhat2=N*tr(WW)*Omega_I
	for (j in 1:length(tseq)){
		for (k in 1:length(tseq)){
			Hhat1=Hhat1+c(WW[j,k])*t(Xmat[1:N+(j-1)*N,])%*%Xmat[1:N+(k-1)*N,]
		}
	}
	###need make sure positivity
	if (min(eigen(Hhat1-nrow(Xmat)/(nrow(Xmat)-1)*Hhat2)$values)>0){
		Hhat=Hhat1-Hhat2
	}else{
		myf<-function(lambda){
			abs(min(eigen(Hhat1-lambda*Hhat2)$values))
		}
		lambda=optimize(myf,c(0,nrow(Xmat)/(nrow(Xmat)-1)))$minimum
		Hhat=Hhat1-(lambda-6/(nrow(Xmat)-1))*Hhat2
	}	
	HHhat=kronecker(solve(Hhat),Sigma_e)
	
	###update U
	SSigma_e=kronecker(Sigma_e_t, Sigma_e)
	WWtilde=solve(SSigma_e)
	Omega_IS=matrix(data=0,nrow=2+p+Q,ncol=K)
	Omega_IS[2+p+(1:Q),]=Sigma_IS
	Uhat=matrix(data=0,nrow=(2+p+Q)*K,ncol=1)
	for (i in 1:N){
		for (j in 1:(length(tseq)*K)){
			for (l in 1:(length(tseq)*K)){
				jj=(j-1)%/%K+1
				ll=(l-1)%/%K+1
				jjj=(j-1)%%K+1
				lll=(l-1)%%K+1
				Uhat=Uhat+c(WWtilde[j,l])*kronecker(Xmat[i+((1:length(tseq))-1)*N,],diag(rep(1,K)))[j,]*Ymat[i+(ll-1)*N,lll]
				if (jj==ll){
					e_lll=rep(0,K)
					e_lll[lll]=1
					e_jjj=rep(0,K)
					e_jjj[jjj]=1
					Uhat=Uhat-c(WWtilde[j,l])*kronecker(Omega_IS%*%e_lll,e_jjj)
				}
			}
		}
	}
	beta_onestep=matrix(data=HHhat%*%Uhat,nrow=K,ncol=2+p+Q,byrow=FALSE)
	beta_onestep
}

###code to estimate parameter related to correlation structure: Sigma_e,Sigma_e_t
## myest_par, that Kernel estimated variance-covariance matrix SSigma_e_hat1 
## is not always positive definite when measurement error is large 
## and sometimes SSigma_e_hat are ill-conditioned

myest_par<-function(data,K,Q,p,tseq,beta_ind,Sigma_err,h,type){
	N=nrow(data)/length(tseq)
	###initial estimation 
	Sigma_I=Sigma_err[1:Q,1:Q]
	Sigma_S=Sigma_err[Q+c(1:K),Q+c(1:K)]
	Sigma_IS=Sigma_err[1:Q,Q+c(1:K)]

	###use one step update, estimate rho_e, Sigma_e
	###get empirical estimation for SSigma_e
	SSigma_e_hat1=matrix(data=0,nrow=K*length(tseq),ncol=K*length(tseq))
	b4=beta_ind[,2+p+c(1:Q)]	
	Xmat=as.matrix(cbind(1,data[,c(2:(2+p+Q))]))
	Ymat=as.matrix(data[,2+p+Q+c(1:K)])
	Ehat=Ymat-Xmat%*%t(beta_ind)
	for (jj in 1:length(tseq)){
		for (kk in 1:length(tseq)){
			upper=matrix(data=0,nrow=K,ncol=K)
			lower=0
			if (jj==kk){
				for (i in 1:N){
					for (j in 1:length(tseq)){
						ehat_ij=Ehat[i+(j-1)*N,]
						w_ij=myK(abs(tseq[j]-tseq[jj])/h)
						upper=upper+w_ij*(ehat_ij%o%ehat_ij)
						lower=lower+w_ij
					}
				}
				SSigma_e_hat1[(jj-1)*K+(1:K),(kk-1)*K+(1:K)]=upper/lower				
			}
			if (jj!=kk){
				for (i in 1:N){
					for (j in 1:length(tseq)){
						for (k in 1:length(tseq)){
							ehat_ij=Ehat[i+(j-1)*N,]
							ehat_ik=Ehat[i+(k-1)*N,]
							if (k==j){
								w_ijk=0
							}else{
								w_ijk=myK(mydist(c(tseq[j],tseq[k]),c(tseq[jj],tseq[kk]))/h)
							}
							upper=upper+w_ijk*(ehat_ij%o%ehat_ik)
							lower=lower+w_ijk
						}
					}
				}
				SSigma_e_hat1[(jj-1)*K+(1:K),(kk-1)*K+(1:K)]=upper/lower	
			}
		}
	}
	###make sure SSigma_e_hat1 is positive definite
	SSigma_e_hat1=nearPD(SSigma_e_hat1)$mat

	SSigma_e_hat2=kronecker(diag(rep(1,length(tseq))),Sigma_S-b4%*%(Sigma_IS)-t(b4%*%(Sigma_IS))+b4%*%Sigma_I%*%t(b4))
	###make positivity hold
	if (min(eigen(SSigma_e_hat1-nrow(Xmat)/(nrow(Xmat)-1)*SSigma_e_hat2)$values)>0){
		SSigma_e_hat=SSigma_e_hat1-SSigma_e_hat2
	}else{
		myf<-function(lambda){
			abs(min(eigen(SSigma_e_hat1-lambda*SSigma_e_hat2)$values))
		}
		lambda=optimize(myf,c(0,nrow(Xmat)/(nrow(Xmat)-1)))$minimum
		while(min(eigen(SSigma_e_hat1-lambda*SSigma_e_hat2)$values)<0){
			lambda=lambda/2
		}
		SSigma_e_hat=SSigma_e_hat1-(lambda-6/(nrow(Xmat)-1))*SSigma_e_hat2
	}	
        ###get initial value
	tmpmat=SSigma_e_hat[1:K,1:K]
	if (length(tseq)>1){
		for (tt in 2:length(tseq)){
			tmpmat=tmpmat+SSigma_e_hat[K*(tt-1)+(1:K),K*(tt-1)+(1:K)]	
		}
		tmpmat=tmpmat/length(tseq)
	}
	par_ini=c(0,0,mat2vec(tmpmat))

	###estimate parameter
	myobj<-function(par,type){
		rho_e=par[2]
		sigma_e=exp(par[1])
		Sigma_e=vec2mat(par[-c(1,2)])
		if (type=="AR"){
			Sigma_e_t=Sigma_AR(tseq,c(sigma_e,rho_e))
		}
		if (type=="CS"){
			Sigma_e_t=Sigma_CS(tseq,c(sigma_e,rho_e))
		}
		SSigma_e=kronecker(Sigma_e_t, Sigma_e)
		sum(log(eigen(SSigma_e)$value))+tr(solve(SSigma_e,SSigma_e_hat))
	}
	par_est=NA
	try({fit=optim(par_ini,myobj,type=type)
	par_est=fit$par})
	if (max(is.na(par_est))==1){
            ## estimate (sigma_e,rho_e), set Sigma_e as initial values
		myobj_rr<-function(par1,par2,type){
			rho_e=par1[2]
			sigma_e=exp(par1[1])
			Sigma_e=vec2mat(par2)
			if (type=="AR"){
				Sigma_e_t=Sigma_AR(tseq,c(sigma_e,rho_e))
			}
			if (type=="CS"){
				Sigma_e_t=Sigma_CS(tseq,c(sigma_e,rho_e))
			}
			SSigma_e=kronecker(Sigma_e_t, Sigma_e)
			sum(log(eigen(SSigma_e)$value))+tr(solve(SSigma_e,SSigma_e_hat))
		}
		try({fit=optim(par_ini[1:2],myobj_rr,type=type,par2=par_ini[-c(1:2)])
		par_est=c(fit$par,par_ini[-c(1:2)])})		
	}
	if (is.na(par_est)[1] || is.na(par_est)[2]){
		par_est=par_ini
	}
	###estimate covariance matrix
	rho_e=par_est[2]
	sigma_e=exp(par_est[1])
	Sigma_e=vec2mat(par_est[-c(1,2)])
	if (type=="AR"){Sigma_e_t=Sigma_AR(tseq,c(sigma_e,rho_e))}
	if (type=="CS"){Sigma_e_t=Sigma_CS(tseq,c(sigma_e,rho_e))}
	list(Sigma_e,Sigma_e_t)
}

###############################################################

myfit<-function(data,K,Q,p,tseq,Sigma_err,h=1,type,Sigma_e_t.true){	
	beta_ols=myfit_ols(data,K,Q,p)
	beta_star=myfit_oracle(data,K,Q,p)
	## first fit beta with no V 
      beta_ind=myfit_ind(data,K,Q,p,Sigma_err)
	## obtain V given beta_ind
      par=myest_par(data,K,Q,p,tseq,beta_ind,Sigma_err,h,type)
	Sigma_e=par[[1]]
	Sigma_e_t=par[[2]]
      ## after getting V, fit beta with V
      ## Oracle method using estimated working correlation
      beta_star_cor=myfit_oracle_cor(data,K,Q,p,tseq,Sigma_e,Sigma_e_t) 
      ## measurement correction method using estimated working correlation (one-iteration)
	beta_onestep=myfit_cor(data,K,Q,p,tseq,Sigma_err,Sigma_e,Sigma_e_t)
      par=myest_par(data,K,Q,p,tseq,beta_onestep,Sigma_err,h,type)
	Sigma_e=par[[1]]
	Sigma_e_t=par[[2]]
      ## after getting V, fit beta with V (second time)
	## measurement correction method using estimated working correlation (two iterations)
      beta_twostep=myfit_cor(data,K,Q,p,tseq,Sigma_err,Sigma_e,Sigma_e_t)
	## measurement correction method using true working correlation
      beta_truecor=myfit_cor(data,K,Q,p,tseq,Sigma_err,Sigma_e.true,Sigma_e_t.true)
	## dimension K rows, (2+p+Q)*betaNum columns, 
      return(c(beta_star,beta_star_cor,beta_ols,beta_ind,beta_truecor,beta_onestep,beta_twostep))
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
Nvec=c(200,400,800)
tseq=c(1:3)
seed=1111
sigma_a=0.2
###genearte a random positive definite matrix from Wishart distribution
#Sigma_a=sigma_a*rWishart(1,2*Q,diag(2*Q))[,,1]
Sigma_a=Sigma_CS(1:(2*Q),c(1,-0.1))
p=2
rho_X=0.1
Sigma_X=Sigma_CS(1:p,c(1,rho_X))

###Define correlation over time
# the correlation between different times for specific domain.
rho=0.3
## rho_eAR=c(0.3,0.5,0.8)
rho_eAR=c(0.5)
rho_eCS=c(0.2,0.3,0.5)
rho_err=0

####generate Sigma randomly from Wishart distribution with a scalar sigma
#sigma=5
#Sigma=diag(diag(sigma*rWishart(1,Q,diag(Q))[,,1]))
## 0.2=correlation between different imaging domains at a specific time.
Sigma=Sigma_AR(1:Q,c(3,0.2))

####generate Sigma_e randomly from Wishart distribution with a scalar sigma_e
#sigma_e=0.3
#Sigma_e=sigma_e*rWishart(1,K,diag(K))[,,1]
## K by K compound symmetric covariance matrix 
## with parameters $\sigma^2=2.5$ and $\rho_e=0.15$
Sigma_e=Sigma_CS(1:K,c(2.5,0.15))
Sigma_e.true=Sigma_e

####generate Sigma_err (include Sigma_I,Sigma_S,Sigma_IS) randomly from Wishart distribution with a scalar sigma_err
## SNR=10,4,2,1 sigma_err=c(0.13,0.34,0.66,1.33) (see sigma_err_Estimate.R)
SNR=c(10,4,2,1)
sigma_err=c(0.13,0.34,0.66,1.33)
Wis=rWishart(1,K+Q,diag(K+Q))[,,1]
#Sigma_err=diag(c(rep(1,K),rep(1.2,Q)))

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

# beta_ (oracle, oracle_cor, ols, ind, truecor, onestep, twostep)
betaNum=7

## within each b=1,...,B, do bootstrap M times
M=100
## M=2

###generate one fitting result
onefit<-function(seed,rho_e,Sigma_e_t.true,type){
	set.seed(seed)
      onefitresult=NULL
      for (sigma_errIndex in 1:length(sigma_err)){
      Sigma_err=sigma_err[sigma_errIndex]*Wis
      for (Nindex in 1:length(Nvec)){
	N=Nvec[Nindex]
      alldata=datagen(N,K,Q,p,tseq,beta0,beta1,beta2,b0,b1,b2,b3,Sigma_a,Sigma_X,Sigma,Sigma_e,Sigma_err,rho, rho_e,rho_err,gamma,seed,type)
	data0=alldata$data0
	sname=paste("S",as.character(1:K),sep="")
	sname1=paste("Sstar",as.character(1:K),sep="")
	iname=paste("I",as.character(1:Q),sep="")
	iname1=paste("Istar",as.character(1:Q),sep="")
	## cannot know gamma and Sigma_err in the true setting, so check the performance of
      ## estimator when plug in estimated gamma and Sigma from normal samples.
      gammahat=rep(NA,K)
	Sigma_errhat=matrix(data=NA,nrow=K+Q,ncol=K+Q)
      ### residual from normal N/2 group
	res=matrix(data=NA,nrow=nrow(data0),ncol=K+Q)
	for (k in 1:K){	
		yy=data0[,sname[k]]
		tt=data0$t
		id=data0$ID
		fit=lm(yy~tt+as.factor(id))
		gammahat[k]=fit$coef[2]
		res[,Q+k]=fit$res
	} 
	for (q in 1:Q){	
		yy=data0[,iname[q]]
		id=data0$ID
		res[,q]=lm(yy~as.factor(id))$res
	} 
      #### degree of freedom excluding fixed intercept (N intercepts)
      ## (N*n_i-N) is used to adjust for degrees of freedom used for estimating 
      ## the fixed intercept for the N individuals.
	Sigma_errhat=crossprod(res)/(nrow(data0)-length(unique(data0$ID)))
	data=alldata$data
	data[,sname]=data[,sname]-data$t%o%gammahat
	data[,sname1]=data[,sname1]-data$t%o%gammahat
	est=myfit(data,K,Q,p,tseq,Sigma_errhat,h=1,type,Sigma_e_t.true)

## for one sample of N subjects (N/2 normal), do M bootstrap
## the Bootstrap code within onefit function so that it report not only
## the estimates but also the estimated SE from Bootstrap (use 100 Bootstraps). 
## Within onefit function, after generate the data, you need add a loop to 
## resampale cases and controls to form new dataset. You don t need to compute 
## bias for Bootstrap, you just need to get SE.
## What you need to do is the following assuming Sample size N, Bootstrap time M, 
## Simulation time B
## Step 1: Simulation B samples each with sample size N cases and N/2 controls
## Step 2: For each simulated sample, get the estimate est[b]
## Step 3: For each simulated sample, get M bootstrap sample and 
## calculate the estimate est[b,m]. Then compute the sample variance over m 
## to obtain se[b].
## Step 4: For each simulated sample, compute whether the 95% confidence interval 
## will cover the truth or not cr[b]= as.numeric(abs(est[b]-true.est)/se[b])<1.96)
## Step 5: Compute bias as the average over b for (est[b]-est.true), 
## ESE as the average over b for se[b], SD is the variance over b for est[b], 
## CR is the average over b for cr[b]

      N=nrow(data)/length(tseq)
      Mest=array(data=NA,dim=c(K,betaNum*(2+p+Q),M))
      for (m in 1:M){
      Msample=sort(sample(1:N,replace=T))
      Msample0=sort(sample(1:(floor(N/3)),replace=T))
      Mdata0=Mdata=NULL
      for (i in 1:length(tseq)){
	Mdata0=rbind(Mdata0,data0[floor(N/3)*(i-1)+Msample0,])
        Mdata=rbind(Mdata,alldata$data[N*(i-1)+Msample,])
      }
      Mgammahat=rep(NA,K)
      MSigma_errhat=matrix(data=NA,nrow=K+Q,ncol=K+Q)
      ### residual from normal group N/2
	Mres=matrix(data=NA,nrow=nrow(Mdata0),ncol=K+Q)
	###redefine ID
	Mdata0$ID=rep(1:(N/3),length(tseq))
	Mdata$ID=rep(1:N,length(tseq))
	for (k in 1:K){	
		yy=Mdata0[,sname[k]]
		tt=Mdata0$t
		id=Mdata0$ID
		fit=lm(yy~tt+as.factor(id))
		Mgammahat[k]=fit$coef[2]
		Mres[,Q+k]=fit$res
	} 
	for (q in 1:Q){	
		yy=Mdata0[,iname[q]]
		id=Mdata0$ID
		Mres[,q]=lm(yy~as.factor(id))$res
	} 
      #### degree of freedom excluding fixed intercept (N intercepts)
	MSigma_errhat=crossprod(Mres)/(nrow(Mdata0)-nrow(Mdata0)/length(tseq))
	Mdata[,sname]=Mdata[,sname]-Mdata$t%o%Mgammahat
	Mdata[,sname1]=Mdata[,sname1]-Mdata$t%o%Mgammahat
	Mest[,,m]=myfit(Mdata,K,Q,p,tseq,MSigma_errhat,h=1,type,Sigma_e_t.true)
      }
      ## compute ESE
      est=matrix(data=est,nrow=K,ncol=betaNum*(2+p+Q),byrow=FALSE)
      Ese=apply(Mest,c(1,2),sd,na.rm=TRUE)
      ## CR (compute whether the 95% confidence interval will cover the truth or not)
      CR=abs(est-kronecker(matrix(rep(1,betaNum),nrow=1,ncol=betaNum),beta_true))/Ese<1.96      
      onefitresult=rbind(onefitresult,est,Ese,CR)   
     } # end Nindex
    } # end sigma_errIndex
     return(onefitresult)
    }

myresult_Sigma_e_t=function(rho_eVec,type){
for (rho_e in rho_eVec){
if (type=="AR") {Sigma_e_t.true=Sigma_AR(tseq,c(1,rho_e))}
if (type=="CS") {Sigma_e_t.true=Sigma_CS(tseq,c(1,rho_e))}	
#### repeat fitting B times (with B replications), and compute bias and SD
B=2
## B=2
result=array(data=NA,dim=c(K*3*length(Nvec)*length(SNR),betaNum*(2+p+Q),B))
## 8888+b or simseed+b
for (b in 1:B){try(result[,,b]<-onefit(simseed+b,rho_e,Sigma_e_t.true,type))}
bias=SD=Ese=CR=NULL
for (index in 1:(length(Nvec)*length(SNR))){
bias=rbind(bias,apply(result[(index-1)*3*K+c(1:K),,],c(1,2),mean,na.rm=TRUE)-kronecker(matrix(rep(1,betaNum),nrow=1,ncol=betaNum),beta_true))
SD=rbind(SD,apply(result[(index-1)*3*K+c(1:K),,],c(1,2),sd,na.rm=TRUE))
Ese=rbind(Ese,apply(result[(index-1)*3*K+K+c(1:K),,],c(1,2),mean,na.rm=TRUE))
CR=rbind(CR,apply(result[(index-1)*3*K+2*K+c(1:K),,],c(1,2),mean,na.rm=TRUE))
}

#### write out the results, (Bias, sd, Ese, CR)
result=data.frame(rbind(bias=round(bias,5),SD=round(SD,5),Ese=round(Ese,5),CR=round(CR,5)))
### rename column names
xname=paste("X",as.character(1:p),sep="")
iname=paste("I",as.character(1:Q),sep="")
inputname=c("intercept","t",xname,iname)
betaname=c("Oracle","Oraclecor","ols","Ind","truecor","onestep","twostep")
names(result)=paste(rep(betaname,each=length(inputname)),inputname,sep=" ")
## rename row names
Nname=paste("N=",Nvec,sep="")
sname=paste("S",as.character(1:K),sep="")
Nsname=paste(rep(Nname,each=length(sname)),sname,sep=" ")
SNRname=paste("SNR=",SNR,sep="")
SNRNsname=paste(rep(SNRname,each=length(Nsname)),Nsname,sep=" ")
outputname=c("Bias","SD","Ese","CR")
rownames(result)=paste(rep(outputname,each=length(SNRNsname)),SNRNsname,sep=" ")
write.csv(result,paste("result1",type,"rho_e=",rho_e,"datagenseed=",as.character(seed),"simseed=",as.character(simseed),".csv"))
} # end rho_e
} # end myresult_Sigma_e_t

myresult_Sigma_e_t(rho_eAR,type="AR")

###############################################################################################################################
###############################################################################################################################
###############################################################################################################################
## revision 2: change Sigma_err

library(MASS)
library(Matrix)

## Simseed These lines cannot be directly run in R, it need to be run using the submit file. Because the code
## R CMD BATCH "--args $SLURM_ARRAY_TASK_ID" Simulation_cor2_est_parallel.R is running with the args value assigned from the name 
## of SLURM_ARRAY_TASK_ID which vary among different tasks.

## args=commandArgs(trailingOnly=TRUE)
## simseed=8888+5*as.numeric(args[1])

####Simulation for step 1 models
####Data generation process

####Define variance process model here
Sigma_AR<-function(tseq,theta){
	Sigma=theta[1]*theta[2]^abs(outer(tseq,tseq,"-"))
	return(as.matrix(Sigma))
}
Sigma_CS<-function(tseq,theta){
	Sigma=theta[1]*theta[2]^(abs(outer(tseq,tseq,"-"))>0)
	return(as.matrix(Sigma))
}


###Data generation code
datagen<-function(N,K,Q,p,tseq,beta0,beta1,beta2,b0,b1,b2,b3,Sigma_a,Sigma_X,Sigma,Sigma_e,Sigma_err,rho, rho_e,rho_err,gamma,seed,type){
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
	if (type=="AR"){
            ###variance-covariance structure over time for imaging domains
		Sigma_t=Sigma_AR(tseq,c(1,rho))
		###variance-covariance structure over time for clinical domains
		Sigma_e_t=Sigma_AR(tseq,c(1,rho_e))
		###variance-covariance structure over time for measurement errors
		Sigma_err_t=Sigma_AR(tseq,c(1,rho_err))
	}
	if (type=="CS"){
            ###variance-covariance structure over time for imaging domains
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
    ## the learning effect is a linear function of how many times you get the test 
    ## (i.e., t). If you use constant, then you are assuming the skills you gained in 
    ## 2 tests are the same as the skills you gained in 3 tests.

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
	
####now generate normal population using half of the sample
	N=floor(N/2)
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
	data0=NULL
	###generate all residuals
	###assume kronecker structure for variance-covariance over time at different domains. Might need more complicated form for the transformed model
	SSigma_t=kronecker(Sigma_t, Sigma)
	epi_t_all=mvrnorm(N,mu=rep(0,nrow(SSigma_t)),Sigma=SSigma_t)
	SSigma_e_t=kronecker(Sigma_e_t, Sigma_e)
	e_t_all=mvrnorm(N,mu=rep(0,nrow(SSigma_e_t)),Sigma=SSigma_e_t)
	SSigma_err_t=kronecker(Sigma_err_t, 0.9*Sigma_err)
	err_t_all=mvrnorm(N,mu=rep(0,nrow(SSigma_err_t)),Sigma=SSigma_err_t)

	###generate data at each time
	for (i in 1:length(tseq)){
		t=tseq[i]
		###Generate imaging domains
		epi_t=epi_t_all[,(i-1)*Q+c(1:Q)]
		for (q in 1:Q){
			mu_I=beta0[q]+beta1[q]*0+X%*%beta2[q,]+a_i[,2*(q-1)+c(1:2)]%*%c(1,0)
			epi_t[,q]=mu_I
		}
		###Generate clinical domains (using I(t))
		e_t=e_t_all[,(i-1)*K+c(1:K)]
		for (k in 1:K){
            ## assume I(t) is constant over t for each individual
			mu_S=b0[k]+b1[k]*0+X%*%b2[k,]
			e_t[,k]=mu_S
		}
		e_t=e_t+rep(1,N)%o%(gamma*t)
		####add measurement error
		err_t=err_t_all[,(i-1)*(Q+K)+c(1:(Q+K))]
		tmp=cbind(1:N,t,X,cbind(epi_t,e_t)+err_t,epi_t,e_t)
		data0=rbind(data0,tmp)
	}
	###rename variables
	###use star to denote the truth here
	data0=data.frame(data0)
	names(data0)=c("ID","t",xname,iname,sname,iname1,sname1)	
	return(list(data=data,data0=data0))	
}  # end datagen

################# some useful functions

### l2 distance between 2 vectors
mydist<-function(u,v){
	sqrt(sum((u-v)^2))
}

### Box kernel
myK<-function(x){
	as.numeric(abs(x)<=1)/2
}

###trace
tr<-function(mat){
	return(sum(diag(mat)))
}

###mat to vec   lower triangle of matrix    
mat2vec<-function(mat){
	d=nrow(mat)
	v=NULL
	for (i in 1:d){
		for (j in 1:i){
			v=c(v,mat[i,j])
		}
	}
	v
}

###vec to mat
vec2mat<-function(v){
	d=floor(sqrt(2*length(v)))
	count=0
	mat=matrix(data=NA,nrow=d,ncol=d)
	for (i in 1:d){
		for (j in 1:i){
			count=count+1
			mat[i,j]=mat[j,i]=v[count]
		}
	}
	mat
}

########################################

###Model fitting code 
###OLS fitting code  (observed data) Xmat=t,X,I   Ymat=S,
## (independent working correlation over time), no measurement error correction 
myfit_ols<-function(data,K,Q,p){
	###OLS fitting
	Xmat=as.matrix(cbind(1,data[,c(2:(2+p+Q))]))
	Ymat=as.matrix(data[,2+p+Q+c(1:K)])
	YY=t(Xmat)%*%Ymat
	XX=t(Xmat)%*%Xmat
	beta_ols=t(solve(XX,YY))
	beta_ols
}
###Oracle fitting code, assuming error-free variables are observed
## independence working correlation
myfit_oracle<-function(data,K,Q,p){
	###oracle fitting  (true data Xmat=(t,X,I*)  Ymat=S*	
	Xmatstar=as.matrix(cbind(1,data[,c(2:(2+p),2+p+Q+K+c(1:Q))]))
	Ymatstar=as.matrix(data[,2+p+Q+K+Q+c(1:K)])
	YYstar=t(Xmatstar)%*%Ymatstar
	XXstar=t(Xmatstar)%*%Xmatstar
	beta_star=t(solve(XXstar,YYstar))
	beta_star
}

###Fitting code assuming independence working correlation, 
## includes measurement error adjustment
myfit_ind<-function(data,K,Q,p,Sigma_err){
	N=nrow(data)/length(tseq)
	###initial estimation 
	Sigma_I=Sigma_err[1:Q,1:Q]
	Sigma_S=Sigma_err[Q+c(1:K),Q+c(1:K)]
	Sigma_IS=Sigma_err[1:Q,Q+c(1:K)]
	###Measurement error correction
	Xmat=as.matrix(cbind(1,data[,c(2:(2+p+Q))]))
	Ymat=as.matrix(data[,2+p+Q+c(1:K)])
	YY=t(Xmat)%*%Ymat
	XX=t(Xmat)%*%Xmat
	####make sure positivity of XX, adjustment for measurement error, S has N*n_i elements
	## although XX is positive definite, after you subtract another positive definite 
      ## matrix lambda* nrow(Xmat)*Sigma_I from it, it can become negative, so the minimum 
      ## eigen can be negative. Here using the absolute value is for the purpose of solving 
      ## the lambda such that the minimum eigen value equals 0.
      
      YY[2+p+c(1:Q),1:K]=YY[2+p+c(1:Q),1:K]-nrow(Xmat)*Sigma_IS
	if (min(eigen(XX[2+p+c(1:Q),2+p+c(1:Q)]-nrow(Xmat)/(nrow(Xmat)-1)*nrow(Xmat)*Sigma_I)$values)>0){
		XX[2+p+c(1:Q),2+p+c(1:Q)]=XX[2+p+c(1:Q),2+p+c(1:Q)]-nrow(Xmat)*Sigma_I
	}else{
		myf<-function(lambda){
			abs(min(eigen(XX[2+p+c(1:Q),2+p+c(1:Q)]-lambda*nrow(Xmat)*Sigma_I)$values))
		}
		lambda=optimize(myf,c(0,nrow(Xmat)/(nrow(Xmat)-1)))$minimum
		XX[2+p+c(1:Q),2+p+c(1:Q)]=XX[2+p+c(1:Q),2+p+c(1:Q)]-(lambda-6/(nrow(Xmat)-1))*nrow(Xmat)*Sigma_I
	}	
	beta_ind=t(solve(XX,YY))
	beta_ind
}

###Oracle fitting code, assuming error-free variables are observed;
###Fitting code assuming specific working correlation matrix 
## either true or one-step working correlation Sigma_e, Sigma_e_t
myfit_oracle_cor<-function(data,K,Q,p,tseq,Sigma_e,Sigma_e_t){
	N=nrow(data)/length(tseq)
	##oracle fitting  (true data Xmatstar=(t,X,I*)  Ymatstar=S*	
	Xmatstar=as.matrix(cbind(1,data[,c(2:(2+p),2+p+Q+K+c(1:Q))]))
	Ymatstar=as.matrix(data[,2+p+Q+K+Q+c(1:K)])
	###estimate parameter
	WW=solve(Sigma_e_t)
	Hhat=matrix(data=0,nrow=2+p+Q,ncol=2+p+Q)
	for (j in 1:length(tseq)){
		for (k in 1:length(tseq)){
			Hhat=Hhat+c(WW[j,k])*t(Xmatstar[1:N+(j-1)*N,])%*%Xmatstar[1:N+(k-1)*N,]
		}
	}
	HHhat=kronecker(solve(Hhat),Sigma_e)
	###update U
	SSigma_e=kronecker(Sigma_e_t, Sigma_e)
	WWtilde=solve(SSigma_e)
	Uhat=matrix(data=0,nrow=(2+p+Q)*K,ncol=1)
	for (i in 1:N){
		for (j in 1:(length(tseq)*K)){
			for (l in 1:(length(tseq)*K)){
				jj=(j-1)%/%K+1
				ll=(l-1)%/%K+1
				jjj=(j-1)%%K+1
				lll=(l-1)%%K+1
				Uhat=Uhat+c(WWtilde[j,l])*kronecker(Xmatstar[i+((1:length(tseq))-1)*N,],diag(rep(1,K)))[j,]*Ymatstar[i+(ll-1)*N,lll]
			}
		}
	}
	beta_star_cor=matrix(data=HHhat%*%Uhat,nrow=K,ncol=2+p+Q,byrow=FALSE)
	beta_star_cor
}

###Fitting code assuming specific working correlation matrix, 
## includes measurement error adjustment
## either true or one-step working correlation Sigma_e, Sigma_e_t
myfit_cor<-function(data,K,Q,p,tseq,Sigma_err,Sigma_e,Sigma_e_t){
	N=nrow(data)/length(tseq)
	###initial estimation 
	Sigma_I=Sigma_err[1:Q,1:Q]
	Sigma_S=Sigma_err[Q+c(1:K),Q+c(1:K)]
	Sigma_IS=Sigma_err[1:Q,Q+c(1:K)]
	###Measurement error correction
	Xmat=as.matrix(cbind(1,data[,c(2:(2+p+Q))]))
	Ymat=as.matrix(data[,2+p+Q+c(1:K)])
	YY=t(Xmat)%*%Ymat
	XX=t(Xmat)%*%Xmat
	###estimate parameter
	WW=solve(Sigma_e_t)
	Hhat1=matrix(data=0,nrow=2+p+Q,ncol=2+p+Q)
	###expanded sigma_I
	Omega_I=matrix(data=0,nrow=2+p+Q,ncol=2+p+Q)
	Omega_I[2+p+(1:Q),2+p+(1:Q)]=Sigma_I
	Hhat2=N*tr(WW)*Omega_I
	for (j in 1:length(tseq)){
		for (k in 1:length(tseq)){
			Hhat1=Hhat1+c(WW[j,k])*t(Xmat[1:N+(j-1)*N,])%*%Xmat[1:N+(k-1)*N,]
		}
	}
	###need make sure positivity
	if (min(eigen(Hhat1-nrow(Xmat)/(nrow(Xmat)-1)*Hhat2)$values)>0){
		Hhat=Hhat1-Hhat2
	}else{
		myf<-function(lambda){
			abs(min(eigen(Hhat1-lambda*Hhat2)$values))
		}
		lambda=optimize(myf,c(0,nrow(Xmat)/(nrow(Xmat)-1)))$minimum
		Hhat=Hhat1-(lambda-6/(nrow(Xmat)-1))*Hhat2
	}	
	HHhat=kronecker(solve(Hhat),Sigma_e)
	
	###update U
	SSigma_e=kronecker(Sigma_e_t, Sigma_e)
	WWtilde=solve(SSigma_e)
	Omega_IS=matrix(data=0,nrow=2+p+Q,ncol=K)
	Omega_IS[2+p+(1:Q),]=Sigma_IS
	Uhat=matrix(data=0,nrow=(2+p+Q)*K,ncol=1)
	for (i in 1:N){
		for (j in 1:(length(tseq)*K)){
			for (l in 1:(length(tseq)*K)){
				jj=(j-1)%/%K+1
				ll=(l-1)%/%K+1
				jjj=(j-1)%%K+1
				lll=(l-1)%%K+1
				Uhat=Uhat+c(WWtilde[j,l])*kronecker(Xmat[i+((1:length(tseq))-1)*N,],diag(rep(1,K)))[j,]*Ymat[i+(ll-1)*N,lll]
				if (jj==ll){
					e_lll=rep(0,K)
					e_lll[lll]=1
					e_jjj=rep(0,K)
					e_jjj[jjj]=1
					Uhat=Uhat-c(WWtilde[j,l])*kronecker(Omega_IS%*%e_lll,e_jjj)
				}
			}
		}
	}
	beta_onestep=matrix(data=HHhat%*%Uhat,nrow=K,ncol=2+p+Q,byrow=FALSE)
	beta_onestep
}

###code to estimate parameter related to correlation structure: Sigma_e,Sigma_e_t
## myest_par, that Kernel estimated variance-covariance matrix SSigma_e_hat1 
## is not always positive definite when measurement error is large 
## and sometimes SSigma_e_hat are ill-conditioned

myest_par<-function(data,K,Q,p,tseq,beta_ind,Sigma_err,h,type){
	N=nrow(data)/length(tseq)
	###initial estimation 
	Sigma_I=Sigma_err[1:Q,1:Q]
	Sigma_S=Sigma_err[Q+c(1:K),Q+c(1:K)]
	Sigma_IS=Sigma_err[1:Q,Q+c(1:K)]

	###use one step update, estimate rho_e, Sigma_e
	###get empirical estimation for SSigma_e
	SSigma_e_hat1=matrix(data=0,nrow=K*length(tseq),ncol=K*length(tseq))
	b4=beta_ind[,2+p+c(1:Q)]	
	Xmat=as.matrix(cbind(1,data[,c(2:(2+p+Q))]))
	Ymat=as.matrix(data[,2+p+Q+c(1:K)])
	Ehat=Ymat-Xmat%*%t(beta_ind)
	for (jj in 1:length(tseq)){
		for (kk in 1:length(tseq)){
			upper=matrix(data=0,nrow=K,ncol=K)
			lower=0
			if (jj==kk){
				for (i in 1:N){
					for (j in 1:length(tseq)){
						ehat_ij=Ehat[i+(j-1)*N,]
						w_ij=myK(abs(tseq[j]-tseq[jj])/h)
						upper=upper+w_ij*(ehat_ij%o%ehat_ij)
						lower=lower+w_ij
					}
				}
				SSigma_e_hat1[(jj-1)*K+(1:K),(kk-1)*K+(1:K)]=upper/lower				
			}
			if (jj!=kk){
				for (i in 1:N){
					for (j in 1:length(tseq)){
						for (k in 1:length(tseq)){
							ehat_ij=Ehat[i+(j-1)*N,]
							ehat_ik=Ehat[i+(k-1)*N,]
							if (k==j){
								w_ijk=0
							}else{
								w_ijk=myK(mydist(c(tseq[j],tseq[k]),c(tseq[jj],tseq[kk]))/h)
							}
							upper=upper+w_ijk*(ehat_ij%o%ehat_ik)
							lower=lower+w_ijk
						}
					}
				}
				SSigma_e_hat1[(jj-1)*K+(1:K),(kk-1)*K+(1:K)]=upper/lower	
			}
		}
	}
	###make sure SSigma_e_hat1 is positive definite
	SSigma_e_hat1=nearPD(SSigma_e_hat1)$mat

	SSigma_e_hat2=kronecker(diag(rep(1,length(tseq))),Sigma_S-b4%*%(Sigma_IS)-t(b4%*%(Sigma_IS))+b4%*%Sigma_I%*%t(b4))
	###make positivity hold
	if (min(eigen(SSigma_e_hat1-nrow(Xmat)/(nrow(Xmat)-1)*SSigma_e_hat2)$values)>0){
		SSigma_e_hat=SSigma_e_hat1-SSigma_e_hat2
	}else{
		myf<-function(lambda){
			abs(min(eigen(SSigma_e_hat1-lambda*SSigma_e_hat2)$values))
		}
		lambda=optimize(myf,c(0,nrow(Xmat)/(nrow(Xmat)-1)))$minimum
		while(min(eigen(SSigma_e_hat1-lambda*SSigma_e_hat2)$values)<0){
			lambda=lambda/2
		}
		SSigma_e_hat=SSigma_e_hat1-(lambda-6/(nrow(Xmat)-1))*SSigma_e_hat2
	}	
        ###get initial value
	tmpmat=SSigma_e_hat[1:K,1:K]
	if (length(tseq)>1){
		for (tt in 2:length(tseq)){
			tmpmat=tmpmat+SSigma_e_hat[K*(tt-1)+(1:K),K*(tt-1)+(1:K)]	
		}
		tmpmat=tmpmat/length(tseq)
	}
	par_ini=c(0,0,mat2vec(tmpmat))

	###estimate parameter
	myobj<-function(par,type){
		rho_e=par[2]
		sigma_e=exp(par[1])
		Sigma_e=vec2mat(par[-c(1,2)])
		if (type=="AR"){
			Sigma_e_t=Sigma_AR(tseq,c(sigma_e,rho_e))
		}
		if (type=="CS"){
			Sigma_e_t=Sigma_CS(tseq,c(sigma_e,rho_e))
		}
		SSigma_e=kronecker(Sigma_e_t, Sigma_e)
		sum(log(eigen(SSigma_e)$value))+tr(solve(SSigma_e,SSigma_e_hat))
	}
	par_est=NA
	try({fit=optim(par_ini,myobj,type=type)
	par_est=fit$par})
	if (max(is.na(par_est))==1){
            ## estimate (sigma_e,rho_e), set Sigma_e as initial values
		myobj_rr<-function(par1,par2,type){
			rho_e=par1[2]
			sigma_e=exp(par1[1])
			Sigma_e=vec2mat(par2)
			if (type=="AR"){
				Sigma_e_t=Sigma_AR(tseq,c(sigma_e,rho_e))
			}
			if (type=="CS"){
				Sigma_e_t=Sigma_CS(tseq,c(sigma_e,rho_e))
			}
			SSigma_e=kronecker(Sigma_e_t, Sigma_e)
			sum(log(eigen(SSigma_e)$value))+tr(solve(SSigma_e,SSigma_e_hat))
		}
		try({fit=optim(par_ini[1:2],myobj_rr,type=type,par2=par_ini[-c(1:2)])
		par_est=c(fit$par,par_ini[-c(1:2)])})		
	}
	if (is.na(par_est)[1] || is.na(par_est)[2]){
		par_est=par_ini
	}
	###estimate covariance matrix
	rho_e=par_est[2]
	sigma_e=exp(par_est[1])
	Sigma_e=vec2mat(par_est[-c(1,2)])
	if (type=="AR"){Sigma_e_t=Sigma_AR(tseq,c(sigma_e,rho_e))}
	if (type=="CS"){Sigma_e_t=Sigma_CS(tseq,c(sigma_e,rho_e))}
	list(Sigma_e,Sigma_e_t)
}

###############################################################

myfit<-function(data,K,Q,p,tseq,Sigma_err,h=1,type,Sigma_e_t.true){	
	beta_ols=myfit_ols(data,K,Q,p)
	beta_star=myfit_oracle(data,K,Q,p)
	## first fit beta with no V 
      beta_ind=myfit_ind(data,K,Q,p,Sigma_err)
	## obtain V given beta_ind
      par=myest_par(data,K,Q,p,tseq,beta_ind,Sigma_err,h,type)
	Sigma_e=par[[1]]
	Sigma_e_t=par[[2]]
      ## after getting V, fit beta with V
      ## Oracle method using estimated working correlation
      beta_star_cor=myfit_oracle_cor(data,K,Q,p,tseq,Sigma_e,Sigma_e_t) 
      ## measurement correction method using estimated working correlation (one-iteration)
	beta_onestep=myfit_cor(data,K,Q,p,tseq,Sigma_err,Sigma_e,Sigma_e_t)
      par=myest_par(data,K,Q,p,tseq,beta_onestep,Sigma_err,h,type)
	Sigma_e=par[[1]]
	Sigma_e_t=par[[2]]
      ## after getting V, fit beta with V (second time)
	## measurement correction method using estimated working correlation (two iterations)
      beta_twostep=myfit_cor(data,K,Q,p,tseq,Sigma_err,Sigma_e,Sigma_e_t)
	## measurement correction method using true working correlation
      beta_truecor=myfit_cor(data,K,Q,p,tseq,Sigma_err,Sigma_e.true,Sigma_e_t.true)
	## dimension K rows, (2+p+Q)*betaNum columns, 
      return(c(beta_star,beta_star_cor,beta_ols,beta_ind,beta_truecor,beta_onestep,beta_twostep))
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
Nvec=c(200,400,800)
tseq=c(1:3)
seed=1111
sigma_a=0.2
###genearte a random positive definite matrix from Wishart distribution
#Sigma_a=sigma_a*rWishart(1,2*Q,diag(2*Q))[,,1]
Sigma_a=Sigma_CS(1:(2*Q),c(1,-0.1))
p=2
rho_X=0.1
Sigma_X=Sigma_CS(1:p,c(1,rho_X))

###Define correlation over time
# the correlation between different times for specific domain.
rho=0.3
## rho_eAR=c(0.3,0.5,0.8)
rho_eAR=c(0.5)
rho_eCS=c(0.2,0.3,0.5)
rho_err=0

####generate Sigma randomly from Wishart distribution with a scalar sigma
#sigma=5
#Sigma=diag(diag(sigma*rWishart(1,Q,diag(Q))[,,1]))
## 0.2=correlation between different imaging domains at a specific time.
Sigma=Sigma_AR(1:Q,c(3,0.2))

####generate Sigma_e randomly from Wishart distribution with a scalar sigma_e
#sigma_e=0.3
#Sigma_e=sigma_e*rWishart(1,K,diag(K))[,,1]
## K by K compound symmetric covariance matrix 
## with parameters $\sigma^2=2.5$ and $\rho_e=0.15$
Sigma_e=Sigma_CS(1:K,c(2.5,0.15))
Sigma_e.true=Sigma_e

####generate Sigma_err (include Sigma_I,Sigma_S,Sigma_IS) randomly from Wishart distribution with a scalar sigma_err
## SNR=10,4,2,1 sigma_err=c(0.13,0.34,0.66,1.33) (see sigma_err_Estimate.R)
SNR=c(10,4,2,1)
sigma_err=c(0.13,0.34,0.66,1.33)
Wis=rWishart(1,K+Q,diag(K+Q))[,,1]
#Sigma_err=diag(c(rep(1,K),rep(1.2,Q)))

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

# beta_ (oracle, oracle_cor, ols, ind, truecor, onestep, twostep)
betaNum=7

## within each b=1,...,B, do bootstrap M times
M=100
## M=2

###generate one fitting result
onefit<-function(seed,rho_e,Sigma_e_t.true,type){
	set.seed(seed)
      onefitresult=NULL
      for (sigma_errIndex in 1:length(sigma_err)){
      Sigma_err=sigma_err[sigma_errIndex]*Wis
      for (Nindex in 1:length(Nvec)){
	N=Nvec[Nindex]
      alldata=datagen(N,K,Q,p,tseq,beta0,beta1,beta2,b0,b1,b2,b3,Sigma_a,Sigma_X,Sigma,Sigma_e,Sigma_err,rho, rho_e,rho_err,gamma,seed,type)
	data0=alldata$data0
	sname=paste("S",as.character(1:K),sep="")
	sname1=paste("Sstar",as.character(1:K),sep="")
	iname=paste("I",as.character(1:Q),sep="")
	iname1=paste("Istar",as.character(1:Q),sep="")
	## cannot know gamma and Sigma_err in the true setting, so check the performance of
      ## estimator when plug in estimated gamma and Sigma from normal samples.
      gammahat=rep(NA,K)
	Sigma_errhat=matrix(data=NA,nrow=K+Q,ncol=K+Q)
      ### residual from normal N/2 group
	res=matrix(data=NA,nrow=nrow(data0),ncol=K+Q)
	for (k in 1:K){	
		yy=data0[,sname[k]]
		tt=data0$t
		id=data0$ID
		fit=lm(yy~tt+as.factor(id))
		gammahat[k]=fit$coef[2]
		res[,Q+k]=fit$res
	} 
	for (q in 1:Q){	
		yy=data0[,iname[q]]
		id=data0$ID
		res[,q]=lm(yy~as.factor(id))$res
	} 
      #### degree of freedom excluding fixed intercept (N intercepts)
      ## (N*n_i-N) is used to adjust for degrees of freedom used for estimating 
      ## the fixed intercept for the N individuals.
	Sigma_errhat=crossprod(res)/(nrow(data0)-length(unique(data0$ID)))
	data=alldata$data
	data[,sname]=data[,sname]-data$t%o%gammahat
	data[,sname1]=data[,sname1]-data$t%o%gammahat
	est=myfit(data,K,Q,p,tseq,Sigma_errhat,h=1,type,Sigma_e_t.true)

## for one sample of N subjects (N/2 normal), do M bootstrap
## the Bootstrap code within onefit function so that it report not only
## the estimates but also the estimated SE from Bootstrap (use 100 Bootstraps). 
## Within onefit function, after generate the data, you need add a loop to 
## resampale cases and controls to form new dataset. You don t need to compute 
## bias for Bootstrap, you just need to get SE.
## What you need to do is the following assuming Sample size N, Bootstrap time M, 
## Simulation time B
## Step 1: Simulation B samples each with sample size N cases and N/2 controls
## Step 2: For each simulated sample, get the estimate est[b]
## Step 3: For each simulated sample, get M bootstrap sample and 
## calculate the estimate est[b,m]. Then compute the sample variance over m 
## to obtain se[b].
## Step 4: For each simulated sample, compute whether the 95% confidence interval 
## will cover the truth or not cr[b]= as.numeric(abs(est[b]-true.est)/se[b])<1.96)
## Step 5: Compute bias as the average over b for (est[b]-est.true), 
## ESE as the average over b for se[b], SD is the variance over b for est[b], 
## CR is the average over b for cr[b]

      N=nrow(data)/length(tseq)
      Mest=array(data=NA,dim=c(K,betaNum*(2+p+Q),M))
      for (m in 1:M){
      Msample=sort(sample(1:N,replace=T))
      Msample0=sort(sample(1:(N/2),replace=T))
      Mdata0=Mdata=NULL
      for (i in 1:length(tseq)){
	Mdata0=rbind(Mdata0,data0[N/2*(i-1)+Msample0,])
        Mdata=rbind(Mdata,alldata$data[N*(i-1)+Msample,])
      }
      Mgammahat=rep(NA,K)
      MSigma_errhat=matrix(data=NA,nrow=K+Q,ncol=K+Q)
      ### residual from normal group N/2
	Mres=matrix(data=NA,nrow=nrow(Mdata0),ncol=K+Q)
	###redefine ID
	Mdata0$ID=rep(1:(N/2),length(tseq))
	Mdata$ID=rep(1:N,length(tseq))
	for (k in 1:K){	
		yy=Mdata0[,sname[k]]
		tt=Mdata0$t
		id=Mdata0$ID
		fit=lm(yy~tt+as.factor(id))
		Mgammahat[k]=fit$coef[2]
		Mres[,Q+k]=fit$res
	} 
	for (q in 1:Q){	
		yy=Mdata0[,iname[q]]
		id=Mdata0$ID
		Mres[,q]=lm(yy~as.factor(id))$res
	} 
      #### degree of freedom excluding fixed intercept (N intercepts)
	MSigma_errhat=crossprod(Mres)/(nrow(Mdata0)-nrow(Mdata0)/length(tseq))
	Mdata[,sname]=Mdata[,sname]-Mdata$t%o%Mgammahat
	Mdata[,sname1]=Mdata[,sname1]-Mdata$t%o%Mgammahat
	Mest[,,m]=myfit(Mdata,K,Q,p,tseq,MSigma_errhat,h=1,type,Sigma_e_t.true)
      }
      ## compute ESE
      est=matrix(data=est,nrow=K,ncol=betaNum*(2+p+Q),byrow=FALSE)
      Ese=apply(Mest,c(1,2),sd,na.rm=TRUE)
      ## CR (compute whether the 95% confidence interval will cover the truth or not)
      CR=abs(est-kronecker(matrix(rep(1,betaNum),nrow=1,ncol=betaNum),beta_true))/Ese<1.96      
      onefitresult=rbind(onefitresult,est,Ese,CR)   
     } # end Nindex
    } # end sigma_errIndex
     return(onefitresult)
    }

myresult_Sigma_e_t=function(rho_eVec,type){
for (rho_e in rho_eVec){
if (type=="AR") {Sigma_e_t.true=Sigma_AR(tseq,c(1,rho_e))}
if (type=="CS") {Sigma_e_t.true=Sigma_CS(tseq,c(1,rho_e))}	
#### repeat fitting B times (with B replications), and compute bias and SD
B=2
## B=2
result=array(data=NA,dim=c(K*3*length(Nvec)*length(SNR),betaNum*(2+p+Q),B))
## 8888+b or simseed+b
for (b in 1:B){try(result[,,b]<-onefit(simseed+b,rho_e,Sigma_e_t.true,type))}
bias=SD=Ese=CR=NULL
for (index in 1:(length(Nvec)*length(SNR))){
bias=rbind(bias,apply(result[(index-1)*3*K+c(1:K),,],c(1,2),mean,na.rm=TRUE)-kronecker(matrix(rep(1,betaNum),nrow=1,ncol=betaNum),beta_true))
SD=rbind(SD,apply(result[(index-1)*3*K+c(1:K),,],c(1,2),sd,na.rm=TRUE))
Ese=rbind(Ese,apply(result[(index-1)*3*K+K+c(1:K),,],c(1,2),mean,na.rm=TRUE))
CR=rbind(CR,apply(result[(index-1)*3*K+2*K+c(1:K),,],c(1,2),mean,na.rm=TRUE))
}

#### write out the results, (Bias, sd, Ese, CR)
result=data.frame(rbind(bias=round(bias,5),SD=round(SD,5),Ese=round(Ese,5),CR=round(CR,5)))
### rename column names
xname=paste("X",as.character(1:p),sep="")
iname=paste("I",as.character(1:Q),sep="")
inputname=c("intercept","t",xname,iname)
betaname=c("Oracle","Oraclecor","ols","Ind","truecor","onestep","twostep")
names(result)=paste(rep(betaname,each=length(inputname)),inputname,sep=" ")
## rename row names
Nname=paste("N=",Nvec,sep="")
sname=paste("S",as.character(1:K),sep="")
Nsname=paste(rep(Nname,each=length(sname)),sname,sep=" ")
SNRname=paste("SNR=",SNR,sep="")
SNRNsname=paste(rep(SNRname,each=length(Nsname)),Nsname,sep=" ")
outputname=c("Bias","SD","Ese","CR")
rownames(result)=paste(rep(outputname,each=length(SNRNsname)),SNRNsname,sep=" ")
write.csv(result,paste("result2",type,"rho_e=",rho_e,"datagenseed=",as.character(seed),"simseed=",as.character(simseed),".csv"))
} # end rho_e
} # end myresult_Sigma_e_t

myresult_Sigma_e_t(rho_eAR,type="AR")

###############################################################################################################################
###############################################################################################################################
###############################################################################################################################
## revision 3: change gamma

library(MASS)
library(Matrix)

## Simseed These lines cannot be directly run in R, it need to be run using the submit file. Because the code
## R CMD BATCH "--args $SLURM_ARRAY_TASK_ID" Simulation_cor2_est_parallel.R is running with the args value assigned from the name 
## of SLURM_ARRAY_TASK_ID which vary among different tasks.

## args=commandArgs(trailingOnly=TRUE)
## simseed=8888+5*as.numeric(args[1])

####Simulation for step 1 models
####Data generation process

####Define variance process model here
Sigma_AR<-function(tseq,theta){
	Sigma=theta[1]*theta[2]^abs(outer(tseq,tseq,"-"))
	return(as.matrix(Sigma))
}
Sigma_CS<-function(tseq,theta){
	Sigma=theta[1]*theta[2]^(abs(outer(tseq,tseq,"-"))>0)
	return(as.matrix(Sigma))
}


###Data generation code
datagen<-function(N,K,Q,p,tseq,beta0,beta1,beta2,b0,b1,b2,b3,Sigma_a,Sigma_X,Sigma,Sigma_e,Sigma_err,rho, rho_e,rho_err,gamma,seed,type){
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
	if (type=="AR"){
            ###variance-covariance structure over time for imaging domains
		Sigma_t=Sigma_AR(tseq,c(1,rho))
		###variance-covariance structure over time for clinical domains
		Sigma_e_t=Sigma_AR(tseq,c(1,rho_e))
		###variance-covariance structure over time for measurement errors
		Sigma_err_t=Sigma_AR(tseq,c(1,rho_err))
	}
	if (type=="CS"){
            ###variance-covariance structure over time for imaging domains
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
    ## the learning effect is a linear function of how many times you get the test 
    ## (i.e., t). If you use constant, then you are assuming the skills you gained in 
    ## 2 tests are the same as the skills you gained in 3 tests.

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
	
####now generate normal population using half of the sample
	N=floor(N/2)
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
	data0=NULL
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
			mu_I=beta0[q]+beta1[q]*0+X%*%beta2[q,]+a_i[,2*(q-1)+c(1:2)]%*%c(1,0)
			epi_t[,q]=mu_I
		}
		###Generate clinical domains (using I(t))
		e_t=e_t_all[,(i-1)*K+c(1:K)]
		for (k in 1:K){
            ## assume I(t) is constant over t for each individual
			mu_S=b0[k]+b1[k]*0+X%*%b2[k,]
			e_t[,k]=mu_S
		}
		e_t=e_t+rep(1,N)%o%(1.1*gamma*t)
		####add measurement error
		err_t=err_t_all[,(i-1)*(Q+K)+c(1:(Q+K))]
		tmp=cbind(1:N,t,X,cbind(epi_t,e_t)+err_t,epi_t,e_t)
		data0=rbind(data0,tmp)
	}
	###rename variables
	###use star to denote the truth here
	data0=data.frame(data0)
	names(data0)=c("ID","t",xname,iname,sname,iname1,sname1)	
	return(list(data=data,data0=data0))	
}  # end datagen

################# some useful functions

### l2 distance between 2 vectors
mydist<-function(u,v){
	sqrt(sum((u-v)^2))
}

### Box kernel
myK<-function(x){
	as.numeric(abs(x)<=1)/2
}

###trace
tr<-function(mat){
	return(sum(diag(mat)))
}

###mat to vec   lower triangle of matrix    
mat2vec<-function(mat){
	d=nrow(mat)
	v=NULL
	for (i in 1:d){
		for (j in 1:i){
			v=c(v,mat[i,j])
		}
	}
	v
}

###vec to mat
vec2mat<-function(v){
	d=floor(sqrt(2*length(v)))
	count=0
	mat=matrix(data=NA,nrow=d,ncol=d)
	for (i in 1:d){
		for (j in 1:i){
			count=count+1
			mat[i,j]=mat[j,i]=v[count]
		}
	}
	mat
}

########################################

###Model fitting code 
###OLS fitting code  (observed data) Xmat=t,X,I   Ymat=S,
## (independent working correlation over time), no measurement error correction 
myfit_ols<-function(data,K,Q,p){
	###OLS fitting
	Xmat=as.matrix(cbind(1,data[,c(2:(2+p+Q))]))
	Ymat=as.matrix(data[,2+p+Q+c(1:K)])
	YY=t(Xmat)%*%Ymat
	XX=t(Xmat)%*%Xmat
	beta_ols=t(solve(XX,YY))
	beta_ols
}
###Oracle fitting code, assuming error-free variables are observed
## independence working correlation
myfit_oracle<-function(data,K,Q,p){
	###oracle fitting  (true data Xmat=(t,X,I*)  Ymat=S*	
	Xmatstar=as.matrix(cbind(1,data[,c(2:(2+p),2+p+Q+K+c(1:Q))]))
	Ymatstar=as.matrix(data[,2+p+Q+K+Q+c(1:K)])
	YYstar=t(Xmatstar)%*%Ymatstar
	XXstar=t(Xmatstar)%*%Xmatstar
	beta_star=t(solve(XXstar,YYstar))
	beta_star
}

###Fitting code assuming independence working correlation, 
## includes measurement error adjustment
myfit_ind<-function(data,K,Q,p,Sigma_err){
	N=nrow(data)/length(tseq)
	###initial estimation 
	Sigma_I=Sigma_err[1:Q,1:Q]
	Sigma_S=Sigma_err[Q+c(1:K),Q+c(1:K)]
	Sigma_IS=Sigma_err[1:Q,Q+c(1:K)]
	###Measurement error correction
	Xmat=as.matrix(cbind(1,data[,c(2:(2+p+Q))]))
	Ymat=as.matrix(data[,2+p+Q+c(1:K)])
	YY=t(Xmat)%*%Ymat
	XX=t(Xmat)%*%Xmat
	####make sure positivity of XX, adjustment for measurement error, S has N*n_i elements
	## although XX is positive definite, after you subtract another positive definite 
      ## matrix lambda* nrow(Xmat)*Sigma_I from it, it can become negative, so the minimum 
      ## eigen can be negative. Here using the absolute value is for the purpose of solving 
      ## the lambda such that the minimum eigen value equals 0.
      
      YY[2+p+c(1:Q),1:K]=YY[2+p+c(1:Q),1:K]-nrow(Xmat)*Sigma_IS
	if (min(eigen(XX[2+p+c(1:Q),2+p+c(1:Q)]-nrow(Xmat)/(nrow(Xmat)-1)*nrow(Xmat)*Sigma_I)$values)>0){
		XX[2+p+c(1:Q),2+p+c(1:Q)]=XX[2+p+c(1:Q),2+p+c(1:Q)]-nrow(Xmat)*Sigma_I
	}else{
		myf<-function(lambda){
			abs(min(eigen(XX[2+p+c(1:Q),2+p+c(1:Q)]-lambda*nrow(Xmat)*Sigma_I)$values))
		}
		lambda=optimize(myf,c(0,nrow(Xmat)/(nrow(Xmat)-1)))$minimum
		XX[2+p+c(1:Q),2+p+c(1:Q)]=XX[2+p+c(1:Q),2+p+c(1:Q)]-(lambda-6/(nrow(Xmat)-1))*nrow(Xmat)*Sigma_I
	}	
	beta_ind=t(solve(XX,YY))
	beta_ind
}

###Oracle fitting code, assuming error-free variables are observed;
###Fitting code assuming specific working correlation matrix 
## either true or one-step working correlation Sigma_e, Sigma_e_t
myfit_oracle_cor<-function(data,K,Q,p,tseq,Sigma_e,Sigma_e_t){
	N=nrow(data)/length(tseq)
	##oracle fitting  (true data Xmatstar=(t,X,I*)  Ymatstar=S*	
	Xmatstar=as.matrix(cbind(1,data[,c(2:(2+p),2+p+Q+K+c(1:Q))]))
	Ymatstar=as.matrix(data[,2+p+Q+K+Q+c(1:K)])
	###estimate parameter
	WW=solve(Sigma_e_t)
	Hhat=matrix(data=0,nrow=2+p+Q,ncol=2+p+Q)
	for (j in 1:length(tseq)){
		for (k in 1:length(tseq)){
			Hhat=Hhat+c(WW[j,k])*t(Xmatstar[1:N+(j-1)*N,])%*%Xmatstar[1:N+(k-1)*N,]
		}
	}
	HHhat=kronecker(solve(Hhat),Sigma_e)
	###update U
	SSigma_e=kronecker(Sigma_e_t, Sigma_e)
	WWtilde=solve(SSigma_e)
	Uhat=matrix(data=0,nrow=(2+p+Q)*K,ncol=1)
	for (i in 1:N){
		for (j in 1:(length(tseq)*K)){
			for (l in 1:(length(tseq)*K)){
				jj=(j-1)%/%K+1
				ll=(l-1)%/%K+1
				jjj=(j-1)%%K+1
				lll=(l-1)%%K+1
				Uhat=Uhat+c(WWtilde[j,l])*kronecker(Xmatstar[i+((1:length(tseq))-1)*N,],diag(rep(1,K)))[j,]*Ymatstar[i+(ll-1)*N,lll]
			}
		}
	}
	beta_star_cor=matrix(data=HHhat%*%Uhat,nrow=K,ncol=2+p+Q,byrow=FALSE)
	beta_star_cor
}

###Fitting code assuming specific working correlation matrix, 
## includes measurement error adjustment
## either true or one-step working correlation Sigma_e, Sigma_e_t
myfit_cor<-function(data,K,Q,p,tseq,Sigma_err,Sigma_e,Sigma_e_t){
	N=nrow(data)/length(tseq)
	###initial estimation 
	Sigma_I=Sigma_err[1:Q,1:Q]
	Sigma_S=Sigma_err[Q+c(1:K),Q+c(1:K)]
	Sigma_IS=Sigma_err[1:Q,Q+c(1:K)]
	###Measurement error correction
	Xmat=as.matrix(cbind(1,data[,c(2:(2+p+Q))]))
	Ymat=as.matrix(data[,2+p+Q+c(1:K)])
	YY=t(Xmat)%*%Ymat
	XX=t(Xmat)%*%Xmat
	###estimate parameter
	WW=solve(Sigma_e_t)
	Hhat1=matrix(data=0,nrow=2+p+Q,ncol=2+p+Q)
	###expanded sigma_I
	Omega_I=matrix(data=0,nrow=2+p+Q,ncol=2+p+Q)
	Omega_I[2+p+(1:Q),2+p+(1:Q)]=Sigma_I
	Hhat2=N*tr(WW)*Omega_I
	for (j in 1:length(tseq)){
		for (k in 1:length(tseq)){
			Hhat1=Hhat1+c(WW[j,k])*t(Xmat[1:N+(j-1)*N,])%*%Xmat[1:N+(k-1)*N,]
		}
	}
	###need make sure positivity
	if (min(eigen(Hhat1-nrow(Xmat)/(nrow(Xmat)-1)*Hhat2)$values)>0){
		Hhat=Hhat1-Hhat2
	}else{
		myf<-function(lambda){
			abs(min(eigen(Hhat1-lambda*Hhat2)$values))
		}
		lambda=optimize(myf,c(0,nrow(Xmat)/(nrow(Xmat)-1)))$minimum
		Hhat=Hhat1-(lambda-6/(nrow(Xmat)-1))*Hhat2
	}	
	HHhat=kronecker(solve(Hhat),Sigma_e)
	
	###update U
	SSigma_e=kronecker(Sigma_e_t, Sigma_e)
	WWtilde=solve(SSigma_e)
	Omega_IS=matrix(data=0,nrow=2+p+Q,ncol=K)
	Omega_IS[2+p+(1:Q),]=Sigma_IS
	Uhat=matrix(data=0,nrow=(2+p+Q)*K,ncol=1)
	for (i in 1:N){
		for (j in 1:(length(tseq)*K)){
			for (l in 1:(length(tseq)*K)){
				jj=(j-1)%/%K+1
				ll=(l-1)%/%K+1
				jjj=(j-1)%%K+1
				lll=(l-1)%%K+1
				Uhat=Uhat+c(WWtilde[j,l])*kronecker(Xmat[i+((1:length(tseq))-1)*N,],diag(rep(1,K)))[j,]*Ymat[i+(ll-1)*N,lll]
				if (jj==ll){
					e_lll=rep(0,K)
					e_lll[lll]=1
					e_jjj=rep(0,K)
					e_jjj[jjj]=1
					Uhat=Uhat-c(WWtilde[j,l])*kronecker(Omega_IS%*%e_lll,e_jjj)
				}
			}
		}
	}
	beta_onestep=matrix(data=HHhat%*%Uhat,nrow=K,ncol=2+p+Q,byrow=FALSE)
	beta_onestep
}

###code to estimate parameter related to correlation structure: Sigma_e,Sigma_e_t
## myest_par, that Kernel estimated variance-covariance matrix SSigma_e_hat1 
## is not always positive definite when measurement error is large 
## and sometimes SSigma_e_hat are ill-conditioned

myest_par<-function(data,K,Q,p,tseq,beta_ind,Sigma_err,h,type){
	N=nrow(data)/length(tseq)
	###initial estimation 
	Sigma_I=Sigma_err[1:Q,1:Q]
	Sigma_S=Sigma_err[Q+c(1:K),Q+c(1:K)]
	Sigma_IS=Sigma_err[1:Q,Q+c(1:K)]

	###use one step update, estimate rho_e, Sigma_e
	###get empirical estimation for SSigma_e
	SSigma_e_hat1=matrix(data=0,nrow=K*length(tseq),ncol=K*length(tseq))
	b4=beta_ind[,2+p+c(1:Q)]	
	Xmat=as.matrix(cbind(1,data[,c(2:(2+p+Q))]))
	Ymat=as.matrix(data[,2+p+Q+c(1:K)])
	Ehat=Ymat-Xmat%*%t(beta_ind)
	for (jj in 1:length(tseq)){
		for (kk in 1:length(tseq)){
			upper=matrix(data=0,nrow=K,ncol=K)
			lower=0
			if (jj==kk){
				for (i in 1:N){
					for (j in 1:length(tseq)){
						ehat_ij=Ehat[i+(j-1)*N,]
						w_ij=myK(abs(tseq[j]-tseq[jj])/h)
						upper=upper+w_ij*(ehat_ij%o%ehat_ij)
						lower=lower+w_ij
					}
				}
				SSigma_e_hat1[(jj-1)*K+(1:K),(kk-1)*K+(1:K)]=upper/lower				
			}
			if (jj!=kk){
				for (i in 1:N){
					for (j in 1:length(tseq)){
						for (k in 1:length(tseq)){
							ehat_ij=Ehat[i+(j-1)*N,]
							ehat_ik=Ehat[i+(k-1)*N,]
							if (k==j){
								w_ijk=0
							}else{
								w_ijk=myK(mydist(c(tseq[j],tseq[k]),c(tseq[jj],tseq[kk]))/h)
							}
							upper=upper+w_ijk*(ehat_ij%o%ehat_ik)
							lower=lower+w_ijk
						}
					}
				}
				SSigma_e_hat1[(jj-1)*K+(1:K),(kk-1)*K+(1:K)]=upper/lower	
			}
		}
	}
	###make sure SSigma_e_hat1 is positive definite
	SSigma_e_hat1=nearPD(SSigma_e_hat1)$mat

	SSigma_e_hat2=kronecker(diag(rep(1,length(tseq))),Sigma_S-b4%*%(Sigma_IS)-t(b4%*%(Sigma_IS))+b4%*%Sigma_I%*%t(b4))
	###make positivity hold
	if (min(eigen(SSigma_e_hat1-nrow(Xmat)/(nrow(Xmat)-1)*SSigma_e_hat2)$values)>0){
		SSigma_e_hat=SSigma_e_hat1-SSigma_e_hat2
	}else{
		myf<-function(lambda){
			abs(min(eigen(SSigma_e_hat1-lambda*SSigma_e_hat2)$values))
		}
		lambda=optimize(myf,c(0,nrow(Xmat)/(nrow(Xmat)-1)))$minimum
		while(min(eigen(SSigma_e_hat1-lambda*SSigma_e_hat2)$values)<0){
			lambda=lambda/2
		}
		SSigma_e_hat=SSigma_e_hat1-(lambda-6/(nrow(Xmat)-1))*SSigma_e_hat2
	}	
        ###get initial value
	tmpmat=SSigma_e_hat[1:K,1:K]
	if (length(tseq)>1){
		for (tt in 2:length(tseq)){
			tmpmat=tmpmat+SSigma_e_hat[K*(tt-1)+(1:K),K*(tt-1)+(1:K)]	
		}
		tmpmat=tmpmat/length(tseq)
	}
	par_ini=c(0,0,mat2vec(tmpmat))

	###estimate parameter
	myobj<-function(par,type){
		rho_e=par[2]
		sigma_e=exp(par[1])
		Sigma_e=vec2mat(par[-c(1,2)])
		if (type=="AR"){
			Sigma_e_t=Sigma_AR(tseq,c(sigma_e,rho_e))
		}
		if (type=="CS"){
			Sigma_e_t=Sigma_CS(tseq,c(sigma_e,rho_e))
		}
		SSigma_e=kronecker(Sigma_e_t, Sigma_e)
		sum(log(eigen(SSigma_e)$value))+tr(solve(SSigma_e,SSigma_e_hat))
	}
	par_est=NA
	try({fit=optim(par_ini,myobj,type=type)
	par_est=fit$par})
	if (max(is.na(par_est))==1){
            ## estimate (sigma_e,rho_e), set Sigma_e as initial values
		myobj_rr<-function(par1,par2,type){
			rho_e=par1[2]
			sigma_e=exp(par1[1])
			Sigma_e=vec2mat(par2)
			if (type=="AR"){
				Sigma_e_t=Sigma_AR(tseq,c(sigma_e,rho_e))
			}
			if (type=="CS"){
				Sigma_e_t=Sigma_CS(tseq,c(sigma_e,rho_e))
			}
			SSigma_e=kronecker(Sigma_e_t, Sigma_e)
			sum(log(eigen(SSigma_e)$value))+tr(solve(SSigma_e,SSigma_e_hat))
		}
		try({fit=optim(par_ini[1:2],myobj_rr,type=type,par2=par_ini[-c(1:2)])
		par_est=c(fit$par,par_ini[-c(1:2)])})		
	}
	if (is.na(par_est)[1] || is.na(par_est)[2]){
		par_est=par_ini
	}
	###estimate covariance matrix
	rho_e=par_est[2]
	sigma_e=exp(par_est[1])
	Sigma_e=vec2mat(par_est[-c(1,2)])
	if (type=="AR"){Sigma_e_t=Sigma_AR(tseq,c(sigma_e,rho_e))}
	if (type=="CS"){Sigma_e_t=Sigma_CS(tseq,c(sigma_e,rho_e))}
	list(Sigma_e,Sigma_e_t)
}

###############################################################

myfit<-function(data,K,Q,p,tseq,Sigma_err,h=1,type,Sigma_e_t.true){	
	beta_ols=myfit_ols(data,K,Q,p)
	beta_star=myfit_oracle(data,K,Q,p)
	## first fit beta with no V 
      beta_ind=myfit_ind(data,K,Q,p,Sigma_err)
	## obtain V given beta_ind
      par=myest_par(data,K,Q,p,tseq,beta_ind,Sigma_err,h,type)
	Sigma_e=par[[1]]
	Sigma_e_t=par[[2]]
      ## after getting V, fit beta with V
      ## Oracle method using estimated working correlation
      beta_star_cor=myfit_oracle_cor(data,K,Q,p,tseq,Sigma_e,Sigma_e_t) 
      ## measurement correction method using estimated working correlation (one-iteration)
	beta_onestep=myfit_cor(data,K,Q,p,tseq,Sigma_err,Sigma_e,Sigma_e_t)
      par=myest_par(data,K,Q,p,tseq,beta_onestep,Sigma_err,h,type)
	Sigma_e=par[[1]]
	Sigma_e_t=par[[2]]
      ## after getting V, fit beta with V (second time)
	## measurement correction method using estimated working correlation (two iterations)
      beta_twostep=myfit_cor(data,K,Q,p,tseq,Sigma_err,Sigma_e,Sigma_e_t)
	## measurement correction method using true working correlation
      beta_truecor=myfit_cor(data,K,Q,p,tseq,Sigma_err,Sigma_e.true,Sigma_e_t.true)
	## dimension K rows, (2+p+Q)*betaNum columns, 
      return(c(beta_star,beta_star_cor,beta_ols,beta_ind,beta_truecor,beta_onestep,beta_twostep))
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
Nvec=c(200,400,800)
tseq=c(1:3)
seed=1111
sigma_a=0.2
###genearte a random positive definite matrix from Wishart distribution
#Sigma_a=sigma_a*rWishart(1,2*Q,diag(2*Q))[,,1]
Sigma_a=Sigma_CS(1:(2*Q),c(1,-0.1))
p=2
rho_X=0.1
Sigma_X=Sigma_CS(1:p,c(1,rho_X))

###Define correlation over time
# the correlation between different times for specific domain.
rho=0.3
## rho_eAR=c(0.3,0.5,0.8)
rho_eAR=c(0.5)
rho_eCS=c(0.2,0.3,0.5)
rho_err=0

####generate Sigma randomly from Wishart distribution with a scalar sigma
#sigma=5
#Sigma=diag(diag(sigma*rWishart(1,Q,diag(Q))[,,1]))
## 0.2=correlation between different imaging domains at a specific time.
Sigma=Sigma_AR(1:Q,c(3,0.2))

####generate Sigma_e randomly from Wishart distribution with a scalar sigma_e
#sigma_e=0.3
#Sigma_e=sigma_e*rWishart(1,K,diag(K))[,,1]
## K by K compound symmetric covariance matrix 
## with parameters $\sigma^2=2.5$ and $\rho_e=0.15$
Sigma_e=Sigma_CS(1:K,c(2.5,0.15))
Sigma_e.true=Sigma_e

####generate Sigma_err (include Sigma_I,Sigma_S,Sigma_IS) randomly from Wishart distribution with a scalar sigma_err
## SNR=10,4,2,1 sigma_err=c(0.13,0.34,0.66,1.33) (see sigma_err_Estimate.R)
SNR=c(10,4,2,1)
sigma_err=c(0.13,0.34,0.66,1.33)
Wis=rWishart(1,K+Q,diag(K+Q))[,,1]
#Sigma_err=diag(c(rep(1,K),rep(1.2,Q)))

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

# beta_ (oracle, oracle_cor, ols, ind, truecor, onestep, twostep)
betaNum=7

## within each b=1,...,B, do bootstrap M times
M=100
## M=2

###generate one fitting result
onefit<-function(seed,rho_e,Sigma_e_t.true,type){
	set.seed(seed)
      onefitresult=NULL
      for (sigma_errIndex in 1:length(sigma_err)){
      Sigma_err=sigma_err[sigma_errIndex]*Wis
      for (Nindex in 1:length(Nvec)){
	N=Nvec[Nindex]
      alldata=datagen(N,K,Q,p,tseq,beta0,beta1,beta2,b0,b1,b2,b3,Sigma_a,Sigma_X,Sigma,Sigma_e,Sigma_err,rho, rho_e,rho_err,gamma,seed,type)
	data0=alldata$data0
	sname=paste("S",as.character(1:K),sep="")
	sname1=paste("Sstar",as.character(1:K),sep="")
	iname=paste("I",as.character(1:Q),sep="")
	iname1=paste("Istar",as.character(1:Q),sep="")
	## cannot know gamma and Sigma_err in the true setting, so check the performance of
      ## estimator when plug in estimated gamma and Sigma from normal samples.
      gammahat=rep(NA,K)
	Sigma_errhat=matrix(data=NA,nrow=K+Q,ncol=K+Q)
      ### residual from normal N/2 group
	res=matrix(data=NA,nrow=nrow(data0),ncol=K+Q)
	for (k in 1:K){	
		yy=data0[,sname[k]]
		tt=data0$t
		id=data0$ID
		fit=lm(yy~tt+as.factor(id))
		gammahat[k]=fit$coef[2]
		res[,Q+k]=fit$res
	} 
	for (q in 1:Q){	
		yy=data0[,iname[q]]
		id=data0$ID
		res[,q]=lm(yy~as.factor(id))$res
	} 
      #### degree of freedom excluding fixed intercept (N intercepts)
      ## (N*n_i-N) is used to adjust for degrees of freedom used for estimating 
      ## the fixed intercept for the N individuals.
	Sigma_errhat=crossprod(res)/(nrow(data0)-length(unique(data0$ID)))
	data=alldata$data
	data[,sname]=data[,sname]-data$t%o%gammahat
	data[,sname1]=data[,sname1]-data$t%o%gammahat
	est=myfit(data,K,Q,p,tseq,Sigma_errhat,h=1,type,Sigma_e_t.true)

## for one sample of N subjects (N/2 normal), do M bootstrap
## the Bootstrap code within onefit function so that it report not only
## the estimates but also the estimated SE from Bootstrap (use 100 Bootstraps). 
## Within onefit function, after generate the data, you need add a loop to 
## resampale cases and controls to form new dataset. You don t need to compute 
## bias for Bootstrap, you just need to get SE.
## What you need to do is the following assuming Sample size N, Bootstrap time M, 
## Simulation time B
## Step 1: Simulation B samples each with sample size N cases and N/2 controls
## Step 2: For each simulated sample, get the estimate est[b]
## Step 3: For each simulated sample, get M bootstrap sample and 
## calculate the estimate est[b,m]. Then compute the sample variance over m 
## to obtain se[b].
## Step 4: For each simulated sample, compute whether the 95% confidence interval 
## will cover the truth or not cr[b]= as.numeric(abs(est[b]-true.est)/se[b])<1.96)
## Step 5: Compute bias as the average over b for (est[b]-est.true), 
## ESE as the average over b for se[b], SD is the variance over b for est[b], 
## CR is the average over b for cr[b]

      N=nrow(data)/length(tseq)
      Mest=array(data=NA,dim=c(K,betaNum*(2+p+Q),M))
      for (m in 1:M){
      Msample=sort(sample(1:N,replace=T))
      Msample0=sort(sample(1:(N/2),replace=T))
      Mdata0=Mdata=NULL
      for (i in 1:length(tseq)){
	Mdata0=rbind(Mdata0,data0[N/2*(i-1)+Msample0,])
        Mdata=rbind(Mdata,alldata$data[N*(i-1)+Msample,])
      }
      Mgammahat=rep(NA,K)
      MSigma_errhat=matrix(data=NA,nrow=K+Q,ncol=K+Q)
      ### residual from normal group N/2
	Mres=matrix(data=NA,nrow=nrow(Mdata0),ncol=K+Q)
	###redefine ID
	Mdata0$ID=rep(1:(N/2),length(tseq))
	Mdata$ID=rep(1:N,length(tseq))
	for (k in 1:K){	
		yy=Mdata0[,sname[k]]
		tt=Mdata0$t
		id=Mdata0$ID
		fit=lm(yy~tt+as.factor(id))
		Mgammahat[k]=fit$coef[2]
		Mres[,Q+k]=fit$res
	} 
	for (q in 1:Q){	
		yy=Mdata0[,iname[q]]
		id=Mdata0$ID
		Mres[,q]=lm(yy~as.factor(id))$res
	} 
      #### degree of freedom excluding fixed intercept (N intercepts)
	MSigma_errhat=crossprod(Mres)/(nrow(Mdata0)-nrow(Mdata0)/length(tseq))
	Mdata[,sname]=Mdata[,sname]-Mdata$t%o%Mgammahat
	Mdata[,sname1]=Mdata[,sname1]-Mdata$t%o%Mgammahat
	Mest[,,m]=myfit(Mdata,K,Q,p,tseq,MSigma_errhat,h=1,type,Sigma_e_t.true)
      }
      ## compute ESE
      est=matrix(data=est,nrow=K,ncol=betaNum*(2+p+Q),byrow=FALSE)
      Ese=apply(Mest,c(1,2),sd,na.rm=TRUE)
      ## CR (compute whether the 95% confidence interval will cover the truth or not)
      CR=abs(est-kronecker(matrix(rep(1,betaNum),nrow=1,ncol=betaNum),beta_true))/Ese<1.96      
      onefitresult=rbind(onefitresult,est,Ese,CR)   
     } # end Nindex
    } # end sigma_errIndex
     return(onefitresult)
    }

myresult_Sigma_e_t=function(rho_eVec,type){
for (rho_e in rho_eVec){
if (type=="AR") {Sigma_e_t.true=Sigma_AR(tseq,c(1,rho_e))}
if (type=="CS") {Sigma_e_t.true=Sigma_CS(tseq,c(1,rho_e))}	
#### repeat fitting B times (with B replications), and compute bias and SD
B=2
## B=2
result=array(data=NA,dim=c(K*3*length(Nvec)*length(SNR),betaNum*(2+p+Q),B))
## 8888+b or simseed+b
for (b in 1:B){try(result[,,b]<-onefit(simseed+b,rho_e,Sigma_e_t.true,type))}
bias=SD=Ese=CR=NULL
for (index in 1:(length(Nvec)*length(SNR))){
bias=rbind(bias,apply(result[(index-1)*3*K+c(1:K),,],c(1,2),mean,na.rm=TRUE)-kronecker(matrix(rep(1,betaNum),nrow=1,ncol=betaNum),beta_true))
SD=rbind(SD,apply(result[(index-1)*3*K+c(1:K),,],c(1,2),sd,na.rm=TRUE))
Ese=rbind(Ese,apply(result[(index-1)*3*K+K+c(1:K),,],c(1,2),mean,na.rm=TRUE))
CR=rbind(CR,apply(result[(index-1)*3*K+2*K+c(1:K),,],c(1,2),mean,na.rm=TRUE))
}

#### write out the results, (Bias, sd, Ese, CR)
result=data.frame(rbind(bias=round(bias,5),SD=round(SD,5),Ese=round(Ese,5),CR=round(CR,5)))
### rename column names
xname=paste("X",as.character(1:p),sep="")
iname=paste("I",as.character(1:Q),sep="")
inputname=c("intercept","t",xname,iname)
betaname=c("Oracle","Oraclecor","ols","Ind","truecor","onestep","twostep")
names(result)=paste(rep(betaname,each=length(inputname)),inputname,sep=" ")
## rename row names
Nname=paste("N=",Nvec,sep="")
sname=paste("S",as.character(1:K),sep="")
Nsname=paste(rep(Nname,each=length(sname)),sname,sep=" ")
SNRname=paste("SNR=",SNR,sep="")
SNRNsname=paste(rep(SNRname,each=length(Nsname)),Nsname,sep=" ")
outputname=c("Bias","SD","Ese","CR")
rownames(result)=paste(rep(outputname,each=length(SNRNsname)),SNRNsname,sep=" ")
write.csv(result,paste("result3",type,"rho_e=",rho_e,"datagenseed=",as.character(seed),"simseed=",as.character(simseed),".csv"))
} # end rho_e
} # end myresult_Sigma_e_t

myresult_Sigma_e_t(rho_eAR,type="AR")