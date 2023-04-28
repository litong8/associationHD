library(MASS)
library(Matrix)

### Batch=number of files from parallel result, B=number of replications in each file
seed=1111
Batch=500
B=2
K=4
p=2
Q=3
Nvec=c(200,400,800)
betaNum=7
SNR=c(10,4,2,1)

rho_eAR=c(0.3,0.5,0.8)
rho_eCS=c(0.2,0.3,0.5)

for (type in c("AR","CS"))
{
if (type=="AR") {rho_eVec=rho_eAR}
else {rho_eVec=rho_eCS}
for (rho_e in rho_eVec) {
result=array(data=NA,dim=c(K*4*length(Nvec)*length(SNR),betaNum*(2+p+Q),Batch))
### obtain result from each file
files = list.files(pattern=paste(type,"rho_e=",rho_e,"datagenseed=",seed, "simseed="))
for (i in 1:Batch) {result[,,i]=as.matrix(read.csv(file = files[i],row.names = 1))}

## compute final result
## For bias, ESE, CR, simply take average over the Batch files
bias=apply(result[1:(K*length(Nvec)*length(SNR)),,],c(1,2),mean,na.rm=TRUE)
## Sd: add up the mean within batch variance and the between batch variance, and finally take square root.
if (B>1) {
SumSqWithin=apply(result[c(1:(K*length(Nvec)*length(SNR)))+K*length(Nvec)*length(SNR),,]^2*(B-1),c(1,2),sum,na.rm=TRUE)
SumSqBetween=apply(sweep(result[1:(K*length(Nvec)*length(SNR)),,],c(1,2),bias,FUN="-")^2*B,c(1,2),sum,na.rm=TRUE)
SD=sqrt((SumSqWithin+SumSqBetween)/(B*Batch-1))
}
else
{ # B=1
SD=apply(result[1:(K*length(Nvec)*length(SNR)),,],c(1,2),sd,na.rm=TRUE)
}
Ese=apply(result[c(1:(K*length(Nvec)*length(SNR)))+2*K*length(Nvec)*length(SNR),,],c(1,2),mean,na.rm=TRUE)
 CR=apply(result[c(1:(K*length(Nvec)*length(SNR)))+3*K*length(Nvec)*length(SNR),,],c(1,2),mean,na.rm=TRUE)

#### write out the Final results, (Bias, sd, Ese, CR)
finalresult=data.frame(rbind(bias=round(bias,3),SD=round(SD,3),Ese=round(Ese,3),CR=round(CR,3)))
### rename column names
xname=paste("X",as.character(1:p),sep="")
iname=paste("I",as.character(1:Q),sep="")
inputname=c("intercept","t",xname,iname)
betaname=c("Oracle","Oraclecor","ols","Ind","truecor","onestep","twostep")
names(finalresult)=paste(rep(betaname,each=length(inputname)),inputname,sep=" ")
## rename row names
Nname=paste("N=",Nvec,sep="")
sname=paste("S",as.character(1:K),sep="")
Nsname=paste(rep(Nname,each=length(sname)),sname,sep=" ")
SNRname=paste("SNR=",SNR,sep="")
SNRNsname=paste(rep(SNRname,each=length(Nsname)),Nsname,sep=" ")
outputname=c("Bias","SD","Ese","CR")
rownames(finalresult)=paste(rep(outputname,each=length(SNRNsname)),SNRNsname,sep=" ")
write.csv(finalresult,paste("finalresult",type,"rho_e=",rho_e,"datagenseed=",as.character(seed),".csv"))
}
}

########################################################################################
## create table for publication purpose
for (type in c("AR","CS"))
{
if (type=="AR") {rho_eVec=rho_eAR}
else {rho_eVec=rho_eCS}
for (rho_e in rho_eVec) {
result=array(data=NA,dim=c(K*4*length(Nvec)*length(SNR),betaNum*(2+p+Q),Batch))
### obtain result from each file
files = list.files(pattern=paste(type,"rho_e=",rho_e,"datagenseed=",seed, "simseed="))
for (i in 1:Batch) {result[,,i]=as.matrix(read.csv(file = files[i],row.names = 1))}
## compute final result
bias=apply(result[1:(K*length(Nvec)*length(SNR)),,],c(1,2),mean,na.rm=TRUE)
## Sd: add up the mean within batch variance and the between batch variance, and finally take square root.
if (B>1) {
SumSqWithin=apply(result[c(1:(K*length(Nvec)*length(SNR)))+K*length(Nvec)*length(SNR),,]^2*(B-1),c(1,2),sum,na.rm=TRUE)
SumSqBetween=apply(sweep(result[1:(K*length(Nvec)*length(SNR)),,],c(1,2),bias,FUN="-")^2*B,c(1,2),sum,na.rm=TRUE)
SD=sqrt((SumSqWithin+SumSqBetween)/(B*Batch-1))
}
else
{ # B=1
SD=apply(result[1:(K*length(Nvec)*length(SNR)),,],c(1,2),sd,na.rm=TRUE)
}
Ese=apply(result[c(1:(K*length(Nvec)*length(SNR)))+2*K*length(Nvec)*length(SNR),,],c(1,2),mean,na.rm=TRUE)
CR=apply(result[c(1:(K*length(Nvec)*length(SNR)))+3*K*length(Nvec)*length(SNR),,],c(1,2),mean,na.rm=TRUE)

for(Sdomain in 1:K)
{
resulttable=NULL

for (Idomain in 1:Q)
{
## select row:Sdomain; column: Idomain first six methods
biasIS=bias[(0:(length(Nvec)*length(SNR)-1))*K+Sdomain,(0:(betaNum-2))*(2+p+Q)+2+p+Idomain]
    SDIS=SD[(0:(length(Nvec)*length(SNR)-1))*K+Sdomain,(0:(betaNum-2))*(2+p+Q)+2+p+Idomain]
  EseIS=Ese[(0:(length(Nvec)*length(SNR)-1))*K+Sdomain,(0:(betaNum-2))*(2+p+Q)+2+p+Idomain]
    CRIS=CR[(0:(length(Nvec)*length(SNR)-1))*K+Sdomain,(0:(betaNum-2))*(2+p+Q)+2+p+Idomain]

## select row: SNR
biasISSNR=SDISSNR=EseISSNR=CRISSNR=c()
for (SNRIndex in 1:length(SNR))
{
biasISSNR=c(biasISSNR,as.vector(biasIS[1:length(Nvec)+(SNRIndex-1)*length(Nvec),])) 
SDISSNR=c(SDISSNR,as.vector(SDIS[1:length(Nvec)+(SNRIndex-1)*length(Nvec),]))
EseISSNR=c(EseISSNR,as.vector(EseIS[1:length(Nvec)+(SNRIndex-1)*length(Nvec),]))
CRISSNR=c(CRISSNR,as.vector(CRIS[1:length(Nvec)+(SNRIndex-1)*length(Nvec),]))
} # end SNRIndex

resulttable=cbind(resulttable,biasISSNR,SDISSNR,EseISSNR,CRISSNR)
} # end Idomain

write.csv(round(resulttable,3),paste("resulttable",type,"rho_e=",rho_e,"Sdomain=",Sdomain,".csv"))
} # end Sdomain 

} # end rho_e
} # end type