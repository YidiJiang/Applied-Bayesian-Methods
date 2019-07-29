
# Run between 10,000 to 30,000 iterations of the MCMC

library(R2OpenBUGS)
cat(
  "smoke obese snore male hypoten n
0 0 0 1 5 60
0 0 0 0 10 149
1 0 0 1 2 17
1 0 0 0 6 16
0 1 0 1 1 12
0 1 0 0 2 9
0 0 1 1 36 187
0 0 1 0 28 138
1 0 1 1 13 85
1 0 1 0 4 39
0 1 1 1 15 51
0 1 1 0 11 28
1 1 1 1 8 23
1 1 1 0 4 12
", file= "SmokeHyperData.txt")
SmokeHyper=read.table("SmokeHyperData.txt",header=TRUE,sep = "")
attach(SmokeHyper)
cat("
model{
for( i in 1:14){
hypoten[i] ~ dbin(mu[i], n[i])
logit(mu[i]) <- b0 + b.smok*smoke[i]+ b.ob*obese[i]+ b.sn*snore[i] +
b.male*male[i] + b.smsn*smoke[i]*snore[i] + b[i]
b[i] ~dnorm(0, tau.b)
}
b.smok ~ dnorm(0, .04) # so, sd =5. exp(5) ~ 148 which is huge
b.ob ~ dnorm(0, .04)
b.sn ~ dnorm(0, .04)
b.male ~ dnorm(0, .04)
b0 ~ dnorm(0, .04)
b.smsn ~dnorm(0, .04)
sd.b ~ dunif(0, 5)
tau.b <- 1/sd.b/sd.b
}
", file="SmokeHyperMod3.txt")
bugM3.dat=list("hypoten", "n", "smoke", "obese", "snore", "male") # variables required in the model
initM3.fun=function(){ list( b=runif(14,-.8,-.2),
                             b0=runif(1,-.8,-.2),
                             b.smok=runif(1,-.8,-.2),b.ob=runif(1,-.8,-.2), b.sn=runif(1,-.8,-.2),
                             b.male=runif(1,-.8,-.2), b.smsn=runif(1, -8,-.2), sd.b=runif(1,.2,.8)
) }
paramsM3=c("b.smok", "b.ob", "sd.b") # variables you want to monitor
SmokeHypeBaseM3=bugs(bugM3.dat, initM3.fun, paramsM3, model.file="SmokeHyperMod3.txt",
                     n.chains=3, n.iter=10000, n.burnin=1,
                     n.thin=1 , debug=TRUE
)
print(SmokeHypeBaseM3,dig=3)
if(T){
  SArray= SmokeHypeBaseM3$sims.array
  vname=attr(SArray,"dimnames")[3][[1]]
  chainL=attr(SArray,"dim")[1][[1]]
  for(i in 1:length(vname)){
    nn=vname[i]
    plot(density(SArray[,,nn]), main=nn)
    acf( SArray[,1,nn], main=paste(nn,"chain 1")) # auto correlation function plot for the 1st chain
    acf( SArray[,2,nn], main=paste(nn,"chain 2")) # auto correlation function plot for the 2nd chain
    acf( SArray[,3,nn], main=paste(nn,"chain 3")) # auto correlation function plot for the 3rd chain
    matplot(1:chainL,SArray[,,nn], main=nn,xlab="index",type="l")
  }
}

# Input the MCMC values into R. Plot trace plot. Remove early values of the chain (throwing away a part that is "burned-in")

burned_in <- 1000
if(T){
  SArray= SmokeHypeBaseM3$sims.array
  vname=attr(SArray,"dimnames")[3][[1]]
  chainL=attr(SArray,"dim")[1][[1]]
  for(i in 1:length(vname)){
    nn=vname[i]
    plot(density(SArray[burned_in:chainL,1,nn]), col=1, main= paste(nn,"after removing burned-in"))
    lines(density(SArray[burned_in:chainL,2,nn]), col=2)
    lines(density(SArray[burned_in:chainL,3,nn]), col=3)
    acf( SArray[burned_in:chainL,1,nn], main=paste(nn,"chain 1"))
    acf( SArray[burned_in:chainL,2,nn], main=paste(nn,"chain 2"))
    acf( SArray[burned_in:chainL,3,nn], main=paste(nn,"chain 3"))
    matplot(burned_in:chainL,SArray[burn_in:chainL,,nn], main= paste(nn,"after removing burned-in") ,xlab="index",type="l")
  }
}

# Estimate the posterior mean of three parameters for each chain, give the Monte Carlo accuracy of the estimation

library(coda)
#Obtain summaries of MCMC chain via coda package
make.mcmc.list=function(x){
  aa=x$sims.array
  zz=list(list())
  for(i in 1:(dim(aa)[2]) ){
    tmp=mcmc(aa[,i,])
    zz=c(zz,list(tmp)) }
  res=mcmc.list(zz[-1])
  res}
# Calculates the batch means estimate of the standard error of the mean when one inputs a vector of MCMC values:
CalcBatchMeans=function(x,Batn=50){
  BigN=length(x)
  BatInc=ceiling( (1:BigN)/(BigN/Batn) )
  BM=tapply(x,BatInc,mean)
  list(MCE=(sd(BM)/sqrt(length(BM))), BM=BM)}
# Calculate the MC standard error for a vector of sampled MCMC values:
CalcAcSe=function(x,lag.max=50){
  autoc=(acf(x,lag.max=lag.max,plot=FALSE))$acf
  sd(x)/sqrt(length(x))*sqrt(-1+2*sum(autoc))}
mcmc_list <- make.mcmc.list(SmokeHypeBaseM3)
for(i in 1:3){
  # For each of the three chain, calculate posterior mean
  mean <- apply(mcmc_list[[i]],2,mean)
  # Use batch mean method to calculate chain accuracy
  batchmean <- batchSE(mcmc_list[[i]])
  # Use Autocorrelation function to calculate autocorrelations
  auto_cor <- apply(mcmc_list[[i]],2,CalcAcSe)
  print(rbind(mean,batchmean,auto_cor))
}

# Use coda package, use the Geweke and Brooks-Gelman-Rubin diagnostic procedures to assess how well the MCMC algorithm has converged

# Geweke diagnose procedure
geweke.diag(mcmc_list)
# Brooks-Gelman-Rubin diagnose procedure;
gelman.diag(mcmc_list)