# Fit models to the data with a Bayesian analysis use Openbugs 

# Fit model 1: Yi is the difference between two analysis (meta analysis with i=12), delta is a normal 
# distribution with mean= theta (true mean: difference between treatment group and control group), there is
# only one true mean for 12 studies.
library(R2OpenBUGS)
DiabetesDrugEffect <- read.table(file.choose(), sep=",", header=T)
attach(DiabetesDrugEffect)
data=list("StudyID", "StudyN", "diff", "Sediff")
init.fun=function(){list(
  delta=rnorm(12,-1.5,0.25),
  theta=rnorm(1,-1.5,0.25),
  sd0=runif(1,0,1)
)}
param=c("delta", "theta", "tau0", "sd0")
cat(" model{
for(i in 1:12){
diff[i]~dnorm(delta[i], tau[i])
delta[i]~dnorm(theta, tau0)
tau[i]<-1/Sediff[i]/Sediff[i]}
theta~dnorm(0,1)
sd0~dunif(0,1)
tau0<-1/sd0/sd0
}", file="DiabetesDrugEffect.txt")
diabetes=bugs(data, init.fun, param
              , model.file="DiabetesDrugEffect.txt",
              n.chains=5, n.iter=30000, n.burnin=10000,
              n.thin=5 ,debug=TRUE)
print(diabetes, digits.summary = 3)
output<-diabetes$sims.array
pdf("ACFPlotBetaS.pdf")
par(mfrow=c(3,2))
for(i in 1:12){
  acf(output[,1,paste0("delta[",i,"]")], main=paste0("delta[",i,"]"))
}
acf(output[,1,"theta"], main="theta")
acf(output[,1,"tau0"], main="tau0")
acf(output[,1,"sd0"], main="sd0") #sd0 is sigma0
dev.off()
# plots of autocorrelation
par(mfrow=c(3,2))
for(i in 1:12){
  acf(output[,1,paste0("delta[",i,"]")], main=paste0("delta[",i,"]"))
}
acf(output[,1,"theta"], main="theta")
acf(output[,1,"tau0"], main="tau0")
acf(output[,1,"sd0"], main="sd0")
dev.off()

# Fit model 2, every study has its own variance
DiabetesDrugEffect <- read.table(file.choose(), sep=",", header=T)
attach(DiabetesDrugEffect)
data=list("StudyID", "StudyN", "diff", "Sediff")
init.fun=function(){list(
  delta=rnorm(12,-1.5,0.25),
  theta=rnorm(1,-1.5,0.25),
  sd0=runif(1,0,1)
)}
parameters=c("sd", "theta", "sd0")
cat(" model{
for(i in 1:12){
diff[i]~dnorm(theta, tau[i])
sd[i]<-pow(sd_0,2)+pow(Sediff[i],2)
tau[i]<-1/sd[i]/sd[i]}
theta~dnorm(0,1)
sd0~dunif(0,1)
}", file="DiabetesDrugEffect.txt")
diabetes=bugs(data, init.fun, parameters
              , model.file="DiabetesDrugEffect.txt",
              n.chains=5, n.iter=30000, n.burnin=10000,
              n.thin=5 ,debug=TRUE)
print(diabetes, digits.summary = 3)
output<-diabetes$sims.array
pdf("ACFPlotBetaS.pdf")
par(mfrow=c(3,2))
for(i in 1:12){
  acf(output[,1,paste0("sd[",i,"]")], main=paste0("sd[",i,"]"))
}
acf(output[,1,"theta"], main="theta")
acf(output[,1,"sd0"], main="sd0")
dev.off()

# fit the model 2: In model 2, every study has its own variance
DiabetesDrugEffect <- read.table(file.choose(), sep=",", header=T)
attach(DiabetesDrugEffect)
data=list("StudyID", "StudyN", "diff", "Sediff")
init.fun=function(){list(
  delta=rnorm(12,-1.5,0.25),
  theta=rnorm(1,-1.5,0.25),
  sd0=runif(1,0,1)
)}
parameters=c("sd", "theta", "sd0")
cat(" model{
for(i in 1:12){
diff[i]~dnorm(theta, tau[i])
sd[i]<-pow(sd_0,2)+pow(Sediff[i],2)
tau[i]<-1/sd[i]/sd[i]}
theta~dnorm(0,1)
sd0~dunif(0,1)
}", file="DiabetesDrugEffect.txt")
diabetes=bugs(data, init.fun, parameters
              , model.file="DiabetesDrugEffect.txt",
              n.chains=5, n.iter=30000, n.burnin=10000,
              n.thin=5 ,debug=TRUE)
print(diabetes, digits.summary = 3)
output<-diabetes$sims.array
pdf("ACFPlotBetaS.pdf")
par(mfrow=c(3,2))
for(i in 1:12){
  acf(output[,1,paste0("sd[",i,"]")], main=paste0("sd[",i,"]"))
}
acf(output[,1,"theta"], main="theta")
acf(output[,1,"sd0"], main="sd0")
dev.off()
# Plots of autocorrelation
par(mfrow=c(3,2))
for(i in 1:12){
  acf(output[,1,paste0("sd[",i,"]")], main=paste0("sd[",i,"]"))
}
acf(output[,1,"theta"], main="theta")
acf(output[,1,"sd0"], main="sd0")
dev.off()