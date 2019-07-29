# The data contains number of lung cancers for different age groups and different histories of smoking
sum(SmokeAgeDeath$pyears)
722844/20  # = 36142.2, so the total person-years per categories is less than 36142.2
ln(36142.2)  #= 10.49522, so rate has to be bigger than 1/36142.2 and log(rate) > -10.49522, so log( base rate) should be between about -11 and 11.1/11/11 is about 0.008264463

# univariate posterior distributions of each level of smoking and age 
beta0 ~ dnorm(0, .008264463)
beta0.adj <- beta0 + mean(b[]) + mean(beta.s[])+ mean(beta.c[])
std ~ dunif(0, 9)
tau <- 1/std/std
std.s ~dunif(0, 5)
tau.s <- 1/std.s/std.s
std.a ~ dunif(0,5)
tau.a <- 1/std.a/std.a

# Fit a poisson model by OpenBug and MCMC is as follows:
library(R2OpenBUGS)
SmokeAgeDeath <- read.table(file.choose(), sep=",", header=T)
attach(SmokeAgeDeath)
data3=list("death", "pyears", "age", "smoke")
init.fun=function(){list(
  beta.a=rnorm(5) ,mu.a=rnorm(1),beta.s=rnorm(4),
  mu.s=rnorm(1),b=rnorm(20,0,1),beta0=rnorm(1),
  std.a=runif(1,1,2),std.s=runif(1,1,2),std=runif(1,1,2))}
params=c("beta0.adj", "beta.a.adj", "beta.s.adj", "std", "std.a", "std.s")
cat(" model{
for(i in 1:20)
{death[i]~dpois(lam[i])
log(lam[i]) <- log(pyears[i])+beta0+beta.a[age[i]] + beta.s[smoke[i]] + b[i]
b[i]~dnorm(0,tau)
b.adj[i]<-b[i]-mean(b[])}
beta0~dnorm(0,0.008264463)
beta0.adj<-beta0+mean(b[])+mean(beta.a[])+mean(beta.s[])
for(ia in 1:5){
beta.a[ia]~dnorm(mu.a,tau.a)
beta.a.adj[ia]<-beta.a[ia]-mean(beta.a[])}
mu.a~dnorm(0,0.0001)
for(is in 1:4){
beta.s[is]~dnorm(mu.s,tau.s)
beta.s.adj[is]<-beta.s[is]-mean(beta.s[])}
mu.s~dnorm(0,0.0001)
std ~ dunif(0, 9)
tau <- 1/std/std
std.a ~dunif(0, 5)
tau.a <- 1/std.a/std.a
std.s ~ dunif(0,5)
tau.s <- 1/std.s/std.s}", file="SmokeAgeDeath.txt")
smokeagedeath0=bugs(data3, init.fun, params
                    , model.file="SmokeAgeDeath.txt",
                    n.chains=5, n.iter=30000, n.burnin=10000,
                    n.thin=5 ,debug=TRUE) #for production
print(smokeagedeath0)

beta.s1 <- smokeagedeath0$summary["beta.s.adj[1]","mean"]#nonsmoker
beta.s4 <- smokeagedeath0$summary["beta.s.adj[4]","mean"]#smokes >20 cigarettes per day
mean <- exp(beta.s4 - beta.s1) # the mean of exponential of the beta difference between two smoking levels
summary(exp(output[,1,"beta.s.adj[4]"] - output[,1,"beta.s.adj[1]"]))
sd(exp(output[,1,"beta.s.adj[4]"] - output[,1,"beta.s.adj[1]"]))
quantile(exp(output[,1,"beta.s.adj[4]"] - output[,1,"beta.s.adj[1]"]),c(0.025,.975))