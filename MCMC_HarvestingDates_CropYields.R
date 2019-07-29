# The first model predicts the population as a linear function of the date 
# The second predicts with a quadratic function

# Compare these two models by the "deviance" measures and the DIC
harvest <- data.frame(x=c(16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46),
                      y=c(2508,2518,3304,3423,3057,3190,3500,3883,3823,3646,3708,
                          3333,3517,3241,3103,2776))
data <- list("x","y")
# The deviance and DIC for model 1
init.fun = function(){
  list(
    b=c(rnorm(1,0,.000001),rnorm(1,0,.01)),
    tau=rgamma(1,.0001,.0001)
  )
}
params<-c("b","tau","dev","dev.rep","dev.pval");
cat("
model{
for(i in 1:16){
y[i]~dnorm(mu[i],tau)
mu[i]<- b[1] + b[2]*(x[i]-31)
# model checking steps are here........getting the residuals for the observed values...
# note: I am deviating from the bugs manual... not getting the moments.
res[i]<-(y[i]-mu[i]) # estimate of the residuals for this model
stdres[i]<-res[i]*sqrt(tau) # for the standardized residuals
# getting a replicated sample..... This is a sample of the predictive distribution
y.rep[i]~dnorm(mu[i],tau)
# likelihood for each observed and replicated data....
# note: need to know the density function of the probability model
loglike[i]<- (0.5)*log(tau/6.283) + (-0.5)*tau*pow((y[i]-mu[i]),2)
loglike.rep[i]<- (0.5)*log(tau/6.283) + (-0.5)*tau*pow((y.rep[i]-mu[i]),2)
}
b[1]~dnorm(0,.000001)
b[2]~dnorm(0,.000001)
tau~dgamma(.0001,.0001)
dev<- -2*sum(loglike[])
dev.rep <- -2*sum(loglike.rep[])
dev.pval<-step(dev-dev.rep)
}
", file="model_1.txt")
model_1 <- bugs(data, init.fun,params ,model.file="model_1.txt",
                n.chains=3, n.iter=30000, n.burnin=10000, n.thin=10,
                debug=TRUE )
x1 <- model_1$sims.list$dev.rep
temp <- c(model_1$mean$dev,quantile(x1,probs=c(0.025,.975)),mean(x1),sd(x1))
names(temp) = c("Deviance","2.5%","97.5%","mean","SD");temp

#The code of Deviance and DIC for model 2 
init.fun <- function(){
  list(
    b=c(rnorm(1,0,.000001),rnorm(1,0,.000001),rnorm(1,0,.01)),
    tau=rgamma(1,.0001,.0001)
  )
}
parameters<-c("b","tau","dev","dev.rep","dev.pval");
cat("
model{
for(i in 1:16){
y[i]~dnorm(mu[i],tau)
mu[i]<- b[1] + b[2]*(x[i]-31)+ b[3]*pow((x[i]-31),2)
# getting a replicated sample..... This is a sample of the predictive distribution
y.rep[i]~dnorm(mu[i],tau)
# likelihood for each observed and replicated data....
# note: need to know the density function of the probability model
loglike[i]<- (0.5)*log(tau/6.283) + (-0.5)*tau*pow((y[i]-mu[i]),2)
loglike.rep[i]<- (0.5)*log(tau/6.283) + (-0.5)*tau*pow((y.rep[i]-mu[i]),2)
}
b[1]~dnorm(0,.000001)
b[2]~dnorm(0,.000001)
b[3]~dnorm(0,.01)
tau~dgamma(.0001,.0001)
dev<- -2*sum(loglike[])
dev.rep <- -2*sum(loglike.rep[])
dev.pval<-step(dev-dev.rep)
}
", file="model_2.txt")
model_2 <- bugs(data, init.fun,parameters ,model.file="model_2.txt",
                n.chains=3, n.iter=30000, n.burnin=10000, n.thin=10,
                debug=TRUE )
x2 <- model_2$sims.list$dev.rep
temp <- c(model_2$mean$dev,quantile(x1,probs=c(0.025,.975)),mean(x2),sd(x2))
names(temp) = c("Deviance","2.5%","97.5%","mean","SD");temp

# Compare two models by calculating Bayes Factor. 
# Run an MCMC algorithm which switches between two models to calculate the Bayes factor
# Use means and standard deviances to center the variable x, (x-31)^2 and y
mean(x)
sd(x)
mean((x-31)^2)
sd((x-31)^2)
mean(y)
sd(y)
init.fun <- function(){
  list(
    b=c(rnorm(1,0,.000001),rnorm(1,0,.01)),
    tau=rgamma(1,.0001,.0001),
    del=rbinom(1,1,0.5)
  )
}
params<-c("b","tau","del","mod");
cat("
model{
for(i in 1:16){
sx1[i] <- (x[i]-31) / 9.521905
sx2[i] <- (pow(x[i] - 31,2) - 85) / 78.05639
sy[i] <- (y[i]-3283.125) / 418.3434
sy[i]~dnorm(mu[i],tau)
mu[i]<- b[1]*sx1[i] + del*b[2]*sx2[i]
}
b[1]~dnorm(0,.000001)
b[2]~dnorm(0,.01)
tau~dgamma(.0001,.0001)
del~dbern(0.5)
for(i in 1:2){
mod[i] <- equals((i-1),del)
}
}
", file="bayes_model.txt")
bayes_model <- bugs(data,init.fun,params ,model.file="bayes_model.txt",
                    n.chains=3, n.iter=30000, n.burnin=10000, n.thin=10,
                    debug=TRUE )
bayes_model$summary

# For model 1, look at the residuals for the model, calculate the residuals, the standardized residuals, 
# and the chance of getting a more extreme observation

init.fun <- function(){
  list(
    b=c(rnorm(1,0,.000001),rnorm(1,0,.01)),
    tau=rgamma(1,.0001,.0001)
  )
}
params<-c("b","tau", "mu","res", "stdres", "res.rep", "stdres.rep", "p.smaller",
          "chidev1.pval", "chidev2.pval", "chidev1.obs", "chidev2.obs","chidev1.rep", "chidev2.rep");
cat("
model{
for(i in 1:16){
y[i]~dnorm(mu[i],tau)
mu[i]<- b[1] + b[2]*(x[i]-31)
# model checking steps are here.........
# getting the residuals for the observed values...
# note: I am deviating from the bugs manual... not getting the moments.
res[i]<-(y[i]-mu[i]) # estimate of the residuals for this model
stdres[i]<-res[i]*sqrt(tau) # for the standardized residuals
dev1.obs[i]<-pow(res[i],2)
dev2.obs[i]<-pow(stdres[i],2)
# getting a replicated sample..... This is a sample of the predictive distribution
y.rep[i]~dnorm(mu[i],tau)
p.smaller[i] <-step(y[i]-y.rep[i]) # check to see the probability of getting a more extreme value
# residual and moments of replicated data....this gives the predicted distribution for these values.
res.rep[i]<- y.rep[i] - mu[i]
stdres.rep[i]<- res.rep[i]*sqrt(tau)
dev1.rep[i]<-pow(res.rep[i],2)
dev2.rep[i]<-pow(stdres.rep[i],2)
}
b[1]~dnorm(0,.000001)
b[2]~dnorm(0,.000001)
tau~dgamma(.0001,.0001)
# summing the diagnostic values
chidev1.obs <- sum(dev1.obs[])
chidev2.obs <- sum(dev2.obs[])
chidev1.rep <- sum( dev1.rep[] )
chidev2.rep <- sum( dev2.rep[] )
chidev1.pval<-step(chidev1.obs-chidev1.rep)
chidev2.pval<-step(chidev2.obs-chidev2.rep)
}
", file="model_1.txt")
model1_res <- bugs(data, init.fun,params ,model.file="model_1.txt",
                   n.chains=3, n.iter=50000, n.burnin=25000, n.thin=10,
                   debug=TRUE )
# Get the residuals and the calibrations:
res <- cbind(model1_res $mean$res,t(apply(model1_res $sims.list$res.rep,2,function(x){c(quantile(x,probs=c(0.025,.975)),mean(x),sd(x))})))
colnames(res)=c("res","2.5%","97.5%","mean","SD");res
# Get the standardized residual and predictive distribution of the standardized residual
stdres <- cbind(model1_res $mean$stdres,t(apply(model1_res $sims.list$stdres.rep,2,function(x){c(quantile(x,probs=c(0.025,.975)),mean(x),sd(x))})))
colnames(stdres)=c("stdres","2.5%","97.5%","mean","SD");stdres
# Draw the Q-Q plot
qqplot(p.smaller[,1],runif(160,0,1), xlab="sample data (more extreme observations) Quantiles",
       ylab="Uniform distribution Quantiles")