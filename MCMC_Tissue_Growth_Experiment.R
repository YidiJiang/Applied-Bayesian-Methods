# Run Markov chain Monte Carlo algorithm with 20,000 iterations
# prior
alpha <- 1
beta <- 1
alpha0 <- 1
beta0 <- 1
mu00 <- 0
tau00 <- 1
# load the data
Temp1 <- c(1.13,1.20,1.00,0.91,1.05)
Temp2 <- c(1.75,1.45,1.55,1.64,1.60)
Temp3 <- c(2.30,2.15,2.25,2.40,2.49)
Temp4 <- c(3.18,3.10,3.28,3.35,3.12)
Temp1bar <- mean(Temp1)
Temp2bar <- mean(Temp2)
Temp3bar <- mean(Temp3)
Temp4bar <- mean(Temp4)
n1 <- length(Temp1)
n2 <- length(Temp2)
n3 <- length(Temp3)
n4 <- length(Temp4)
data <- cbind(Temp1,Temp2,Temp3,Temp4)
# set the starting values
mu <- 0
mu0 <- 0
tau <- 1
tau0 <- 1
alpha <- 1
alpha0 <- 1
beta <- 1
beta0 <- 1
muStart<-1
mu0Start<-1
mu1Start<-1
mu2Start<-1
mu3Start<-1
mu4Start<-1
tauStart<-0.0025
tau0Start<-0.0025
# The posterior distributions are as follows:
mu0p <- function(data, mu0, tau0, tau){
  n <- length(data)
  meandata <- mean(data)
  mu0p <- (tau0*mu0 + n*tau*meandata) / (tau0 + n*tau)
  return(mu0p)
}
mu00p <- function(mu, tau0, mu00, tau00){
  mu_mean <- mean(mu)
  mu_n <- length(mu)
  mu00p <- (mu_n*mu_mean*tau0 + mu00*tau00) / (mu_n*tau0 + tau00)
  return(mu00p)
}
tau0p <- function(data, tau0, tau){
  n <- length(data)
  tau0p <- tau0 + n*tau
  return(tau0p)
}
tau00p <- function(mu, tau0, tau00){
  mu_n <- length(mu)
  tau00p <- mu_n*tau0 + tau00
  return(tau00p)
}
alphap <- function(data, alpha) {
  data_n <- apply(data,c(2),function(x){length(x)})
  alphap <- alpha + sum(data_n)/2
  return(alphap)
}
alpha0p <- function(mu, alpha0){
  mu_n <- length(mu)
  alpha0p <- alpha0 + mu_n/2
  return(alpha0p)
}
betap <- function(data, beta, mu){
  datasum <- sum(sum((data[,1]-mu[1])^2),sum((data[,2]-mu[2])^2),sum((data[,3]-mu[3])^2),sum((data[,4]-mu[4])^2))
  beta_p <- beta + 0.5*datasum
  return(beta_p)
}
beta0p <- function(mu, beta0, mu0){
  mu_n <- length(mu)
  beta0p <- beta0 + 0.5*sum((mu-mu0)^2)
  return(beta0p)
}
# MCMC parameters
nbig<-20000;
set.seed(333);
# Initialize chain
data_simulation <- matrix(NA, nrow=nbig, ncol=21)
colnames(data_simulation) <- c("mu0","mu1","mu2","mu3","mu4",
                               "mu00p","mu01p","mu02p","mu03p","mu04p",
                               "tau","tau0","tau00p","tau01p","tau02p","tau03p","tau04p",
                               "alphap","betap","alpha0p","beta0p")
data_simulation[1,"mu0"] <- mu0
data_simulation[1,"mu1"] <- mu
data_simulation[1,"mu2"] <- mu
data_simulation[1,"mu3"] <- mu
data_simulation[1,"mu4"] <- mu
data_simulation[1,"mu00p"] <- mu00
data_simulation[1,"mu01p"] <- mu0
data_simulation[1,"mu02p"] <- mu0
data_simulation[1,"mu03p"] <- mu0
data_simulation[1,"mu04p"] <- mu0
data_simulation[1,"tau"] <- tau
data_simulation[1,"tau0"] <- tau00
data_simulation[1,"tau00p"] <- mu00
data_simulation[1,"tau01p"] <- tau0
data_simulation[1,"tau02p"] <- tau0
data_simulation[1,"tau03p"] <- tau0
data_simulation[1,"tau04p"] <- tau0
data_simulation[1,"alphap"] <- alpha
data_simulation[1,"betap"] <- beta
data_simulation[1,"alpha0p"] <- alpha0
data_simulation[1,"beta0p"] <- beta0
for(i in 2:nbig){
  # sample mu1,mu2,mu3,mu4
  data_simulation[i,"mu01p"] <- mu0p(data=data[,"Temp1"], mu0=data_simulation[i-1,"mu0"], tau0=data_simulation[i-1,"tau0"], tau=data_simulation[i-1,"tau"])
  data_simulation[i,"tau01p"] <- tau0p(data=data[,"Temp1"], tau0=data_simulation[i-1,"tau0"], tau=data_simulation[i-1,"tau"])
  data_simulation[i,"mu1"] <- rnorm(1, data_simulation[i,"mu01p"], 1/sqrt(data_simulation[i,"tau01p"]))
  data_simulation[i,"mu02p"] <- mu0p(data=data[,"Temp2"], mu0=data_simulation[i-1,"mu0"], tau0=data_simulation[i-1,"tau0"], tau=data_simulation[i-1,"tau"])
  data_simulation[i,"tau02p"] <- tau0p(data=data[,"Temp2"], tau0=data_simulation[i-1,"tau0"], tau=data_simulation[i-1,"tau"])
  data_simulation[i,"mu2"] <- rnorm(1, data_simulation[i,"mu02p"], 1/sqrt(data_simulation[i,"tau02p"]))
  data_simulation[i,"mu03p"] <- mu0p(data=data[,"Temp3"], mu0=data_simulation[i-1,"mu0"], tau0=data_simulation[i-1,"tau0"], tau=data_simulation[i-1,"tau"])
  data_simulation[i,"tau03p"] <- tau0p(data=data[,"Temp3"], tau0=data_simulation[i-1,"tau0"], tau=data_simulation[i-1,"tau"])
  data_simulation[i,"mu3"] <- rnorm(1, data_simulation[i,"mu03p"], 1/sqrt(data_simulation[i,"tau03p"]))
  data_simulation[i,"mu04p"] <- mu0p(data=data[,"Temp4"], mu0=data_simulation[i-1,"mu0"], tau0=data_simulation[i-1,"tau0"], tau=data_simulation[i-1,"tau"])
  data_simulation[i,"tau04p"] <- tau0p(data=data[,"Temp4"], tau0=data_simulation[i-1,"tau0"], tau=data_simulation[i-1,"tau"])
  data_simulation[i,"mu4"] <- rnorm(1,data_simulation[i,"mu04p"], 1/sqrt(data_simulation[i,"tau04p"]))
  # sample mu0
  data_simulation[i,"mu00p"] <- mu00p(mu=data_simulation[i,c("mu1","mu2","mu3","mu4")],
                                      tau0=data_simulation[i-1,"tau0"], mu00=mu00, tau00=tau00)
  data_simulation[i,"tau00p"] <- tau00p(mu=data_simulation[i,c("mu1","mu2","mu3","mu4")],
                                        tau0=data_simulation[i-1,"tau0"], tau00=tau00)
  data_simulation[i,"mu0"] <- rnorm(1,data_simulation[i,"mu00p"], 1/sqrt(data_simulation[i,"tau00p"]))
  # sample tau
  data_simulation[i,"alphap"] <- alphap(data=data, alpha=alpha)
  data_simulation[i,"betap"] <- betap(data=data, beta=beta, mu=data_simulation[i,c("mu1","mu2","mu3","mu4")])
  data_simulation[i,"tau"] <- rgamma(1, data_simulation[i,"alphap"], data_simulation[i,"betap"])
  # sample tau0
  data_simulation[i,"alpha0p"] <- alpha0p(mu=data_simulation[i,c("mu1","mu2","mu3","mu4")], alpha0=alpha0)
  data_simulation[i,"beta0p"] <- beta0p(mu=data_simulation[i,c("mu1","mu2","mu3","mu4")], beta0=beta0, mu0=data_simulation[i,"mu0"])
  data_simulation[i,"tau0"] <- rgamma(1, data_simulation[i,"alpha0p"], data_simulation[i,"beta0p"])
}

# Posterior distribution for the mean of each temperature effect for ??1-??4:
p1 <- plot(density(data_simulation[,"mu1"]),xlab='mu1',ylab='Density (posterior distribution)',main='the postrior density of mu1')
p2 <- plot(density(data_simulation[,"mu2"]),xlab='mu2',ylab='Density (posterior distribution)',main='the postrior density of mu2')
p3 <- plot(density(data_simulation[,"mu3"]),xlab='mu3',ylab='Density (posterior distribution)',main='the postrior density of mu3')
p4 <- plot(density(data_simulation[,"mu4"]),xlab='mu4',ylab='Density (posterior distribution)',main='the postrior density of mu4')
# Combine the plots
library(ggplot2)
data1 <- data.frame("Mean of each temperature effect" = c(data_simulation[,"mu1"],data_simulation[,"mu2"],data_simulation[,"mu3"],data_simulation[,"mu4"]),"mu"=c(rep("mu1",nbig),rep("mu2",nbig),rep("mu3",nbig),rep("mu4",nbig)))
ggplot(data1, aes(x = Mean.of.each.temperature.effect, color = mu)) +
  geom_density(alpha = 0.5, size = 1.2) + scale_color_grey()+
  theme_bw(base_size = 12) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_continuous() +
  ylab("Density (posterior distribution)") + xlab("Mean of each temperature effect")
# Summary statistics
summary(data_simulation[,"mu1"])
sd(data_simulation[,"mu1"])
quantile(data_simulation[,"mu1"],prob=c(.025,.975))
summary(data_simulation[,"mu2"])
sd(data_simulation[,"mu2"])
quantile(data_simulation[,"mu2"],prob=c(.025,.975))
summary(data_simulation[,"mu3"])
sd(data_simulation[,"mu3"])
quantile(data_simulation[,"mu3"],prob=c(.025,.975))
summary(data_simulation[,"mu4"])
sd(data_simulation[,"mu4"])
quantile(data_simulation[,"mu4"],prob=c(.025,.975))

# Posterior distribution of the difference in number of cells at the temperature of 40 versus 80
data2 <- data.frame("Mean of each temperature effect" = c(data_simulation[,"mu1"] - data_simulation[,"mu3"]),
                    "mu"=c(rep("mu1-mu3",nbig)))
p6 <- plot(density(data2[,"Mean.of.mu1.mu3"]),xlab="Mean of mu1-mu3",ylab="Density (posterior distribution)",main="mu1-mu3")
# Summary statistics
diff <- data_simulation[,"mu1"] - data_simulation[,"mu3"]
mean(diff)
sd(diff)
quantile(diff,probs=c(.025,.975))
1-pnorm(q=0,mean=-1.226428, sd=0.2353597) # posterior distribution




