
# The first project
# Project Background: Predict the results of an election. There are three select-persons, each select-person 
# is elected in one of three town districts and both parties have candidates running in each district. There
# are 5001 citizens voting in each district. 


# Prior distributions
# Info: In the past, each party would get between 40% to 60% of the votes.
# Informative prior
library(LearnBayes)
quantile1=list(p=0.025, x=0.4)
quantile2=list(p=0.975, x=0.6)
beta.select(quantile1,quantile2)
qbeta(c(0.025,0.975), 47.36,47.36)


# Posterior distribution for the voting percentage for each district
# Info: Random sample of the citizen of each district is asked whom they plan to vote for. 
# District 1: 53 purple, 45 brown; District 2: 72 purple, 78 brown; District 3: 18 purple, 22 brown.

# Non-informative prior
non_info <- cbind("District"=c(1,2,3),
                  "VoteForPurple"=c(53,72,18),
                  "VoteAll"=c(53+45,72+78,18+22),
                  "PriorAlpha"=c(1,1,1),
                  "PriorBeta"=c(1,1,1))
non_info_dist <- cbind(non_info,"PosteriorAlpha"=non_info[,"PriorAlpha"]+non_info[,"VoteForPurple"],
                       "PosteriorBeta"=non_info[,"PriorBeta"]+non_info[,"VoteAll"]-non_info[,"VoteForPurple"])

# Calculate the statistics for Beta distributions for 3 districts
non_info_stats <- cbind("District"=non_info_dist[,"District"],
                        "Posterior_??"=non_info_dist[,"PosteriorAlpha"],
                        "Posterior_??"=non_info_dist[,"PosteriorBeta"],
                        "PosteriorMean"=non_info_dist[,"PosteriorAlpha"]/(non_info_dist[,"PosteriorAlpha"]+non_info_dist[,"PosteriorBeta"]),
                        "95%CR_lowerbound"=qbeta(c(.025), non_info_dist[,"PosteriorAlpha"], non_info_dist[,"PosteriorBeta"]),
                        "95%CR_Upperbound"=qbeta(c(.975), non_info_dist[,"PosteriorAlpha"], non_info_dist[,"PosteriorBeta"]))

# Informative prior
info <- cbind("District"=c(1,2,3),
              "VoteForPurple"=c(53,72,18),
              "VoteAll"=c(53+45,72+78,18+22),
              "PriorAlpha"=c(47.36,47.36,47.36),
              "PriorBeta"=c(47.36,47.36,47.36))
info_dist <- cbind(info, "PosteriorAlpha"=info[,"PriorAlpha"]+info[,"VoteForPurple"],
                     + "PosteriorBeta"=info[,"PriorBeta"]+info[,"VoteAll"]-info[,"VoteForPurple"])
View(info_dist)
# Calculate the statistics for Beta distributions for 3 districts
info_stats <- cbind("District"=info_dist[,"District"],
                    "Posterior_??"=info_dist[,"PosteriorAlpha"],
                    "Posterior_??"=info_dist[,"PosteriorBeta"],
                    "PosteriorMean"=info_dist[,"PosteriorAlpha"]/(info_dist[,"PosteriorAlpha"]+info_dist[,"PosteriorBeta"]),
                    "95%CR_lowerbound"=qbeta(c(.025), info_dist[,"PosteriorAlpha"], info_dist[,"PosteriorBeta"]),
                    "95%CR_Upperbound"=qbeta(c(.975), info_dist[,"PosteriorAlpha"], info_dist[,"PosteriorBeta"]))


# Simulate 10000 elections. For each simulated election, generate a probability (theta i) of a citizen 
# voting for the purple party in district i, and the number of citizens voting for purple by district.

# Non-informative prior
n1 <-rbeta(10000,54,46)
n2 <-rbeta(10000,73,79)
n3 <-rbeta(10000,19,23)
n <- cbind(n1,n2,n3)
non_sim <- as.data.frame(matrix(NA,nrow=10000,ncol=3))
for (i in 1:10000){
  for (j in 1:3) {
    non_sim[i,j]=rbinom(1,5001,n[i,j])
  }
  }
non_sim_2 <- as.data.frame(matrix(NA,nrow=10000,ncol=3))

for (i in 1:10000){
  for (j in 1:3) {
    if (non_sim[i,j]>2500){
      non_sim_2[i,j] = 1 # purple party win
      } else {
        non_sim_2[i,j] = 0 # purple party lose
      }
  }
  }
# Calculate probabilities of a citizen voting for purple party in each district
win_1<- sum(non_sim_2[,1])
win_1/10000
win_2<- sum(non_sim_2[,2])
win_2/10000
win_3<- sum(non_sim_2[,3])
win_3/10000
non_sim_2$total <- rep(NA, times=10000)

for (i in 1:10000){
  non_sim_2[i,"total"]=sum(non_sim_2[i,1:3])
  if (non_sim_2[i,"total"]>=2){
    non_sim_2[i,"total"] =1
    } else {
      non_sim_2[i,"total"] =0
    }
  }
sum(non_sim_2$total)
sum(non_sim_2$total)/10000

# Informative prior
nn1 <-rbeta(10000,100.36,92.36)
nn2 <-rbeta(10000,119.36,125.36)
nn3 <-rbeta(10000,65.36,69.36)
nn <- cbind(nn1,nn2,nn3)
info_sim <- as.data.frame(matrix(NA,nrow=10000,ncol=3))
info_sim <- as.data.frame(matrix(NA,nrow=10000,ncol=3))
for (i in 1:10000){
  for (j in 1:3) {
    info_sim[i,j]=rbinom(1,5001,nn[i,j])
  }
  }
info_sim_2 <- as.data.frame(matrix(NA,nrow=10000,ncol=3))
for (i in 1:10000){
  for (j in 1:3) {
    if (info_sim[i,j]>2500){
      info_sim_2[i,j] = 1
      } else {
        info_sim_2[i,j] = 0
      }
  }
}
win_11 <- sum(info_sim_2[,1])
win_11/10000
win_22<- sum(info_sim_2[,2])
win_22/10000
win_33<- sum(info_sim_2[,3])
win_33/10000
info_sim_2$total <- rep(NA, times=10000)
for (i in 1:10000){
  info_sim_2[i,"total"]=sum(info_sim_2[i,1:3])
  if (info_sim_2[i,"total"]>=2){
    info_sim_2[i,"total"] =1
    } else {
      info_sim_2[i,"total"] =0
    }
  }
sum(info_sim_2$total)
sum(info_sim_2$total)/10000




# The second project
# Basic Bayesian methods: Calculate major statistics for each of three priors and for two datasets (6 combinations in total) 
dataset <- cbind(c(53,49,63,72,55,65),c(28,27,36,42,25,35))
colnames(dataset) <- c("dataset1","dataset2")
dataset <- rbind(dataset,c(mean(dataset[,1]),mean(dataset[,2])))
rownames(dataset)[7] <- c("mean")
prior1 <- rbind(c((4/9*66+6*1/36*dataset["mean","dataset1"])/(4/9+6*1/36),
                  (4/9*66+6*1/36*dataset["mean","dataset2"])/(4/9+6*1/36)),
                c((4/9+6*1/36),(4/9+6*1/36)))
prior1 <- rbind(prior1, sqrt(1/prior1[2,]),
                cbind(qnorm(c(0.025,0.975), mean=prior1[1,1], sd =sqrt(1/prior1[2,1])),
                      qnorm(c(0.025,0.975), mean=prior1[1,2], sd =sqrt(1/prior1[2,2]))))
rownames(prior1) <- c("Post_Mean","Post_tau","Post_sd","PostLowerCI","PostUpperCI")

Prior2 <- rbind(c((4*66+6*dataset["mean","dataset1"])/(4+6),(4*66+6*dataset["mean","dataset2"])/(4+6)),
                c(4+6,4+6),c(1+6/2,1+6/2))
rownames(Prior2) <- c("Post_Mean","Post_Theta","Post_Alpha")
dataset <- rbind(dataset,c(sum((dataset[,1]-mean(dataset[,1]))^2),sum((dataset[,2]-mean(dataset[,2]))^2)))
rownames(dataset)[8] <- c("SumSquare")
Prior2 <- rbind(Prior2, c((36+0.5*dataset["SumSquare","dataset1"]+(4*6*(dataset["mean","dataset1"]-66)^2)/(2*(4+6))),
                            + (36+0.5*dataset["SumSquare","dataset2"]+(4*6*(dataset["mean","dataset2"]-66)^2)/(2*(4+6)))))
rownames(Prior2)[4] <- "Post_Beta"
Prior2 <- rbind(Prior2,c(Prior2["Post_Alpha",1]*Prior2["Post_Theta",1]/Prior2["Post_Beta",1],
                           Prior2["Post_Alpha",2]*Prior2["Post_Theta",2]/Prior2["Post_Beta",2]))
rownames(Prior2)[5] <- "Post_Tau"
Prior2 <- rbind(Prior2,c(sqrt((1/Prior2["Post_Tau",])*(2*Prior2["Post_Alpha",])/
                                  (2*Prior2["Post_Alpha",]-2))))
rownames(Prior2)[6] <- "Post_SD"
Prior2 <- rbind(Prior2,cbind(qt(c(0.025,0.975), df=2*Prior2["Post_Alpha",1]),qt(c(0.025,0.975), df=2*Prior2["Post_Alpha",2])))
Prior2[7,] <- Prior2[7,]*sqrt(Prior2["Post_Beta",]/Prior2["Post_Theta",]/Prior2["Post_Alpha",])+ Prior2["Post_Mean",]
Prior2[8,] <- Prior2[8,]*sqrt(Prior2["Post_Beta",]/Prior2["Post_Theta",]/Prior2["Post_Alpha",])+ Prior2["Post_Mean",]
rownames(Prior2)[7:8] <- c("95%CIlower","95%CIupper")
colnames(Prior2) <- c("dataset1","dataset2")

Prior3 <- rbind(c((0.1*66+6*dataset["mean","dataset1"])/(0.1+6),(0.1*66+6*dataset["mean","dataset2"])/(0.1+6)),
                c(0.1+6,0.1+6),c(0.001+6/2,0.001+6/2))
rownames(Prior3) <- c("Post_Mean","Post_Theta","Post_Alpha")
Prior3 <- rbind(Prior3,c((0.001+0.5*dataset["SumSquare","dataset1"]+(0.1*6*(dataset["mean","dataset1"]-66)^2)/(2*0.1+2*6)),
                         (0.001+0.5*dataset["SumSquare","dataset2"]+(0.1*6*(dataset["mean","dataset2"]-66)^2)/(2*0.1+2*6))))
rownames(Prior3)[4] <- "Post_Beta"
Prior3 <- rbind(Prior3,c(Prior3["Post_Alpha",1]*Prior3["Post_Theta",1]/Prior3["Post_Beta",1],
                           Prior3["Post_Alpha",2]*Prior3["Post_Theta",2]/Prior3["Post_Beta",2]))
rownames(Prior3)[5] <- "Post_Tau"
Prior3 <- rbind(Prior3,c(sqrt((1/Prior3["Post_Tau",])*(2*Prior3["Post_Alpha",])/(2*Prior3["Post_Alpha",]-2))))
rownames(Prior3)[6] <- "Post_SD"
Prior3 <- rbind(Prior3,cbind(qt(c(0.025,0.975), df=2*Prior3["Post_Alpha",1]),qt(c(0.025,0.975), df=2*Prior3["Post_Alpha",2])))
Prior3[7,] <- Prior3[7,]*sqrt(Prior3["Post_Beta",]/Prior3["Post_Theta",]/Prior3["Post_Alpha",])+Prior3["Post_Mean",]
Prior3[8,] <- Prior3[8,]*sqrt(Prior3["Post_Beta",]/Prior3["Post_Theta",]/Prior3["Post_Alpha",]) +Prior3["Post_Mean",]
rownames(Prior3)[7:8] <- c("95%CIlower","95%CIupper")
colnames(Prior3) <- c("dataset1","dataset2")

# Give the density plots for six different posterior densities (all of the 6 combinations on the same plot)
# Can be btained by simulating 10000 samples using posterior distributions
# Color: Red-Prior1, Green-Prior2, Blue-Prior3; Line: solid- Data1, dash- Data2

n <- 10000
data <- cbind(rnorm(n,Prior1["Post_Mean",1],Prior1["post_sd",1]),
              rnorm(n,Prior1["Post_Mean",2],Prior1["post_sd",2]),
              rt(n,df=2*Prior2["Post_Alpha",1]),
              rt(n,df=2*Prior2["Post_Alpha",2]),
              rt(n,df=2*Prior3["Post_Alpha",1]),
              rt(n,df=2*Prior3["Post_Alpha",2]))
data[,3] <- data[,3]*sqrt(Prior2["Post_Beta",1]/Prior2["Post_Theta",1]/Prior2["Post_Alpha",1])+Prior2["Post_Mean",1]
data[,4] <- data[,4]*sqrt(Prior2["Post_Beta",2]/Prior2["Post_Theta",2]/Prior2["Post_Alpha",2])+Prior2["Post_Mean",2]
data[,5] <- data[,5]*sqrt(Prior3["Post_Beta",1]/Prior3["Post_Theta",1]/Prior3["Post_Alpha",1])+Prior3["Post_Mean",1]
data[,6] <- data[,6]*sqrt(Prior3["Post_Beta",2]/Prior3["Post_Theta",2]/Prior3["Post_Alpha",2])+Prior3["Post_Mean",2]
plot(density(data[,1]), xlim=c(0,100),col = "red", xlab = "mean", ylab = "density", main="", lwd = 2, lty=1)
lines(density(data[,2]), col = "red", lwd = 2, lty=2)
lines(density(data[,3]), xlim=c(0,100), col = "green", lwd = 2, lty=1)
lines(density(data[,4]), xlim=c(0,100), col = " green", lwd = 2, lty=2)
lines(density(data[,5]), xlim=c(0,100), col = "blue", lwd = 2, lty=1)
lines(density(data[,6]), xlim=c(0,100), col = "blue", lwd = 2, lty=2)
legend("topleft",legend =c("Prior1Data1","Prior1Data2","Prior2Data1","Prior2Data2","Prior3Data1","Prior3Data2"),
       col = c("red", "red", "green", "green", "blue", "blue"),
       lty = c(1,2,1,2,1,2), lwd = 2)

# Use posterior distribution to get new predictive distribution of a new person from this population 
# Use data 1 and prior 3
# Simulate the mean parameter and the precision parameter 10000 
data2 <- rt(n=100000,df=2*Prior3["Post_Alpha",1])
data2 <- data2*sqrt(Prior3["Post_Beta",1]/Prior3["Post_Theta",1]/Prior3["Post_Alpha",1])+Prior3["Post_Mean",1]
data2 <- cbind(mean1=data2, tau1=rgamma(n=100000,shape=Prior3["Post_Alpha",1],rate=Prior3["Post_Beta",1]))
data2 <- cbind(data2, Height=apply(data2,c(1),function(x){rnorm(1,mean=x[1],sd=1/sqrt(x[2]))}))
plot(density(data2[,"Height"]), main="Predicitve Distribution for new observations", xlab="Height",ylab="Density")
summary(data2[,"Height"])
data2 <- rt(n=10000,df=2*Prior3["Post_Alpha",1])
data2 <- data2*sqrt(Prior3["Post_Beta",1]/Prior3["Post_Theta",1]/Prior3["Post_Alpha",1])+Prior3["Post_Mean",1]
data2 <- cbind(mean1=data2, tau1=rgamma(n=10000,shape=Prior3["Post_Alpha",1],rate=Prior3["Post_Beta",1]))
data2 <- cbind(data2, Height=apply(data2,c(1),function(x){rnorm(1,mean=x[1],sd=1/sqrt(x[2]))}))
plot(density(data2[,"Height"]), xlim=c(-50,150), ylim=c(0,0.005), main="Predicitve Distribution for
     new observations", xlab="Height", ylab="Density")
summary(data2[,"Height"])
sd(data2[,"Height"])
mean(data2[,"Height"])+1.96*sd(data2[,"Height"])
mean(data2[,"Height"])-1.96*sd(data2[,"Height"])



# The third project
# Basic, bare-bones, fixed effect meta analysis
# If do a randomized clinical trial with 209 subjects and you get a value for (Y1; S1) of (-1.82; 0.21)
# The posterior belief in theta (average amount that the drug lowers the value of HbA1c) 
mean_0 <- 0
tau_0 <- 0
n <- 1
x_bar <- -1.82
Si <- 0.21
tau <- 1/(Si)^2
mean_1 <- (tau_0*mean_0+n*tau*x_bar)/(tau_0+n*tau)
tau_1 <- tau_0 + n*tau
S_1 <- 1/sqrt(tau_1)
mean_1
tau_1
S_1
qnorm(c(0.025,0.975),mean=mean_1,sd=S_1)

# New Info: 79 subjects and got a value for (Y2; S2) of (-1.02; 0.28). Starting with the prior belief from
# previous info and update the posterior belief using new data
mean_1 <- -1.82
tau_1 <- 22.67574
n <- 1
x_bar <- -1.02
Si <- 0.28
tau <- 1/(Si)^2
mean_2 <- (tau_1*mean_1+n*tau*x_bar)/(tau_1+n*tau)
tau_2 <- tau_1 + n*tau
S_2 <- 1/sqrt(tau_2)
mean_2
tau_2
S_2
qnorm(c(0.025,0.975),mean=mean_2,sd=S_2)

# Posterior belief with new info (without knowing the old info)
mean_0 <- 0
tau_0 <- 0
n <- 1
x_bar <- -1.02
tau <- 1/(0.28)^2
mean_1 <- (tau_0*mean_0+n*tau*x_bar)/(tau_0+n*tau)
tau_1 <- tau_0 + n*tau
S_1 <- 1/sqrt(tau_1)
mean_1
tau_1
S_1
qnorm(c(0.025,0.975),mean=mean_1,sd=S_1)
# Posterior belief with new info (with knowing the old info)
mean_0 <- -1.02
tau_0 <- 12.7551
n <- 1
x_bar <- -1.82
tau <- 1/(0.21)^2
mean_1 <- (tau_0*mean_0+n*tau*x_bar)/(tau_0+n*tau)
tau_1 <- tau_0 + n*tau
S_1 <- 1/sqrt(tau_1)
mean_1
tau_1
S_1
qnorm(c(0.025,0.975),mean=mean_1,sd=S_1)