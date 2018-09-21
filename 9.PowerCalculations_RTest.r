### This script performs power calulations for a continuous phenotype with a range of expected correlations
# Then for each condition find the sample size needed for increasing power


# Setting up
setwd("../Results")
library(pwr)


# Caculating powers for each correlation of interest (increasing n (x-axis) and plotting power (y-axis))
n<-seq(4,5000,1) #n must be greater than 4 for calcuation
cors<-c(0.05, 0.1, 0.2)

powers_all<-vector("list", length(cors))
names(powers_all)<-paste("r =", cors)

for(i in 1:length(cors)){
	powers<-rep(NA, length(n))
	for(j in 1:length(n)){
		powers[j]<-pwr.r.test(n=n[j], r=cors[i], sig.level=9.42e-8, power=NULL)$power
	}
	powers_all[[i]]<-powers
}


# Saving
save(powers_all, n, cors, file="CorPowers.rdata")