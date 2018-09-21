### This script performs power calulations for every DNA methylation site for a range of sample sizes and mean differences
# Then for each condition calculates the percentage of sites with increasing power


# Setting up
setwd("../Results")
library(pwr)
load("../../Data/US_Betas_Pheno_Filtered.rdata")
ls()
#[1] 


# Calculating Cohen's d for each site for two options (mean difference / standard deviation)
sds<-apply(betas, 1, sd)
d2<-0.02/sds
d5<-0.05/sds

#Saving the standard deviations to load into the shiny app
save(sds, file="SiteSDs.rdat")

# Choosing sample sizes to test
n<-c(100,200,500,1000,2000,5000)

# Creating object to save the power for each site - for each sample size
powers_d2<-vector("list", length(n))
powers_d5<-vector("list", length(n))

# Creating a vector of powers from 0 to 1
inc<-0.01
x<-seq(0,1,inc)

# Creating object to save the percentage of sites with each power in x - for each sample size
percpowered_d2<-vector("list", length(n))
names(percpowered_d2)<-n
percpowered_d5<-vector("list", length(n))
names(percpowered_d5)<-n


# Looping through each sample size in n (h) and performing power calculations for each site (i)
# then looping through each power in x (j) and calculating the percentage of sites with > that power
for(h in 1:length(n)){

	powers_tmp<-rep(NA, length(sds))
	names(powers_tmp)<-rownames(betas)
	
	for(i in 1:length(sds)){
		ptest<-pwr.t.test(n=n[h]/2, d=d2[i], sig.level=9.42e-8, type="two.sample", alternative = "two.sided")
		powers_tmp[i]<-unlist(ptest$power)
	}
	
	powers_d2[[h]]<-powers_tmp

	npowered_tmp<-rep(NA, length(x))
	names(npowered_tmp)<-paste(">",x, sep="")

	for(j in 1:length(x)){
		npowered_tmp[j]<-sum(powers_d2[[h]] > x[j])
	}
	
	percpowered_d2[[h]]<-npowered_tmp/nrow(betas)
}

# Repeating for d5
for(h in 1:length(n)){

	powers_tmp<-rep(NA, length(sds))
	names(powers_tmp)<-rownames(betas)
	
	for(i in 1:length(sds)){
		ptest<-pwr.t.test(n=n[h]/2, d=d5[i], sig.level=9.42e-8, type="two.sample", alternative = "two.sided")
		powers_tmp[i]<-unlist(ptest$power)
	}
	
	powers_d5[[h]]<-powers_tmp

	npowered_tmp<-rep(NA, length(x))
	names(npowered_tmp)<-paste(">",x, sep="")

	for(j in 1:length(x)){
		npowered_tmp[j]<-sum(powers_d2[[h]] > x[j])
	}
	
	percpowered_d5[[h]]<-npowered_tmp/nrow(betas)
}


# Saving
save(percpowered_d2, percpowered_d5, n, x, sds, file="PercentPowered.rdata")

