### This script runs 1000 null EWAS permuting case control status 


# Setting up
setwd("../../Data")
library(doParallel)
load("US_Betas_Pheno_Filtered.rdata")
ls()
#"betas"    "pheno"


# Open 30 processors
cl<-makeCluster(30)
registerDoParallel()
clusterExport(cl,betas)

# Define function to run a single linear model testing for association
DMP_model<-function(row, status, pheno){
	model<-lm(row ~ status + pheno$confage + factor(pheno$nsex) + pheno$CD8T + pheno$CD4T + pheno$NK + pheno$Bcell + pheno$Mono + pheno$Gran)
	return(summary(model)$coefficients[2,4] )
}

# Create a binary status vector (i.e Case/Control) to permute in loop
status<-rep(0,nrow(pheno))
status[1:(nrow(pheno)/2)]<-1

# Run 1000 EWAS - saving p-values in 10 batches of 100 permutations
for(num in 1:10){
	nPerms<-100
	res<-matrix(data = NA, nrow = nrow(betas), ncol = nPerms)
	rownames(res)<-rownames(betas)

	for(i in 1:nPerms){
		status<-sample(status)
		res[,i]<-parRapply(cl, betas, DMP_model, status, pheno)

	}

	save(res, file = paste("EWASPermutations/EWASPermutations", nPerms, "Num", num, ".rdata", sep = ""))
}
