### This script runs an EWAS of age, testing each site's model for violations of the assumptions of linear regression

# Setting up
setwd("../Results")
library(gvlma)
load("US_Betas_Pheno_Filtered.rdata")
ls()
#"betas"    "pheno"


# Creating an empty matrix to save results
res<-matrix(data = NA, nrow = nrow(betas), ncol = 5)
colnames(res)<-c("GlobalTestP", "SkewnessTestP", "KurtosisTestP", "LinkFunctionTestP", "HeteroskedasticityTestP")
rownames(res)<-rownames(betas)


# Running an EWAS of age on this dataset, and using gvlma to test each regression model
for(i in 1:nrow(betas)){
	model<-lm(betas[i,] ~ as.numeric(pheno$confage) + as.factor(pheno$nsex) + as.factor(pheno$MethArray_Chip) + pheno$CD8T + pheno$CD4T + pheno$NK + pheno$Bcell + pheno$Mono + pheno$Gran)
	gvmodel <- gvlma(model)
	res[i,1]<-gvmodel$GlobalTest$GlobalStat4$pvalue
	res[i,2]<-gvmodel$GlobalTest$DirectionalStat1$pvalue
	res[i,3]<-gvmodel$GlobalTest$DirectionalStat2$pvalue
	res[i,4]<-gvmodel$GlobalTest$DirectionalStat3$pvalue
	res[i,5]<-gvmodel$GlobalTest$DirectionalStat4$pvalue
}

# Save results
write.csv(res, "TestingLinearRegressionAssumptionsResults.csv")