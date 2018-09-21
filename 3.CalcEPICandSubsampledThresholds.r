### This script loads the EWAS permutations and finds the significance threshold 
# for the full EPIC array and for subsamples of sites (in order to extrapolate later)

# Setting up
setwd("../Results")

filenames<-paste("../../Data/EWASPermutations/EWASPermutations100Num", seq(1,10,1), ".rdata", sep="")
allres<-matrix(ncol=0, nrow=804826)

for (f in 1:10){
	filename<-filenames[f]
	load(filename)
	allres<-cbind(allres, res)
}
dim(allres)
#[1] 804826   1000


# Finding the threshold for the full EPIC array (5th percentile of minimum p-values)
FullSimMinP<-apply(allres, 2, min)
quantileFullSimMinP<-quantile(FullSimMinP, 0.05)

# Choosing proportions to subsample
percents<-seq(from = 0.05,  to = 0.95, by = 0.1)
nprobes<-round(percents*nrow(allres))

# Creating an empty matrix to save the thresholds from all permutations
nRepeats<-100
allThres<-matrix(NA, nrow = length(percents), ncol=nRepeats)
rownames(allThres)<-percents
colnames(allThres)<-paste("Rep", seq(1, nRepeats, 1), "Pt")

# Looping through proportion subsampled (i), repeating 100 times (k), and finding the threshold (5th percentile of minimum p-values)
for(i in 1:length(nprobes)){
	print(paste0("Subsampling procedure being performed for ", percents[i]*100, "% of probes"))

	for(k in 1:nRepeats){
		print(paste("Repeat", k))
		
		#repIndex is the n randomly sampled rows
		repIndex<-sample(1:nrow(allres), nprobes[i])
		
		#minP takes the minimum P value in each of the 1000 EWAS permutations
		minP<-apply(allres[repIndex,], 2, min)
		
		#the record the takes the 5th percentile point of minP
		allThres[i,k]<-quantile(minP, 0.05)
	}
}

# Saving all results
save(allThres, file="AllSubsampledThresholds.rdat")

# Converting all thresholds to an effective number of independant tests (using Bonferroni in reverse)
allNTests<-0.05/allThres

# For each proportion subsampled, finding the mean threshold, mean number of tests, and 95% CIs
Thres<-matrix(data = NA, nrow = length(nprobes), ncol = 8)
colnames(Thres)<-c("x", "nSites", "Pt", "p_lt", "p_ut", "m", "m_lt", "m_ut")
Thres[,1]<-percents
Thres[,2]<-nprobes
Thres[,3]<-apply(allThres, 1, mean)
Thres[,4]<-apply(allThres, 1, quantile, 0.025)
Thres[,5]<-apply(allThres, 1, quantile, 0.975)
Thres[,6]<-0.05/Thres[,3]
Thres[,7]<-apply(allNTests, 1, quantile, 0.025)
Thres[,8]<-apply(allNTests, 1, quantile, 0.975)

# Adding the full EPIC threshold as the final row and the % of tests that are independant as the final column
Thres<-rbind(Thres, c(1, nrow(allres), rep(quantileFullSimMinP, 3), rep(0.05/quantileFullSimMinP, 3)))
Thres<-cbind(Thres, round(100*Thres[,6]/Thres[,2] ,1))
colnames(Thres)<-c("x", "nSites", "Pt", "p_lt", "p_ut", "m", "m_lt", "m_ut", "PC")

# Save results
write.csv(Thres, file="PermuationsPvalResults_WithCI.csv", row.names=F)


