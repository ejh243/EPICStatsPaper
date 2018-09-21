### This script takes the 1000 EWAS permutations and for each site calculates its average ranking


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



# Looping through each EWAS and ranking the sites by their p-value
ranks<-matrix(NA, ncol=ncol(allres), nrow=nrow(allres))
rownames(ranks)<-rownames(allres)

for(i in 1:ncol(allres)){
	perm<-allres[,i]
	ranks[,i]<-order(perm)
}

# Finding average rank for each site
avrank<-rowMeans(ranks)

# Saving 
write.csv(avrank, "AvPermuationsRank.csv")