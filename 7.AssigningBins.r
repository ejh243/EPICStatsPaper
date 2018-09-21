### This script calculates the mean and variability of the DNA methylation sites
# It then assigns sites to bins on these values and on their average EWAS ranking
# This can be loaded and combined with TestingLinearRegressionAssumptionsResults.csv later


# Setting up
setwd("../Results")
load("US_Betas_Pheno_Filtered.rdata")
avrank<-read.csv("AvPermuationsRank.csv", row.names=1)
ls()
#[1] "avrank" "betas"  "pheno"
identical(rownames(avrank), rownames(betas))
#[1] TRUE
avrank<-avrank$x


# Calculate mean sd and diffs (range of middle 80th percentile P90-910)
means<-rowMeans(betas)
sds<-apply(betas, 1, sd)
lt<-apply(betas, 1, quantile, 0.1, na.rm=TRUE)
ut<-apply(betas, 1, quantile, 0.9, na.rm=TRUE)
diffs<-(ut-lt)


# Writing to a table
res<-as.data.frame(cbind(means, sds, diffs, avrank))
rownames(res)<-rownames(betas)


# Plotting histograms of these values to help decide boundarys
pdf("../Figures/HistogramsOfBins.pdf")
par(mfrow=c(2,2))
hist(means, main="Histogram Of Mean Methylation", xlab="Mean DNA Methylation")
hist(sds, main="Histogram Of Standard Deviations", xlab="Standard Deviation")
hist(diffs, main="Histogram Of Range of Middle 80%", xlab="Range of Middle 80%")
hist(avrank$, main="Histogram Of Mean Rank in 1000 EWAS", xlab="Average EWAS Ranking")
dev.off()


# Assigning sites to 10 bins by their mean methylation
res$MethBin<-""
for (b in 1:10){
	probesInBin <- which(means > (0.1*b-0.1) & means <= (0.1*b))
	res$MethBin[probesInBin] <- paste(0.1*b-0.1, 0.1*b, sep="-")
}


# Assigning sites to 11 bins by their standard deviation (one bin needed for the tail)
res$SDBin<-""
for (b in 1:10){
    probesInBin <- which(sds > (0.01*(b-1)) & sds <= (0.01*b))
	res$SDBin[probesInBin] <- paste((0.01*(b-1)), (0.01*b), sep="-")
}
res$SDBin[which(sds > 0.1)] <- "> 0.1"


# Assigning sites to 3 bins by their range of middle 80% 
res$VarBin<-""
res$VarBin[which(diffs <= 0.05)]<-"< 0.05"
res$VarBin[which(diffs > 0.05 & diffs <= 0.1)]<-"0.05-0.1"
res$VarBin[which(diffs > 0.1)]<-"> 0.1"


# Assigning sites to 9 bins by their average ranking (bin needed for the tail on either end)
res$RankBin<-""
for (b in 1:7){
	probesInBin <- which(avrank> (5000*b+380000) & avrank <= (5000*b+385000))
	res$RankBin[probesInBin] <- paste(5000*b+380000, 5000*b+385000, sep="-")
}
res$RankBin[which(avrank > 420000)]<-">420000"
res$RankBin[which(avrank < 385000)]<-"<385000"


# Saving the extended results table
write.csv(res, file="SitesMeanVariabilityRankBins.csv")
