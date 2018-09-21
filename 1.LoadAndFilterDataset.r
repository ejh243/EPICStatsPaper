### This script is to load the already QC'ed Understanding Society dataset, and filter DNA methylation sites

# Setting up
setwd("../../Data")
load("../../US_Betas_Pheno.rda")
ls()
#[1] "dat"  "pheno"
identical(colnames(dat), pheno$barcode)
#[1] TRUE
betas<-dat
dim(betas)
#[1] 857130   1175


# Remove Sample IDs
pheno<-pheno[,c(5:16)] 

# Loading reference files
crosshyb<-read.table("/mnt/data1/EPIC_reference/CrossHydridisingProbes_McCartney.txt", stringsAsFactors = FALSE)
snpProbes<-read.table("/mnt/data1/EPIC_reference/SNPProbes_McCartney.txt", stringsAsFactors = FALSE, header = TRUE)
snpProbes<-snpProbes[which(snpProbes$EUR_AF >= 0.05 & snpProbes$EUR_AF <= 0.95),]
 
# Removing cross hybridising probes and SNPs
betas<-betas[!(rownames(betas) %in% crosshyb[,1]), ]
betas<-betas[!(rownames(betas) %in% unique(snpProbes$IlmnID)), ]
betas<-betas[-grep("rs", rownames(betas)),]

# Using the EPIC manifest to remove sites annotated to the Y chromosome 
epicManifest<-read.csv("/mnt/data1/EPIC_reference/MethylationEPIC_v-1-0_B2.csv", skip = 7)
betas<-betas[!(rownames(betas) %in% epicManifest$IlmnID[epicManifest$CHR == "Y"]),]


# Saving filtered dataset
dim(betas)
[1] 804826    1175
save(betas, pheno, file="US_Betas_Pheno_Filtered.rdata")