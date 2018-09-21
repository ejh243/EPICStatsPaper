### This script loads the subsampled EPIC array thresholds and fits a Monod model to extrapolated to infinite density 

# Setting up
setwd("../Results")
library(FME)

Obs<-read.csv("PermuationsPvalResults_WithCI.csv")

# Defining the Monod Model as a function of x and p
    #x is the actual number of probes
    #p is a vector of parameters p=(u, k)
    #f returned is the y-axis (effective number of independant tests)
    
MonodModel <- function(x, p){
	f <- p[1]*x/(x+p[2])
	return(data.frame(x, f))
}

# Defining a function to calculate the residuals - so that we can call it into modFit to minimise
Residuals <- function(p){
	r <- Obs$m - MonodModel(Obs$x, p)$f
}

# Using modfit to minimise the residuals by changing the parameters p
P <- modFit(f = Residuals, p = c(0.1, 1))
sP <- summary(P)

# The results of the optimised parameters are 
(optimisedP<-P$par)
#[1] 5.803113e+06 9.875127e+00

#The asymptote of the Monod Model occurs at u, so our maximum number of independant tests is: 
(maxM=optimisedP[1])
#[1] 5803067
#With bonferroni correction this gives a genome-wide significance threshold of:
(maxY=0.05/maxM)
#[1] 8.616134e-9


# Applying the monod model to the lower and upper 95% CI limits of m
low.Residuals <- function(p){
	r <- Obs$m_lt - MonodModel(Obs$x, p)$f
}
low.P <- modFit(f = low.Residuals, p = c(0.1, 1))
(low.optimisedP<-low.P$par)
#[1] 4.025465e+13 8.081151e+07
(low.maxM<-low.optimisedP[1])
#[1] 3.360691e+13
(low.maxY<-0.05/low.maxM)
#[1] 1.487789e-15

up.Residuals <- function(p){
	r <- Obs$m_ut - MonodModel(Obs$x, p)$f
}
up.P <- modFit(f = up.Residuals, p = c(0.1, 1))
(up.optimisedP<-up.P$par)
#[1] 1.686186e+06 2.046253e+00
(up.maxM<-up.optimisedP[1])
#[1] 1686227
(up.maxY<-0.05/up.maxM)
#[1] 2.965199e-08


# Saving parameters from fitted models
save(optimisedP, low.optimisedP, up.optimisedP, file="MonodParameters.rdata")

