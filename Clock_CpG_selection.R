########################################################
#	This script contains the statistical  	  #
#          pipeline used to select the		  #
# set of 67 BE clock CpGs used in the Bayesian model   #
########################################################


### Use functions from package 'minfi' for analysis of methylation data
source("http://bioconductor.org/biocLite.R")
biocLite("methylumi")
biocLite("minfi")
biocLite("siggenes")