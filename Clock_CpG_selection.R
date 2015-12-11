########################################################
##	This script contains the statistical  	  #
##          pipeline used to select the		  #
## set of 67 BE clock CpGs used in the Bayesian model   #
########################################################


## Use functions from package 'minfi' for analysis of methylation data
source("http://bioconductor.org/biocLite.R")
biocLite("methylumi")
biocLite("minfi")
biocLite("siggenes")


## Load these packages in package Manager
library(methylumi)
library(minfi)
library(siggenes)

## Load longitudinal data sets (SWAN normalized and filtered) to be used for
## 1) selecting significantly drifting CpGs among 450K 
## 2) deriving drift rates via regression to be used in relaxed molecular clock
########################################################
## 1) The following data set D1 was used for SELECTION of marker set as the bigger of the two serial sets
########################################################

load('gset.nfilt.rda',verbose=T)


mval = getM(gset)
pheno=pData(gset)
age = as.character(pheno$Age.at.Biopsy); age = as.numeric(age)
sampleID = dimnames(mval)[[1]]
## Ten patients in D1 with between 2-5 serial samples
indivID = c(1,1,2,2,3,3,3,3,4,4,5,5,6,6,6,6,6,7,7,7,8,8,8,8,9,9,9,10,10,11,rep(0,18))
mval = mval[,1:30]
## ===========================================================================================
## longitudinal regressions
## Require population drift in longitudinal data, along with indiv significant drift and obtain slopes

m.test=t(mval)
dat = m.test   
age.long = age[1:30]
indivID=indivID[1:30]
M = ncol(dat)
age.reg = age.long

indiv = (1:10); np = length(indiv); Age = numeric(np)  

for( j in indiv) {
    ids  = (1:30)[indivID==j]
    Age[j] = min(age.long[ids])
  for (m in 1:M) {
  	## indiv regressions for normalization
    aux = lm(dat[ids,m] ~ age.long[ids])
    intercept = aux$coef[1]; slope = aux$coef[2]
    dat[ids,m] = dat[ids,m]-(intercept+slope*Age[j])
  }
}
for(i in indiv) {
  age.reg[indivID==i] = age.long[indivID==i] - min(age.long[indivID==i])
}
age.reg = age.reg[-30]; dat = dat[-30,]

## run regressions on adjusted/normalized data
pval = intercept = slope = trend = p.trend = numeric(M)
for (m in 1:M) {
    aux1 = summary(lm(dat[,m] ~ age.reg))
    aux2 = summary(lm(m.test[,m] ~ age.long))
    pval[m] = aux1$coef[2,4]; intercept[m] = aux1$coef[1,1]; slope[m] = aux1$coef[2,1]; trend[m] = aux2$coef[2,1]; p.trend[m] = aux2$coef[2,4]

}

b_sd = numeric(M)
for (m in 1:M) {
    aux1 = summary(lm(dat[,idstest[m]] ~ age.reg))
    b_sd[m] = aux1$coef[2,2]; 
}

pi0<-pi0.est(pval)$p0
tmp = qvalue.cal(pval,p0=pi0); summary(tmp)
## FDR q-value of .2 to loosely control for false discoveries in signifcant longitudinal drift
signif = .2

lnegCpGs = dimnames(m.test)[[2]][tmp<signif & slope < 0  & p.trend < 0.01]
lposCpGs = dimnames(m.test)[[2]][tmp<signif & slope > 0  & p.trend < 0.01]
longCpGs = dimnames(m.test)[[2]][tmp<signif & p.trend < 0.01]

lposSlopes = slope[tmp<signif & slope > 0 & p.trend < 0.01]
lposTrend =  trend[tmp<signif & slope > 0 & p.trend < 0.01]
lSlopes = slope[tmp<signif & p.trend < 0.01]


## ===========================================================================================
## ANCOVA using cross-sectional data
## Require differential drift between matched SQ (Normal Squamous) & BE in population drift in cross-sectional data
########################################################
## The following data set D2 was also used for SELECTION of marker set in conjunction with the longitudinal set D1
########################################################
load("SQ.gset.nfilt.rda")
load("BE.gset.nfilt.rda")

sqM <- getM(SQ.gset.nfilt)
beM <- getM(BE.gset.nfilt)

colnames(sqM) <- SQ.gset.nfilt$Sample_Name
colnames(beM) <- BE.gset.nfilt$Sample_Name

pd <- pData(SQ.gset.nfilt)
age = pd$Age.at.biopsy
Age.BE = Age.SQ = age

## CpG names
CpGs = dimnames(sqM)[[1]]

ido = 1:length(CpGs)

#### ANCOVA SQ-BE ######

age = Age.SQ  

M = nrow(sqM) 
ido = ids = 1:M 
rateSQ = rateBE = pvalSQ = pvalBE = pval.intcptBE = epsBE = numeric(M)

mstart=1;mend=M
  for (m in mstart:mend) {
    mp = ids[m]  #.drift[m]
    tmp = data.frame(age=c(age,age),hist=c(rep("SQ",30),rep("BE",30)),methyl=c(sqM[mp,],beM[mp,]))

    mfit = lm(methyl ~ age*hist, data=tmp)
    aux = summary(mfit)

    rateSQ[m] = aux$coef[2,1]; rateBE[m] = aux$coef[4,1]+rateSQ[m]
    pvalSQ[m] = aux$coef[2,4]; pvalBE[m] = aux$coef[4,4]
    epsBE[m] = aux$coef[1,1]+aux$coef[3,1]
    pval.intcptBE[m] = aux$coef[3,4]
  }

aov.SQ = list(rate=rateSQ,pval=pvalSQ)
aov.BE = list(rate=rateBE,pval=pvalBE,eps=epsBE,pval.intcpt=pval.intcptBE)

## Match CpG names of serial and cross-sectional data sets
ido.lpos = na.omit(match(lposCpGs,CpGs))
ido.lneg = na.omit(match(lnegCpGs,CpGs))

## Eventually used the up drifting CpGs as there were many more with this behavior (see Main Text)
ido.long = ido.lpos 

mset.max = apply(sqM,1,max); mset.min = apply(sqM,1,min)
mtheta = logit2(.25)
confl = 0.05; intcptl=0.05
##for pos
ido.hypo = ido.long[mset.max[ido.long] < mtheta]
ido.diff =  ido.long[mset.max[ido.long] < mtheta & (aov.BE$pval[ido.long] < confl)] 


#### PCA analysis ####
iset1=ido.diff
ibeM = matrix(0,nrow=length(iset1),ncol=length(age))
for(i in 1:length(age)) {ibeM[,i]= aov.BE$eps[iset1]+aov.BE$rate[iset1]*age[i]}
q =(beM[iset1,]-ibeM) 
eig = eigen(cor(q)); pc = q %*% eig$vectors

### For publication, excluded LHS of PCA that all have BE pop slopes < BE individual slopes
up_drift_rates=aov.BE$rate[iset1]
ids1 = which(up_drift_rates>aov.SQ$rate[iset1])   
CpGs.D1.67=NULL
CpGs.D1.67$cpgn = CpGs[ido.diff[ids1]] 

########################################################
## 2)  We can derive drift RATES from data set D1 or an independent validation set DV.
## Here I will save the D1 rates to use in the BE clock model prior distribution 
########################################################
## Note: could use same regression scheme above for DV and obtain drift rates to use in prior

CpGs.D1.67$b.slope = up_drift_rates[ids1]
save(CpGs.D1.67, "CpGs.D1.67.rda")

## This completes the script for "cpgn" marker name selection (67 total in above pipeline)
## + drift rates for each serived from regression over longitudinal data, all kept in CpGs.D1.67




