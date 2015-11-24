## This is the BE clock Gibbs sampler

gibbs <- function(x,L=100,K=100,skip=10,plot=T) {
  if (!is.list(x)) {
	cat("x is not a list! see help file", "\n")
	return()
      }
  names(x)[1] <- "label"
  names(x)[2] <- "est"
  names(x)[3] <- "low"
  names(x)[4] <- "upp"

  npar <- length(x$est)

  if (npar <= 0) stop('no. of parameters must be >= 1')

  # *** initialize graphical output
  if(plot==TRUE) {
      # par(mfrow=c(npar+1,1),bg="grey")
      par(mfrow=c(8,1),bg="grey")
  }

  ## starting point
  y = x$est
  ## data-specific setup
  ## provide: ids, ID, M, a 
  xj = beM[ids,ID]; s = y[M+1]
  
  x.mon <- matrix(0,ncol=npar,nrow=L)

  for (l in 1:L) {
    
    ## full conditionals for b_j's: x[1:M]
    for(m in 1:M) {
      xm = xj[m]
      SD = y[M+2]
      mb = mean_b_i[m] # mean b prior
      sd = sd_b_i[m]   # sd b prior
      r = (xm - alpha_SQ[m] - b_SQ[m]*s)/(a-s)
      denom = (SD^2+((a-s)*sd)^2)
      b.mean = (SD*SD*mb+sd*sd*r*(a-s)^2)/denom
      b.sigm = sd*SD/sqrt(denom)
      y[m] = rnorm(1,mean = b.mean, sd = b.sigm)
    }

    ## full conditional for s: x[M+1]
    test=a+1
    while( test >a || test < 0){
      D = sum((y[1:M]-b_SQ)*(y[1:M]-b_SQ))
      g = alpha_SQ + y[1:M]*a - xj
      E = sum(g*(y[1:M]-b_SQ))
      test = rnorm(1, mean = E/D, sd = y[M+2]/sqrt(D))
    }
    y[M+1] =test 
    s = y[M+1]

    ## full conditional for sigma_BE: x[M+2]
    gam_1=0.6075155
    gam_2=0.1173199
    gam.shape = 0.5*M
    gam.rate = 0.5*sum((xj-alpha_SQ-b_SQ*s-y[1:M]*(a-s))^2)
    dum = rgamma(1, shape=gam_1+gam.shape, rate=gam_2+gam.rate)
    y[M+2] = 1/sqrt(dum)

    x.mon[l,] <- y
    # PLOTTING OF RUNS ##
    if(l%%(100*skip) == 0) {
      if(plot==TRUE) {
        if(l < K) {brncol <- 3} else {brncol <- 2}
        n.skip.1 <- seq(skip,l,skip)
        n.skip.2 <- seq(skip,min(l,K),skip)
        
        par(mar=c(0, 5, 0.6, 4) + 0.1)
        plot(x.mon[n.skip.1,1], type='l', xlab = " ", ylab = x$label[1],col=2, main = ID)
        lines(x.mon[n.skip.2,1],col=brncol) # burn-in cycles

        # for (i in 1:(npar-1)) {
        for (i in 2:6) {
          par(mar=c(0, 5, 0, 4) + 0.1)
          plot(x.mon[n.skip.1,i], type='l', xlab = " ", ylab = x$label[i], col=2)
          lines(x.mon[n.skip.2,i],col=brncol) #pilot cycles
        }
        par(mar=c(0, 5, 0, 4) + 0.1)
        plot(x.mon[n.skip.1,npar-1], type='l', xlab = " ", xaxt='n', ylab = x$label[npar-1], col=2)
        lines(x.mon[n.skip.2,npar-1],col=brncol) #pilot cycles
        
        par(mar=c(0.1, 5, 0, 4) + 0.1)
        plot(x.mon[n.skip.1,npar], type='l', xlab = " ", xaxt='n', ylab = x$label[npar], col=2)
        lines(x.mon[n.skip.2,npar],col=brncol) #pilot cycles
      }
    }
  }
  return(x.mon)
}

#######################
####    SETUP      ####
#######################

#source("http://bioconductor.org/biocLite.R")
#biocLite("methylumi")
#biocLite("minfi")

### Load these packages in package Manager
library(methylumi)
library(minfi)
library(Bhat)
## Load cross-sectional data
load("./Illum450kEsophagus20131205/data/SQ.gset.nfilt.rda")
load("./Illum450kEsophagus20131205/data/BE.gset.nfilt.rda")
load('./Illum450kEsophagus20131205/data/FBE89set.rda')
load('./Illum450kEsophagus20131205/data/BEwHGD89set.rda')
load('./Illum450kEsophagus20131205/data/FBE70set_copy.rda')
load('./Illum450kEsophagus20131205/data/BEwHGD70set_copy.rda')
load('./Illum450kEsophagus20131205/CpGs.D1.70.rda')


CpGs.D1.70 = CpGs.D1.70[id_cpgn,]
#################
## 1) ALL data ##
#################
sqM <- getM(SQ.gset.nfilt)
#bval = getBeta(SQ.gset.nfilt)
#sqM = logit(bval,0,1)
pd <- pData(SQ.gset.nfilt)
beM <- getM(BE.gset.nfilt)
#bval = getBeta(BE.gset.nfilt)
#beM = logit(bval,0,1)
colnames(sqM) <- SQ.gset.nfilt$Sample_Name
colnames(beM) <- BE.gset.nfilt$Sample_Name
M = 67
age = pd$Age.at.biopsy
ids= na.omit(match(CpGs.D1.70$cpgn[1:M],dimnames(beM)[[1]]))
id = ids
slope.BE = intercept.BE = slope.SQ = intercept.SQ = numeric(length(id))
k=1
for(i in id) {
tmp = lm(beM[i,] ~ age) #; cat('slopes:',tmp$coef[2],slopes[k],'\n')
intercept.BE[k] = tmp$coef[1]; slope.BE[k] = tmp$coef[2]
tmp = lm(sqM[i,] ~ age) #; cat('slopes:',tmp$coef[2],slopes[k],'\n')
intercept.SQ[k] = tmp$coef[1]; slope.SQ[k] = tmp$coef[2]
k=k+1
}
########################
## 2) group specific ###
########################
## general population data
ages=pd$Age.at.biopsy
age=ages
#mean_b_i = rep(mean(CpGs.D1.70$b.slope),67)
#mean_b_i = b.slope_v
##  regressions of CCF with 9 patients, 67 markers
mean_b_i=rep(mean(slope.ccf),67)
#sd_b_i = rep(sd(CpGs.D1.70$b.slope),67)
#sd_b_i=b.sd_v
## regressions of CCF
sd_b_i =rep(sd(slope.ccf),67) 
N=length(age)
sum_matrixgen_joe = matrix(rep(0,5*N),ncol=N)


## fBE data
pd_total = pData(mset.swan)
FBE_id=which(pd_total$Tissue_Type=="FBE")
## From Andrew
pd_total$Age_at_biopsy[FBE_id[17:22]]=c(70,48,73,72,70,52)
Age.FBE=pd_total$Age_at_biopsy[FBE_id]
fbeMset = getM(mset.swan)[,FBE_id]
id_fbe_cpgs = match(CpGs.D1.67$cpgn,rownames(fbeMset))
ages=Age.FBE
age=ages
#id = which(!is.na(fbeMset[,1]))
ids=id_fbe_cpgs

intercept.SQ=intercept.SQ
slope.SQ=slope.SQ
beM = fbeMset
#mean_b_i = rep(mean(CpGs.D1.67$b.slope),67)
#mean_b_i = b.slope_v[id]
mean_b_i = rep(mean(slope.ccf),67)
#sd_b_i = rep(sd(CpGs.D1.67$b.slope),67)
#sd_b_i = b.sd_v[id]
sd_b_i = rep(sd(slope.ccf),67)
N=length(age)
sum_matrixFBE_ccf = matrix(rep(0,5*N),ncol=N)

## BE w HGD data
ages=Age.BEwHGD
age=ages
id = which(!is.na(mBEwHGD[,1]))
ids=id
intercept.SQ=intercept.SQ[id]
slope.SQ=slope.SQ[id]
beM = mBEwHGD
mean_b_i = CpGs.D1.70$b.slope[id]
sd_b_i = CpGs.D1.70$b.sd[id]
N=length(age)
sum_matrixBEwHGD = matrix(rep(0,5*N),ncol=N)
#################
## 3) ALL data ##
#################
runs = 100000
M=length(id)
alpha_SQ = intercept.SQ
b_SQ = slope.SQ
#slop_low = rep(-0.02,M)
#slop_upp = rep(0.16,M)
slop_low = rep(-0.04,M)
slop_upp = rep(0.3,M)
sum_matrix = matrix(rep(0,5*N),ncol=N)
out_temp = matrix(rep(0,runs*N),ncol=N)
for (i in 2:N){
  ## choose individual
  ID=i
  ## x is the list of parameters to be estimated (s_j, b_i, sigma_SQ, sigma_BE) I haven't added labels for those yet
  x = NULL
  a = age[i]
  # uniform distribution mean for onset times
  x$label = seq(1:(M+2))
  x$label[1:M]=paste(rep("b",M),1:M,sep='')
  x$label[M+1]="s"
  x$label[M+2]="sig"

  x$est[M+1] = a/2
  x$low[M+1] = 0 
  x$upp[M+1] = a
  ## Could provide/use guess at BE_start ages that's better than age/2, obtained from just drawing back with slope going to age/2
  ## start slopes at mean value 
  x$est[1:M] = mean_b_i
  x$low[1:M] = slop_low[ids]
  x$upp[1:M] = slop_upp[ids]
  ## estimates of standard deviations for dnorm, implicit uniform prior
  x$est[M+2] = 1.2 #sigma_BEstart
  x$low[M+2] = 0.1
  x$upp[M+2] = 3
  print(i)
  out = gibbs(x,L=runs,K=500,skip=10,plot=T)
  out_temp[,i] = out[,(M+1)]
  s = quantile(out[,(M+1)],probs=c(.025,.25,.5,.75,.975))
  sum_matrix[,i]=s 
}

##############################
##  4) group specific       ##
##############################
# general BE data
sum_matrixgen_joe= sum_matrix
outgen_joe = out_temp
# fBE data
sum_matrixFBE_ccf = sum_matrix
outFBE_ccf = out_temp
# BE w HGD data
sum_matrixBEwHGD = sum_matrix
out_BEwHGD=out_temp


#################################
## saving vectors with onsets ###
#################################
##### general BE data 
N=30
onsetsgen = rep(0,N)
sum_matrixgen_round = signif(sum_matrixgen_ccf,4)
for (i in 1:N){ onsetsgen[i]=paste(sum_matrixgen_round[3,i],"(",sum_matrixgen_round[1,i],",",sum_matrixgen_round[5,i],")")}
#write.table(onsetsgen,"BEonsetsgen_w95.txt",quote=F, row.names=F)
##### FBE data
N=22
onsetsFBE = rep(0,N)
sum_matrixFBE_round = signif(sum_matrixFBE_ccf,4)
for (i in 1:N){ onsetsFBE[i]=paste(sum_matrixFBE_round[3,i],"(",sum_matrixFBE_round[1,i],",",sum_matrixFBE_round[5,i],")")}
#write.table(onsetsFBE,"BEonsetsFBE_w95.txt",quote=F, row.names=F)
##### BE w HGD data
N=10
onsetsBEwHGD = rep(0,N)
sum_matrixgBEwHGD_round = signif(sum_matrixBEwHGD,4)
for (i in 1:N){ onsetsBEwHGD[i]=paste(sum_matrixBEwHGD_round[3,i],"(",sum_matrixBEwHGD_round[1,i],",",sum_matrixBEwHGD_round[5,i],")")}
#write.table(onsetsBEwHGD,"BEonsetsBEwHGD_w95.txt",quote=F, row.names=F)




#######################
## PLOTTING BOXPLOTS ##
#######################
id.agegen=order(pd$Age.at.biopsy)
id.ageBEwHGD=order(Age.BEwHGD)
id.agefBE=order(Age.FBE)

pd_sex_orderedgen=pd$gender[id.agegen]
pd_sex_orderedfBE = Sex.FBE[id.agefBE]
pd_sex_orderedBEwHGD = Sex.BEwHGD[id.ageBEwHGD]


##########################
##   group specific     ##
##########################
## general BE data
N=30
ages=Age.BE
id.age = order(ages)
pd_sex_ordered=pd_sex_orderedgen 
sum_matrixtemp = sum_matrixgen_ccf
## fBE data
N=21
ages=Age.FBE[-1]
id.age = order(ages)
pd_sex_ordered=pd_sex_orderedfBE[-1]
sum_matrixtemp = sum_matrixFBE_ccf[,-1]
## BE w HGD data
N=10
ages=Age.BEwHGD
id.age = order(ages)
pd_sex_ordered=pd_sex_orderedBEwHGD
sum_matrixtemp = sum_matrixBEwHGD

### PLOT ###
quartz("Quartz", width=6.9,height=3.4)
par(ps=10,mar=c(2,2,.5,.5), mgp=c(1,.25,0))
temp=boxplot(sum_matrixtemp[,id.age],names=as.character(1:length(ages)))
bxp(temp,ylim=c(0,90),xlab="Patient ID",ylab="BE onset age",cex.lab=1,las=2,axes=F, outline=F, border=fadeColor("forestgreen",'33'),boxfill=fadeColor("forestgreen",'22'))
abline(h=0,lty=2)
red_id= which(pd_sex_ordered=="Male")
red_id=c(red_id,which(pd_sex_ordered=="M"))
for (i in red_id){
    lines((i-.5):(i+.5),rep(ages[id.age[i]],2), col="red",lwd=1.5)
}
blue_id= which(pd_sex_ordered=="Female")
blue_id=c(blue_id,which(pd_sex_ordered=="F"))
for (i in blue_id){
    lines((i-.5):(i+.5),rep(ages[id.age[i]],2), col="blue",lwd=1.5)
}
legend('topleft',c("Age at biopsy: male BE patient","Age at biopsy: female BE patient"), col=c("red","blue"), lty=1, bty='n',cex=1,lwd=1.5)
axis(1, at=1:N,tck = -.009, lab=1:N, cex.axis=1, las=2)
axis(2, at=seq(0,90,10),tck = -.009, lab=seq(0,90,10), cex.axis=1, las=2,ylab="BE onset age")
points(1:N,sum_matrixtemp[3,id.age],col="forestgreen", pch=18,cex=.75)
box()


points(1:N,sum_matrixgen_ccf[3,id.age],col="purple", pch=18,cex=.75)
for (i in 1:30){
    lines(c((i-.25),(i+.25)),rep(sum_matrixgen_ccf[1,id.age[i]],2), col=fadeColor("purple",'44'),lwd=1)
    lines(c((i-.25),(i+.25)),rep(sum_matrixgen_ccf[5,id.age[i]],2), col=fadeColor("purple",'44'),lwd=1)
    lines(c(i,i),c(sum_matrixgen_ccf[5,id.age[i]],sum_matrixgen_ccf[1,id.age[i]]),, col=fadeColor("purple",'44'),lwd=1,lty=2)
    polygon(c(i-.4,i+.4,i+.4,i-.4),c(sum_matrixgen_ccf[2,id.age[i]],sum_matrixgen_ccf[2,id.age[i]],sum_matrixgen_ccf[4,id.age[i]],sum_matrixgen_ccf[4,id.age[i]]),col=fadeColor("purple",'33'),border=fadeColor("purple",'33'))
}


length(which(sum_matrixgen_ccf[3,] <= sum_matrixgen_joe[5,]& sum_matrixgen_ccf[3,] >= sum_matrixgen_joe[1,]))

x_val = seq(-.2,.2,.001)
quartz("Quartz", width=3.5,height=3.6)
par(ps=10,mar=c(2,2,.5,.5), mgp=c(1,.25,0))
i=1; plot(x_val,dnorm(x_val,CpGs.D1.70$b.slope[i], CpGs.D1.70$b.sd[i]),col=fadeColor("forestgreen","33"),lty=4, type='l',ylim=c(0,30), xlab="Drift rate (M-value/year)",ylab="Density",axes=F)
for (i in 2:67){ lines(x_val,dnorm(x_val,CpGs.D1.70$b.slope[i], CpGs.D1.70$b.sd[i]),col=fadeColor("forestgreen",'33'),lty=4)}

legend('topleft',c("D1 drift prior", "DV drift prior"),col=c("forestgreen","purple"),lwd=6, bty='n',cex=1)
for (i in 1:67){ lines(x_val,dnorm(x_val,slope.ccf[i], b_sd_ccf[i]),col=fadeColor("purple",'33'),lty=4)}
lines(x_val,dnorm(x_val,mean(CpGs.D1.70$b.slope), sd(CpGs.D1.70$b.slope)),col='forestgreen',lwd=6)
lines(x_val,dnorm(x_val,mean(slope.ccf), sd(slope.ccf)),col='purple',lwd=6)
abline(v=0,lty=2)
axis(1, tck = -.02, cex.axis=1, at= seq(-.2,.2,.05),cex.lab=1)
axis(2, tck = -.02,  cex.axis=1, at=seq(0,30,5),cex.lab=1)
box()

## Plot histogram and densities of estimated b_i values
tempout=out
N=1
M=67
#id = which(!is.na(mBEwHGD[,1]))
#bj_slopes= slop70[id]
bj_slopes= mean_b_i

quartz("Quartz", width=3.5,height=3.6)
par(ps=10,mar=c(2,2,.5,.5), mgp=c(1,.25,0))
plot(density(apply(tempout[,1:M],2,mean)),col="red",xlab="BE patient drift slopes",main="", axes=F, cex=1,cex.lab=1,ylim=c(0,50),lwd=1.5)
lines(density(bj_slopes),col="blue",lwd=1.5)
lines(seq(slop_low[1],slop_upp[1],.0001),dnorm(seq(slop_low[1],slop_upp[1],.0001),mean=mean(bj_slopes), sd=sd(bj_slopes)),lty=2,col="deepskyblue",lwd=1.5)
legend("topleft", c( "BE clock drift slopes", "Combined normal prior","MCMC means"), col=c("blue","deepskyblue", "red"), lty=c(1,2,1), bty="n", cex=1,lwd=1.5)
axis(1, tck = -.02, cex.axis=1,lab = seq(0,0.10, by=.01), at= seq(0,0.10, by=.01),cex.lab=1)
axis(2, at=seq(0,40,10),tck = -.02, lab=seq(0,40,10), cex.axis=1,ylab="Density")
box()



## BAYES FACTOR TESTING ##
N2= 30
N3= 22
N4= 10
burnend = runs
burnin = 0
sum1 = sum2 = sum3=rep(0,(burnend-burnin))
count=1
for ( k in 1:(burnend-burnin)){
  temp = sum((pd$Age.at.biopsy - outgen_ccf[k,1:N2])/(pd$Age.at.biopsy))
  sum1[count] = temp/N2
   temp = sum((Age.FBE[-3] - outFBE_ccf[k,(1:N3)[-3]])/(Age.FBE[-3]))
   sum2[count] = temp/(N3-1)
  # temp = sum((Age.BEwHGD - out_BEwHGD[k,1:N4])/(Age.BEwHGD))
  # sum3[count] = temp/N4
  count=count+1
}


count=1
sum1_p = sum2_p = sum3_p=rep(0,(burnend-burnin))

for ( k in 1:(burnend-burnin)){
  prior_temp = runif(30,0,pd$Age.at.biopsy)
  temp = sum((pd$Age.at.biopsy - prior_temp)/(pd$Age.at.biopsy))
  sum1_p[count] = temp/N2
  prior_temp = runif(21,0,Age.FBE)
  temp = sum((Age.FBE[-3] - prior_temp)/(Age.FBE[-3]))
  sum2_p[count] = temp/N3
  #prior_temp = runif(10,0,Age.BEwHGD)
  #temp = sum((Age.BEwHGD -prior_temp)/(Age.BEwHGD))
  #sum3_p[count] = temp/N4
  count=count+1
}

null = sum2>sum1
post_H0 = (length(null[null==TRUE]))/length(sum1)
post_H1 = 1 - post_H0

nullp = sum2_p>sum1_p
prior_H0 = (length(nullp[nullp==TRUE]))/length(sum2_p)
prior_H1 = 1 - prior_H0

BF2over1=(post_H0/prior_H0)/(post_H1/prior_H1)
print(BF2over1)

null = sum3>sum1
post_H0 = (length(null[null==TRUE]))/length(sum1)
post_H1 = 1 - post_H0

nullp = sum3_p>sum1_p
prior_H0 = (length(nullp[nullp==TRUE]))/length(sum3_p)
prior_H1 = 1 - prior_H0

BF3over1=(post_H0/prior_H0)/(post_H1/prior_H1)

null = sum3>sum2
post_H0 = (length(null[null==TRUE]))/length(sum2)
post_H1 = 1 - post_H0

nullp = sum3_p>sum2_p
prior_H0 = (length(nullp[nullp==TRUE]))/length(sum3_p)
prior_H1 = 1 - prior_H0

BF3over2=(post_H0/prior_H0)/(post_H1/prior_H1)


### Violin Plot 
quartz("Quartz", width=6.9,height=3.4)
par(ps=10,mfrow= c(1,2),mar=c(2,2,.5,.5), mgp=c(1,.25,0))
plot(0:1,0:1,type="n",xlim=c(0.5,3.5),ylim=c(0,90),axes=FALSE,ann=FALSE)

vioplot(apply(outgen[1:(burnend-burnin),id.agegen],2,mean),apply(out_fBE[1:(burnend-burnin),id.agefBE],2,mean),apply(out_BEwHGD[1:(burnend-burnin),id.ageBEwHGD],2,mean), names=c('BE','fBE','HGD/EAC'), col="aquamarine3", border="aquamarine3",ylim=c(0,90),add=T)
boxplot(pd$Age.at.biopsy,Age.FBE, Age.BEwHGD, add=T,border="darkgrey", axes=F)
axis(1, tck = -.02, cex.axis=1, at=c(1,2,3), lab=c('BE', "fBE", "HGD/EAC"), cex.lab=1)
axis(2, at=seq(0,90,10),tck = -.009, lab=seq(0,90,10), cex.axis=1, cex.lab=1)
title(ylab="Estimated BE onset time (age)", cex.lab=1)
plot(0:1,0:1,type="n",xlim=c(0.5,3.5),ylim=c(0,1),axes=FALSE,ann=FALSE)
vioplot(EACrisk_givenBE,EACrisk_givenBE_fBE,EACrisk_givenBE_BEwHGD,names=c("BE","fBE","HGD/EAC"),border= "aquamarine3",col="aquamarine3",add=T)
title(ylab="Associated EAC risk (probability)", cex.lab=1)
axis(1, tck = -.02, cex.axis=1, at=c(1,2,3), lab=c('BE', "fBE", "HGD/EAC"), cex.lab=1)
axis(2, at=seq(0,1,.1),tck = -.009, lab=seq(0,1,.10), cex.axis=1, cex.lab=1)



## EAC Risk until age T_final for given tau BE onset age for methylation data
## 4-stage model WITH BE CONVERSION 
phidot = function(s, phi, parms){
  phidot=numeric(6)
  with(as.list(parms),{
  tau = t - s
  RR=5
  GERDpr=1-exp(-gerd1*min(gerd3,tau)-gerd2*max(0,tau-gerd3))
  nu = nu0*((1-GERDpr)+RR*GERDpr)
  phidot[1]=betam -(alpham+betam+rho)*phi[1]+alpham*(phi[1]**2)
  phidot[2]=2*alpham*phi[1]*phi[2]- (alpham+betam+rho)*phi[2]
  phidot[3]=betap+ mu2*phi[1]*phi[3]-(alphap+betap+mu2)*phi[3]+alphap*(phi[3]**2)
  phidot[4]=2*alphap*phi[3]*phi[4]+mu2*(phi[4]*phi[1]+phi[3]*phi[2]) - (alphap+betap+mu2)*phi[4]
  phidot[5]=mu1*phi[5]*(phi[3]-1)
  phidot[6]=mu1*(phi[6]*(phi[3]-1)+phi[5]*phi[4])
  phidot[7]=mu0*phi[7]*(phi[5]-1)
  phidot[8]=mu0*(phi[8]*(phi[5]-1)+phi[7]*phi[6])
  ## A-D transition
  phidot[9]=nu*(phi[7]-phi[9])
  phidot[10]=nu*(phi[8]-phi[10])
  list(c(phidot))
  })
}


source('am_parameter_list.R')
source('af_parameter_list.R')
s_MSCE <- function(parms,t){
  ages<- t
  besurv.ode = NULL
  for (a in ages) {
      # initial condition
      parms['t'] = a
      phi0 = c(phi1=1, phi2=-allmales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0,phi9=1,phi10=0)
      times = c(0,a)
      ode = as.data.frame(lsoda(phi0,times,func=phidot,parms,atol=1.e-9,rtol=1.e-9))
      besurv.ode = c(besurv.ode,ode$phi7[2])
   }
  return(besurv.ode)
}  

maleparms = c(nu=allmales$nu, alphap=allmales$alphaP,betap=allmales$betaP,alpham=allmales$alphaM,betam=allmales$betaM,mu0= X*allmales$mu0,mu1=allmales$mu1,mu2=allmales$mu2,rho=allmales$rho,t=0)
femaleparms = c(nu=allfemales$nu, alphap=allfemales$alphaP,betap=allfemales$betaP,alpham=allfemales$alphaM,betam=allfemales$betaM,mu0= X*allfemales$mu0,mu1=allfemales$mu1,mu2=allfemales$mu2,rho=allfemales$rho,t=0)

T_final=88
tau_ages= 1:T_final

#tau_ages=apply(outgen$mcmc,2,mean)[1:30]
tau_ages=sum_matrixgen_ccf[3,]
tau_ageslow=sum_matrixgen_ccf[5,]
tau_ageshigh=sum_matrixgen_ccf[1,]
#tau_ages=apply(outgen[1:(burnend-burnin),id.agegen],2,mean)
risk_test=risk_testlow=risk_testhigh=rep(0,length(tau_ages))
EACrisk_givenBE = rep(0,length(tau_ages))
EACrisk_givenBElow = rep(length(tau_ages))
EACrisk_givenBEhigh = rep(length(tau_ages))
red_id= which(pd$gender=="Male")
blue_id=which(pd$gender=="Female")
for (i in red_id){
  EACrisk_givenBE[i] = 1-s_MSCE(maleparms,T_final-tau_ages[i])
  EACrisk_givenBElow[i]= 1-s_MSCE(maleparms,T_final-tau_ageslow[i])
  EACrisk_givenBEhigh[i]=1-s_MSCE(maleparms,T_final-tau_ageshigh[i])
  risk_test[i] = (s_MSCE(maleparms, age[i]-tau_ages[i])-s_MSCE(maleparms, T_final-tau_ages[i]))/(s_MSCE(maleparms,age[i]-tau_ages[i]))
  risk_testlow[i] = (s_MSCE(maleparms, age[i]-tau_ageslow[i])-s_MSCE(maleparms, T_final-tau_ageslow[i]))/(s_MSCE(maleparms,age[i]-tau_ageslow[i]))
  risk_testhigh[i] = (s_MSCE(maleparms, age[i]-tau_ageshigh[i])-s_MSCE(maleparms, T_final-tau_ageshigh[i]))/(s_MSCE(maleparms,age[i]-tau_ageshigh[i]))

}
for (i in blue_id){
  EACrisk_givenBE[i] = 1-s_MSCE(femaleparms,T_final-tau_ages[i])
  EACrisk_givenBElow[i]= 1-s_MSCE(femaleparms,T_final-tau_ageslow[i])
  EACrisk_givenBEhigh[i]=1-s_MSCE(femaleparms,T_final-tau_ageshigh[i])
  risk_test[i] = (s_MSCE(femaleparms, age[i]-tau_ages[i])-s_MSCE(femaleparms, T_final-tau_ages[i]))/(s_MSCE(femaleparms,age[i]-tau_ages[i]))
  risk_testlow[i] = (s_MSCE(femaleparms, age[i]-tau_ageslow[i])-s_MSCE(femaleparms, T_final-tau_ageslow[i]))/(s_MSCE(femaleparms,age[i]-tau_ageslow[i]))
  risk_testhigh[i] = (s_MSCE(femaleparms, age[i]-tau_ageshigh[i])-s_MSCE(femaleparms, T_final-tau_ageshigh[i]))/(s_MSCE(femaleparms,age[i]-tau_ageshigh[i]))

}

#tau_ages= apply(out_fBE$mcmc,2,mean)[1:22]
tau_ages = sum_matrixFBE_ccf[3,]
tau_ageslow=sum_matrixFBE_ccf[5,]
tau_ageshigh=sum_matrixFBE_ccf[1,]
risk_testFBE=risk_testFBElow=risk_testFBEhigh=rep(0,length(tau_ages))
#tau_ages= apply(out_fBE[1:(burnend-burnin),id.agefBE],2,mean)
EACrisk_givenBE_fBE= rep(0,length(tau_ages))
EACrisk_givenBE_fBElow = rep(length(tau_ages))
EACrisk_givenBE_fBEhigh = rep(length(tau_ages))
red_id= which(Sex.FBE=="M")
blue_id=which(Sex.FBE=="F")
for (i in red_id){
  EACrisk_givenBE_fBE[i] = 1-s_MSCE(maleparms,T_final-tau_ages[i])
  EACrisk_givenBE_fBElow[i]= 1-s_MSCE(maleparms,T_final-tau_ageslow[i])
  EACrisk_givenBE_fBEhigh[i]=1-s_MSCE(maleparms,T_final-tau_ageshigh[i])
  risk_testFBE[i] = (s_MSCE(maleparms, Age.FBE[i]-tau_ages[i])-s_MSCE(maleparms, T_final-tau_ages[i]))/(s_MSCE(maleparms,Age.FBE[i]-tau_ages[i]))
  risk_testFBElow[i] = (s_MSCE(maleparms, Age.FBE[i]-tau_ageslow[i])-s_MSCE(maleparms, T_final-tau_ageslow[i]))/(s_MSCE(maleparms,Age.FBE[i]-tau_ageslow[i]))
  risk_testFBEhigh[i] = (s_MSCE(maleparms, Age.FBE[i]-tau_ageshigh[i])-s_MSCE(maleparms, T_final-tau_ageshigh[i]))/(s_MSCE(maleparms,Age.FBE[i]-tau_ageshigh[i]))
}
for (i in blue_id){
  EACrisk_givenBE_fBE[i] = 1-s_MSCE(femaleparms,T_final-tau_ages[i])
  EACrisk_givenBE_fBElow[i]= 1-s_MSCE(femaleparms,T_final-tau_ageslow[i])
  EACrisk_givenBE_fBEhigh[i]=1-s_MSCE(femaleparms,T_final-tau_ageshigh[i])
  risk_testFBE[i] = (s_MSCE(femaleparms, Age.FBE[i]-tau_ages[i])-s_MSCE(femaleparms, T_final-tau_ages[i]))/(s_MSCE(femaleparms,Age.FBE[i]-tau_ages[i]))
  risk_testFBElow[i] = (s_MSCE(femaleparms, Age.FBE[i]-tau_ageslow[i])-s_MSCE(femaleparms, T_final-tau_ageslow[i]))/(s_MSCE(femaleparms,Age.FBE[i]-tau_ageslow[i]))
  risk_testFBEhigh[i] = (s_MSCE(femaleparms, Age.FBE[i]-tau_ageshigh[i])-s_MSCE(femaleparms, T_final-tau_ageshigh[i]))/(s_MSCE(femaleparms,Age.FBE[i]-tau_ageshigh[i]))

}

###  ONLY BE VS FBE VIOPLOT
quartz("Quartz", width=6.9,height=3.4)
par(ps=10,mfrow= c(1,2),mar=c(2,2,.5,.5), mgp=c(1,.25,0))
plot(0:1,0:1,type="n",xlim=c(0.5,2.5),ylim=c(0,90),axes=FALSE,ann=FALSE)
vioplot(apply(outgen_ccf[1:(burnend-burnin),id.agegen],2,median),apply(outFBE_ccf[1:(burnend-burnin),id.agefBE],2,median), names=c('BE','FBE'), col="aquamarine3", border="aquamarine3",ylim=c(0,90),add=T)
boxplot(Age.BE,Age.FBE, add=T,border="darkgrey", axes=F)
title(ylab="Estimated BE onset time (age)", cex.lab=1)
axis(1, tck = -.02, cex.axis=1, at=c(1,2), lab=c('BE', "FBE"), cex.lab=1)
axis(2, at=seq(0,90,10),tck = -.009, lab=seq(0,90,10), cex.axis=1, cex.lab=1)
plot(0:1,0:1,type="n",xlim=c(0.5,2.5),ylim=c(0,1),axes=FALSE,ann=FALSE)
#vioplot(EACrisk_givenBE,EACrisk_givenBE_fBE,names=c("BE","FBE"),border= "aquamarine3",col="aquamarine3",add=T)
vioplot(risk_test,risk_testFBE,names=c("BE","FBE"),border= "aquamarine3",col="aquamarine3",add=T)
title(ylab="Associated EAC risk (probability)", cex.lab=1)
axis(1, tck = -.02, cex.axis=1, at=c(1,2), lab=c('BE', "FBE"), cex.lab=1)
axis(2, at=seq(0,1,.1),tck = -.009, lab=seq(0,1,.10), cex.axis=1, cex.lab=1)





write.table(paste(signif(EACrisk_givenBE,3), "(", signif(EACrisk_givenBElow,3),",",signif(EACrisk_givenBEhigh,3), ")"),file='BEriskgen_w95.txt', quote=F, row.names=F,col.names=F)
write.table(paste(signif(EACrisk_givenBE_fBE,3), "(", signif(EACrisk_givenBE_fBElow,3),",",signif(EACrisk_givenBE_fBEhigh,3), ")"),file='BEriskFBE_w95.txt', quote=F, row.names=F,col.names=F)

## risks given no EAC at time of biopsy a
write.table(paste(signif(risk_test,3), "(", signif(risk_testlow,3),",",signif(risk_testhigh,3), ")"),file='BEriskgen_w95.txt', quote=F, row.names=F,col.names=F)
write.table(paste(signif(risk_testFBE,3), "(", signif(risk_testFBElow,3),",",signif(risk_testFBEhigh,3), ")"),file='BEriskFBE_w95.txt', quote=F, row.names=F,col.names=F)



