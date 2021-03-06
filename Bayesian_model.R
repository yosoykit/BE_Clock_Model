

#######################################################
##    This script runs the Bayesian BE Clock model   ##
##      on the cross-sectional patients (BE & FBE).  ##
##    It requires the function 'gibbs' for inference ##
###################################################################
##  SETUP for MCMC (see Main TexT and S1 Text) for expressions)  ##
###################################################################
source("Gibbs_sampler.R")
CpGs.D1.67 = load("CpGs.D1.67.rda", verbose=T)
CpGs.DV.67.CCF = load("CpGs.DV.67.CCF.rda", verbose=T)
##################################
## 1) ALL data use slopes and intercepts from 30 patient cross-sectional
##      population level drift derived from linear regressoin
###################################
source("mset.paper.rda",verbose=T)
pd_total = phenoData(mset)
# pheno data names for D2, first is BE second is matched SQ
source("D2_pd.rda", verbose=T)
# pheno data names for D3
source("D3_pd.rda", verbose=T)

D2_id = match(D2_pd$DNA.soln, pd_total$DNA.soln)
be_id = which(pd_total$Tissue_Type[D2_id]=="BE")
beM_D2 =getM(mset)[,D2_id[be_id]]
sqM_D2 =getM(mset)[,D2_id[-be_id]]

M = 67
Age.BE = as.integer(pd_total$Age_at_biopsy[D2_id[be_id]])
ids= na.omit(match(CpGs.DV.67.CCF$cpgn[1:M],dimnames(mset)[[1]]))
id = ids
slope.BE = intercept.BE = slope.SQ = intercept.SQ = numeric(length(id))
k=1
for(i in id) {
tmp = lm(beM_D2[i,] ~ Age.BE) 
intercept.BE[k] = tmp$coef[1]; slope.BE[k] = tmp$coef[2]
tmp = lm(sqM_D2[i,] ~ Age.BE)
intercept.SQ[k] = tmp$coef[1]; slope.SQ[k] = tmp$coef[2]
k=k+1
}

## for robustness study using matched vs. un-matched samples
#alpha_SQ_1=alpha_SQ
#b_SQ_1=b_SQ
#alpha_SQ=rep(0,67)
#b_SQ=rep(0,67)

################################################
## 2) group specific runs to be done separately###
################################################

## general population data (sporadic BE cases)
ages=Age.BE
age=ages
## mean and sd for normal prior on BE drift slopes derived from serial data sets (either D1 or DV)
mean_b_i = rep(mean(CpGs.DV.67.CCF$b.slope),67)
sd_b_i = rep(sd(CpGs.DV.67.CCF$b.slope),67)
N=length(age)
sum_matrixgen = matrix(rep(0,5*N),ncol=N)
beM_temp=beM_D2
patients = 1:N
## For robustness study considering matched vs. un-matched samples
#sum_matrixgen_blind = matrix(rep(0,5*N),ncol=N)


## FBE data 
FBE_id=which(pd_total$Tissue_Type=="FBE")
Age.FBE=as.integer(pd_total$Age_at_biopsy[FBE_id])
## missing 3555-FAM BE age = 44 M
Age.FBE[which(pd_total$DNA.soln[FBE_id]=="3555-FAM BE")]=44
pd_total$Gender[which(pd_total$DNA.soln=="3555-FAM BE")]="M"
fbeMset = getM(mset)[,FBE_id]
ages=Age.FBE
age=ages
mean_b_i = rep(mean(CpGs.DV.67.CCF$b.slope),67)
sd_b_i = rep(sd(CpGs.DV.67.CCF$b.slope),67)
N=length(age)
beM_temp=fbeMset
sum_matrixFBE = matrix(rep(0,5*N),ncol=N)
## Patient DNA.soln 4280 has < 0 predicted onset age with Prob ~ 1
patients = 1:22
patients=patients[-which(pd_total$DNA.soln[FBE_id]=="4280")]

##############################################################
## 3) ALL data MCMC LOOP OVER ALL INDIVIDUALS IN DATA SET ##
##      Records medians, quartiles, and 95% CI            ##
##      for BE onset ages inferred for each patient       ##
##############################################################
## Number of MCMC cycles to run for each patient
runs = 100000
M=length(id)
alpha_SQ = intercept.SQ
b_SQ = slope.SQ
sum_matrix = matrix(rep(0,5*N),ncol=N)
out_temp = matrix(rep(0,runs*N),ncol=N)
for (i in patients){
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
  ## start BE onset age at half of age-at-biopsy
  x$est[M+1] = a/2
  ## start slopes at mean values from regressions
  x$est[1:M] = mean_b_i

  ## estimates of standard deviations for dnorm, will be given gamma prior distribution
  x$est[M+2] = 1.2  #sigma_BEstart
  # For robustness analysis
  # x$est[M+2] = 1.2 *sqrt(2)
  print(i)
  ## Call function gibbs that runs gibbs sampler for this patient for runs # of cycles
  out = gibbs(x, beM_temp,L=runs,K=500,skip=10,plot=F)
  out_temp[,i] = out[,(M+1)]
  s = quantile(out[,(M+1)],probs=c(.025,.25,.5,.75,.975))
  sum_matrix[,i]=s 
}

##############################
##  4) group specific       ##
##############################
# general BE data
sum_matrixgen= sum_matrix
outgen= out_temp
onsetsgen = rep(0,N)
sum_matrixgen_round = signif(sum_matrixgen,4)
for (i in 1:N){ onsetsgen[i]=paste(sum_matrixgen_round[3,i],"(",sum_matrixgen_round[1,i],",",sum_matrixgen_round[5,i],")")}
#write.table(onsetsgen,"BEonsetsgen_w95_170316.txt",quote=F, row.names=F)
# blind data
#sum_matrixgen_blind = sum_matrix
#outgen_blind=out_temp

# fBE data
sum_matrixFBE= sum_matrix
outFBE = out_temp
onsetsFBE = rep(0,N)
sum_matrixFBE_round = signif(sum_matrixFBE,4)
for (i in 1:N){ onsetsFBE[i]=paste(sum_matrixFBE_round[3,i],"(",sum_matrixFBE_round[1,i],",",sum_matrixFBE_round[5,i],")")}
write.table(onsetsFBE,"BEonsetsFBE_w95_170316.txt",quote=F, row.names=F)
write.table(Age.FBE,"AgesFBE_170316.txt",quote=F, row.names=F)
write.table(pd_total$DNA.soln[FBE_id],"DNAsolnFBE_170316.txt",quote=F, row.names=F)

#######################
## PLOTTING BOXPLOTS ##
#######################
id.agegen=order(Age.BE)
id.agefBE=order(Age.FBE)

pd_sex_orderedgen=pd_total$Gender[D2_id[be_id[id.agegen]]]
pd_sex_orderedfBE = pd_total$Gender[FBE_id[id.agefBE]]


##########################
##   group specific     ##
##########################
## general BE data
N=30
ages=Age.BE
id.age = order(ages)
pd_sex_ordered=pd_sex_orderedgen 
sum_matrixtemp = sum_matrixgen
sum_matrixtemp=sum_matrixgen_blind
## fBE data
N=22
ages=as.integer(Age.FBE)
id.age = order(ages)
pd_sex_ordered=pd_sex_orderedfBE
sum_matrixtemp = sum_matrixFBE


### Boxplots for 2 independent cross-sectional data sets ###
#library(openintro)
quartz("Quartz", width=6.9,height=3.4)
par(ps=10,mar=c(2,2,.5,.5), mgp=c(1,.25,0))
temp=boxplot(sum_matrixtemp[,id.age],names=as.character(1:length(ages)))
bxp(temp,ylim=c(0,90),xlab="Patient ID",ylab="BE onset age",cex.lab=1,las=2,axes=F, outline=F, border=fadeColor("gray"),boxfill=fadeColor("gray"))
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
points(1:N,sum_matrixtemp[3,id.age],col="black", pch=18,cex=.85)
box()

##########################################################################################
##                    BAYES FACTOR TESTING: Bayesian hypothesis test to see if         ##
## FBE patients have lived longer with BE at time of diagnosis than sporadic BE cases  ##
##########################################################################################
N2= 30
N3= 22
burnend = runs
burnin = 0
sum1 = sum2 = sum3=rep(0,(burnend-burnin))
count=1
for ( k in 1:(burnend-burnin)){
  temp = sum((Age.BE- outgen[k,1:N2])/(Age.BE))
  sum1[count] = temp/N2
   temp = sum((Age.FBE - outFBE[k,(1:N3)])/(Age.FBE))
   #temp = sum((Age.FBE[-19] - outFBE[k,c(1:18,20:22)])/(Age.FBE[-19]))
   sum2[count] = temp/N3
  count=count+1
}


count=1
sum1_p = sum2_p = sum3_p=rep(0,(burnend-burnin))

for ( k in 1:(burnend-burnin)){
  prior_temp = runif(N2,0,Age.BE)
  temp = sum((Age.BE - prior_temp)/(Age.BE))
  sum1_p[count] = temp/N2
  prior_temp = runif(N3,0,Age.FBE)
  temp = sum((Age.FBE - prior_temp)/(Age.FBE))
  #temp = sum((Age.FBE[-19] - prior_temp)/(Age.FBE[-19]))
  sum2_p[count] = temp/N3
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

###############################################################################################
##              EAC Patient-specific Risk: predicted until age T_final for given 
##                tau BE onset ages inferred above for methylation data using a
## 4-stage multistage model (stochastic, multi-type branching process, see Main Text and S1 Fig) 
###############################################################################################
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
  ## Armitage-Doll transition
  phidot[9]=nu*(phi[7]-phi[9])
  phidot[10]=nu*(phi[8]-phi[10])
  list(c(phidot))
  })
}

library(deSolve)

source('MSCE_EAC_params.R')
## survival function s_MSCE
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

tau_ages=sum_matrixgen[3,]
tau_ageslow=sum_matrixgen[5,]
tau_ageshigh=sum_matrixgen[1,]
risk_test=risk_testlow=risk_testhigh=rep(0,length(tau_ages))
red_id= which(pd_total$Gender[D2_id[be_id]]=="Male")
blue_id=which(pd_total$Gender[D2_id[be_id]]=="Female")
for (i in red_id){
  #EACrisk_givenBE[i] = 1-s_MSCE(maleparms,T_final-tau_ages[i])
  #EACrisk_givenBElow[i]= 1-s_MSCE(maleparms,T_final-tau_ageslow[i])
  #EACrisk_givenBEhigh[i]=1-s_MSCE(maleparms,T_final-tau_ageshigh[i])
  risk_test[i] = (s_MSCE(maleparms, Age.BE[i]-tau_ages[i])-s_MSCE(maleparms, T_final-tau_ages[i]))/(s_MSCE(maleparms,Age.BE[i]-tau_ages[i]))
  #risk_testlow[i] = (s_MSCE(maleparms, age[i]-tau_ageslow[i])-s_MSCE(maleparms, T_final-tau_ageslow[i]))/(s_MSCE(maleparms,age[i]-tau_ageslow[i]))
  #risk_testhigh[i] = (s_MSCE(maleparms, age[i]-tau_ageshigh[i])-s_MSCE(maleparms, T_final-tau_ageshigh[i]))/(s_MSCE(maleparms,age[i]-tau_ageshigh[i]))
}
for (i in blue_id){
  #EACrisk_givenBE[i] = 1-s_MSCE(femaleparms,T_final-tau_ages[i])
  #EACrisk_givenBElow[i]= 1-s_MSCE(femaleparms,T_final-tau_ageslow[i])
  #EACrisk_givenBEhigh[i]=1-s_MSCE(femaleparms,T_final-tau_ageshigh[i])
  risk_test[i] = (s_MSCE(femaleparms, Age.BE[i]-tau_ages[i])-s_MSCE(femaleparms, T_final-tau_ages[i]))/(s_MSCE(femaleparms,Age.BE[i]-tau_ages[i]))
  #risk_testlow[i] = (s_MSCE(femaleparms, age[i]-tau_ageslow[i])-s_MSCE(femaleparms, T_final-tau_ageslow[i]))/(s_MSCE(femaleparms,age[i]-tau_ageslow[i]))
  #risk_testhigh[i] = (s_MSCE(femaleparms, age[i]-tau_ageshigh[i])-s_MSCE(femaleparms, T_final-tau_ageshigh[i]))/(s_MSCE(femaleparms,age[i]-tau_ageshigh[i]))
}

tau_ages = sum_matrixFBE[3,]
tau_ageslow=sum_matrixFBE[5,]
tau_ageshigh=sum_matrixFBE[1,]
risk_testFBE=risk_testFBElow=risk_testFBEhigh=rep(0,length(tau_ages))
red_id= c(which(pd_total$Gender[FBE_id]=="M"),which(pd_total$Gender[FBE_id]=="Male"))
blue_id=which(pd_total$Gender[FBE_id]=="F")
for (i in red_id){
  #EACrisk_givenBE_fBE[i] = 1-s_MSCE(maleparms,T_final-tau_ages[i])
  #EACrisk_givenBE_fBElow[i]= 1-s_MSCE(maleparms,T_final-tau_ageslow[i])
  #EACrisk_givenBE_fBEhigh[i]=1-s_MSCE(maleparms,T_final-tau_ageshigh[i])
  risk_testFBE[i] = (s_MSCE(maleparms, Age.FBE[i]-tau_ages[i])-s_MSCE(maleparms, T_final-tau_ages[i]))/(s_MSCE(maleparms,Age.FBE[i]-tau_ages[i]))
  #risk_testFBElow[i] = (s_MSCE(maleparms, Age.FBE[i]-tau_ageslow[i])-s_MSCE(maleparms, T_final-tau_ageslow[i]))/(s_MSCE(maleparms,Age.FBE[i]-tau_ageslow[i]))
  #risk_testFBEhigh[i] = (s_MSCE(maleparms, Age.FBE[i]-tau_ageshigh[i])-s_MSCE(maleparms, T_final-tau_ageshigh[i]))/(s_MSCE(maleparms,Age.FBE[i]-tau_ageshigh[i]))
}
for (i in blue_id){
  #EACrisk_givenBE_fBE[i] = 1-s_MSCE(femaleparms,T_final-tau_ages[i])
  #EACrisk_givenBE_fBElow[i]= 1-s_MSCE(femaleparms,T_final-tau_ageslow[i])
  #EACrisk_givenBE_fBEhigh[i]=1-s_MSCE(femaleparms,T_final-tau_ageshigh[i])
  risk_testFBE[i] = (s_MSCE(femaleparms, Age.FBE[i]-tau_ages[i])-s_MSCE(femaleparms, T_final-tau_ages[i]))/(s_MSCE(femaleparms,Age.FBE[i]-tau_ages[i]))
  #risk_testFBElow[i] = (s_MSCE(femaleparms, Age.FBE[i]-tau_ageslow[i])-s_MSCE(femaleparms, T_final-tau_ageslow[i]))/(s_MSCE(femaleparms,Age.FBE[i]-tau_ageslow[i]))
  #risk_testFBEhigh[i] = (s_MSCE(femaleparms, Age.FBE[i]-tau_ageshigh[i])-s_MSCE(femaleparms, T_final-tau_ageshigh[i]))/(s_MSCE(femaleparms,Age.FBE[i]-tau_ageshigh[i]))

}

###  ONLY BE VS FBE VIOPLOT
# library(vioplot)
quartz("Quartz", width=6.9,height=3.4)
par(ps=10,mfrow= c(1,2),mar=c(2,2,.5,.5), mgp=c(1,.25,0))
plot(0:1,0:1,type="n",xlim=c(0.5,2.5),ylim=c(0,90),axes=FALSE,ann=FALSE)
vioplot(apply(outgen[1:(burnend-burnin),id.agegen],2,median),apply(outFBE[1:(burnend-burnin),id.agefBE],2,median), names=c('BE','FBE'), col="aquamarine3", border="aquamarine3",ylim=c(0,90),add=T)
boxplot(Age.BE,Age.FBE, add=T,border="darkgrey", axes=F)
title(ylab="Estimated BE onset time (age)", cex.lab=1)
axis(1, tck = -.02, cex.axis=1, at=c(1,2), lab=c('BE', "FBE"), cex.lab=1)
axis(2, at=seq(0,90,10),tck = -.009, lab=seq(0,90,10), cex.axis=1, cex.lab=1)
plot(0:1,0:1,type="n",xlim=c(0.5,2.5),ylim=c(0,1),axes=FALSE,ann=FALSE)
vioplot(risk_test,risk_testFBE,names=c("BE","FBE"),border= "aquamarine3",col="aquamarine3",add=T)
title(ylab="Associated EAC risk (probability)", cex.lab=1)
axis(1, tck = -.02, cex.axis=1, at=c(1,2), lab=c('BE', "FBE"), cex.lab=1)
axis(2, at=seq(0,1,.1),tck = -.009, lab=seq(0,1,.10), cex.axis=1, cex.lab=1)


