# BE_Clock_Model

These are the scripts and function necessary to implement the Bayesian BE clock model developed by Kit Curtius 
with contributions by Georg Luebeck.

Clock_CpG_selection.R
- contains script for statistical pipeline used to select BE Clock CpG markers from 450K methylation marker sites
- allows for loading any longitudinal or cross-sectional methylation data the user chooses
- builds a list of CpG names that will be used in the BE Clock model to estimate tissue age 
    and the corresponding drift rates for each CpG derived from regressions across serial samples 
    stored in data frame CpGs.D1.67 (because 67 CpGs were selected for through our pipeline and data sets)
    
Bayesian_model.R
- sets up matrices to contain Markov Chain Monte Carlo (MCMC) results for parameters of interest (namely BE onset ages s)
- contains a loop over individuals in specified BE patient data set in which the Gibbs sampler
  function is called ("Gibbs_sampler.R") to perform MCMC and infer patient-specific parameteres for specified # of runs
- produces boxplots for posterior estimates of BE onset ages for each individual across a data set
- performs Bayes factor testing comparing sporadic BE patients with FBE patients to test if FBE has lived a larger percentage of 
  life with BE than the sporadic BE patients
- computes the EAC patient-specific lifetime risk (i.e., until some user input age) for each patient 
  given a posterior estimate drawn by MCMC (e.g., the median) assuming the Multistage Clonal Expansion 
  for EAC Model given BE onset (see Curtius et al., 2015 PLoS CompBio).  
- produces a violin plot of the medians BE onset age estimate for both sporadic BE cases and FBE cases and also
  a violin plot of the associated EAC lifetime risk for the onset ages
  
Gibbs_sampler.R
- takes user input on starting parameter values
- performs Gibbs sampling using full conditionals to estimate patient-specific BE onset ages, BE clock CpG drift rates, and 
  standard deviation of methylation M-value measurement 
- can also plot the MCMC chains 

MSCE_EAC_params.R
- contains the male and female, all race cellular kinetic parameters that were derived by fitting the MSCE-EAC
  model hazard function to Surveillance, Epidemiology, and End Results (SEER) level incidence data. (See Kong et al., 2014 CEBP)
- may be replaced for user defined parameters for gender-specific mutation rates, dysplastic growth rates, etc.


    
