# Compare joint model with extended Cox regression model

R codes for the paper "Joint models to estimate risk of tooth loss from time-varying covariates". We provided example data sets that are simulated to illustrate the models in the paper instead of the original data set. The results from example data sets are different to the results in the paper. 

## vadls_jm_personlevel.csv 

Example data set for person level analysis

## vadls_jm_toothlevel.csv

Example data set for tooth level analysis

## personlevel_excox_jm.R

R code for running person level data analysis

The code includes: 

1. create the data for extended Cox model 

2. create the data for joint model 

3. run the extended Cox model

4. run the joint model 

## toothlevel_excox_jm.R

R code for running tooth level data analysis

The code includes: 

1. create the longitudinal and survival data

2. create the data for extended Cox model 

3. run the extended Cox model

4. run the joint model 

## figure_percentpdlong.R

R code for creating longitudinal trajectories of periodontal measurements of four patients

## simulation_endogenous.R

R code for running simulation study assuming endogenous time-varying covariate

The code includes:

1. generate simulated data set assuming a lme model for the endogenous covariate and a frailty model with positive stable distributed frailty term for the time-to-event outcome

2. run the simulation 1000 times

3. compare the following 4 models: a coxph model assuming that the endogenous covariate is constant, a coxme model assuming that the endogenous covariate is constant, an extended model, and a joint model

## simulation_exogenous.R

R code for running simulation study assuming exogenous time-varying covariate

The code includes:

1. generate simulated data set assuming a frailty model with positive stable distributed frailty term and an exogenous covariate for the time-to-event outcome

2. run the simulation 1000 times

3. compare the following 4 models: a coxph model assuming that the exogenous covariate is constant, a coxme model assuming that the exogenous covariate is constant, an extended model, and a joint model



