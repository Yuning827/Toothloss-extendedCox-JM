# Compare joint model compared with extended Cox regression model

R codes for the paper "Modeling time-varying risk factors of tooth loss: results from joint model compared to extended Cox regression model". We provided example datasets to illustrate the models in the paper instead of the original dataset. The results from example datasets are different to the results in the paper. 

## vadls_jm_personlevel.csv 

Example dataset for person level analysis

## vadls_jm_toothlevel.csv

Example dataset for tooth level analysis

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
