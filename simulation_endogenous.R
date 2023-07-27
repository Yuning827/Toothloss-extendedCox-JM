setwd("/Users/f.x.y/Library/Mobile Documents/com~apple~CloudDocs/UofT/Toothloss-extendedCox-JM")


args <- (commandArgs(TRUE))
print(args)

if (length(args)==0){
  print("No arguments supplied.")
  
  ## supply default values
  ncores <- 10
  nsim_core <- 150
  B = 200 # number of bootstrapping to get robust se
  
  visits <- 10
  t_mean_gap <- 1
  t_vary <- 0.3 # time_gap = t_mean_gap + runif(-t_vary, t_vary)
  max_time <- (t_mean_gap+t_vary) * (visits+1)
  
  N <- 100 # total_patients
  teeth_min <- 22
  teeth_max <- 28
 
  # fixed effect 
  delta = 1
  gamma_0 = -1
  gamma_1 = 0.6
  gamma_2 = 0.4
  alpha = 0.4
  beta_0 = 0.5
  beta_1 = 0.1
  beta_2 = 0.2
  beta_3 = 0.2
  
  # covariance matrix and mean of random effects
  b_covmat <- matrix(c(0.6, -0.2, -0.2, 0.1), 2, 2)
  b_means <- rep(0, 2)
  residual_sd <- 0.5
  
  # frailty variance 
  frailty_alpha <- 0.8
  
  # x1 is a binary covariate, x2 is a normal covariate
  x1_p <- 0.5
  x2_sigma2 <- 1
} else{
  for (i in (1:length(args))) eval(parse(text=args[[i]]))
}

#### load packages ####
## before load packages for the first time, install packages in terminal
library(survival)
library(Matrix)
library(lme4)
library(nlme)
library(splines)
library(statmod)
library(JointModel)
library(MASS)
library(JM)
library(Rcpp)
library(here)
library(tidyverse)
library(ggridges)
library(rstanarm)
library(simsurv)
library(parallel)
library(coxme)
options(warn=-1)

simtotal <- ncores * nsim_core

#task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
task_id <- 1
simgroup_min <- (task_id-1) * simtotal 
simgroup_max <- task_id * simtotal

simulation.func <- function(x){
  from <- seq(1, simtotal, by = nsim_core)[x]
  to <- seq(nsim_core, simtotal, by = nsim_core)[x]
  basecox_coef_result <- matrix(NA, nsim_core, 3+1)
  basecox_sd_result <- matrix(NA, nsim_core, 3+1)
  coxme_coef_result <- matrix(NA, nsim_core, 3+1)
  coxme_sd_result <- matrix(NA, nsim_core, 3+1)
  td_coef_result <- matrix(NA, nsim_core, 3+1)
  td_sd_result <- matrix(NA, nsim_core, 3+1)
  jm_coef_result <- matrix(NA, nsim_core, 13+1)
  jm_sd_result <- matrix(NA, nsim_core, 13+1)
  
  for (s in 1:nsim_core){
    sim <- from + s-1 + simgroup_min
    simseed <- from + s-1 + simgroup_min
    
    set.seed(simseed)
    out <- tryCatch({
      
      #### generate data #### 
      M <- ceiling(runif(N, teeth_min, teeth_max)) # number of teeth per subject
      MN <- sum(M) # total number of teeth in the data
      i <- rep(1:N, M)
      j <- unlist(lapply(1:N, FUN = function(x) 1:M[x]))
      time_gap_list <- lapply(1:N, FUN = function(x) 
        t_mean_gap + runif(visits, -t_vary, t_vary))
      time_gap <- data.frame(matrix(unlist(time_gap_list), visits, N))
      time <- rbind(rep(0, N), cumsum(time_gap))
      time_ij_list <- lapply(1:N, FUN = function(x) rep(time[,x], M[x]))
      time_ij <- data.frame(matrix(unlist(time_ij_list), visits+1, MN))
      
      # weibull hazard function
      haz <- function(t, x, betas, ...) {
        betas[["delta"]] * (t ^ (betas[["delta"]] - 1)) *  betas[["frailty_i"]] * exp(
          betas[["gamma_0"]] +
            betas[["gamma_1"]] * x[["x1"]] + 
            betas[["gamma_2"]] * x[["x2"]] + 
            betas[["alpha"]] * (
              betas[["beta_0i"]] + betas[["beta_1i"]] * t +
                betas[["beta_2"]]  * x[["x1"]] +
                betas[["beta_3"]]  * x[["x2"]]))
      }
      
      # a data frame of parameters
      betas <- data.frame(delta = rep(delta, MN), gamma_0 = rep(gamma_0, MN),
                          gamma_1 = rep(gamma_1, MN), gamma_2 = rep(gamma_2, MN), 
                          alpha = rep(alpha, MN),
                          beta_0 = rep(beta_0, MN), beta_1 = rep(beta_1, MN), 
                          beta_2 = rep(beta_2, MN), beta_3 = rep(beta_3, MN))
      
      # covariance matrix and mean of random effects
      b <- MASS::mvrnorm(n = MN, mu = b_means, Sigma = b_covmat)
      # subject level frailty 
      frailty_e <- rexp(N, 1)
      frailty_theta <- pi * runif(N, 0, 1)
      frailty_a <- (sin(frailty_alpha*frailty_theta)/sin(frailty_theta))^(1/(1-frailty_alpha)) * 
        (sin((1-frailty_alpha)*frailty_theta)/sin(frailty_alpha*frailty_theta))
   
      frailty <- (frailty_a/frailty_e)^((1-frailty_alpha)/frailty_alpha)
      
      frailty_i <- rep(frailty, M)
      # generate the beta_0i, beta_1i, frailty_i
      betas$beta_0i <- betas$beta_0 + b[, 1]
      betas$beta_1i <- betas$beta_1 + b[, 2]
      betas$frailty_i <- frailty_i
      # a data frame of covariates
      x1 <- stats::rbinom(N, 1, x1_p)
      x1 <- rep(x1, M)
      x2 <- stats::rnorm(N, 0, x2_sigma2)
      x2 <- rep(x2, M)
      covdat <- data.frame(x1,x2)

      
      # generate the tooth-level survival time
      times <- simsurv(hazard = haz, x = covdat, betas = betas, maxt=max_time)
      
      # generate the random error for biomarker
      epsilon_ijk <- matrix(rnorm((visits+1)*MN, 0, residual_sd), (visits+1), MN)
      # generate a data frame of biomarker at each visit
      y_ij <-  lapply(1:(visits+1), FUN = function(x) 
        betas[["beta_0i"]] + betas[["beta_1i"]] * time_ij[x,] +
          betas[["beta_2"]]  * covdat[["x1"]] +
          betas[["beta_3"]]  * covdat[["x2"]] + epsilon_ijk[x,])
      
      y.dat <- as.data.frame(matrix(unlist(y_ij), MN,  visits+1))
      colnames(y.dat) <- paste0("y_ij" ,0:visits)
      
      # create the wide format data, assume same observation time for all patients
      y.dat_wide <- data.frame(i, j, covdat, y.dat, eventtime = times$eventtime)
      
      # long format data
      sim_data_long <- y.dat_wide %>%
        pivot_longer(
          cols=starts_with("y_ij"),
          names_to = "visit",
          values_to = "Y"
        ) %>% 
        mutate(obstime = unlist(time_ij_list),
               visit = as.integer(substr(visit, 5, nchar(visit))),
               status = case_when(eventtime >= obstime ~ 0,
                                  eventtime < obstime ~ 1),
               ij = i*100+j) %>%
        group_by(ij) %>% 
        mutate(nrow_to_be_saved = case_when(visits - sum(status) < visits ~ visits - sum(status) + 2,
                                            visits - sum(status) == visits ~ visits+1)) %>% 
        filter(visit < nrow_to_be_saved) %>%
        mutate(visits_j = max(visit)) %>%
        #filter(!visits_j <= 1) %>%
        mutate_at(vars(Y), ~replace(., status == 1, NA))
      
      ## survival data
      sim_data_surv <- sim_data_long %>%
        dplyr::group_by(i, j) %>%
        mutate(baseY = Y[1L]) %>%
        arrange(obstime) %>%
        slice(n()) %>%
        ungroup() 
      
      ## data for extended Cox model 
      sim_tdsurv <- tmerge(sim_data_surv, sim_data_surv, id = ij, endpt = event(eventtime, status) )
      sim_tdsurv <- tmerge(sim_tdsurv, sim_data_long, id = ij,
                           Y = tdc(obstime, Y),
                           x1 = tdc(obstime, x1),
                           x2 = tdc(obstime, x2))
      
      #### model fitting #### 
      
      # base Cox model
      basecox <- coxph(Surv(obstime, status) ~  baseY + x1 + x2,
                       data = sim_data_surv, cluster = i)
      
      # base cox random effects
      coxmemod <- coxme(Surv(obstime, status) ~  baseY + x1 + x2 + (1|i),
                        data = sim_data_surv)
      
      
      ### extended Cox model
      sim_tdcox <- coxph(Surv(tstart, tstop, endpt) ~ Y + x1 + x2,
                         data = sim_tdsurv, cluster = i)
      
      ### joint model
      fitjm <- function(bootlong, bootsurv){
        sim_jmlong <- lme(Y ~ obstime + x1 + x2,
                          random = ~ obstime | ij,
                          data = bootlong,
                          #                        control = ctrl.lme
                          na.action = na.omit)
        #    ctrl.lme <- lmeControl(opt="optim", maxIter=100, msMaxIter=100)
        
        sim_jmcox <- coxph(Surv(eventtime, status) ~ x1 + x2,
                           cluster = i,
                           data = bootsurv,
                           x = TRUE)
        sim_jmtloss <- jointModel(sim_jmlong, sim_jmcox,
                                  timeVar = "obstime", 
                                  method = "weibull-PH-GH")
        return(sim_jmtloss)
      }
      sim_jmtloss <- fitjm(sim_data_long, sim_data_surv)
      jm_var <- solve(sim_jmtloss$Hessian)
      #### model summary ####
      basecox_coef <- basecox$coefficients
      coxme_coef <- coxmemod$coefficients
      td_coef <- sim_tdcox$coefficients
      jm_coef <-  c(sim_jmtloss$coefficients$betas, sim_jmtloss$coefficients$sigma,
                    sim_jmtloss$coefficients$gammas, sim_jmtloss$coefficients$alpha,
                    sim_jmtloss$coefficients$sigma.t, 
                    sim_jmtloss$coefficients$D[c(1,2,4)])
      
      basecox_sd <- sqrt(diag(basecox$var))
      coxme_sd <- tail(sqrt(diag(coxmemod$var)),3)
      td_sd <- sqrt(diag(sim_tdcox$var))
      jm_sd <- sqrt(diag(jm_var))
    }, error=function(e) e)
    #### result ####
    if (inherits(out, "error")) {
      basecox_coef_result[s,] <- c(sim, rep(NA, 3))
      basecox_sd_result[s,] <- c(sim, rep(NA, 3))
      coxme_coef_result[s,] <- c(sim, rep(NA, 3))
      coxme_sd_result[s,] <- c(sim, rep(NA, 3))
      td_coef_result[s,] <- c(sim, rep(NA, 3))
      td_sd_result[s, ] <- c(sim, rep(NA, 3))
      jm_coef_result[s, ] <- c(sim, rep(NA, 13))
      jm_sd_result[s, ] <- c(sim, rep(NA, 13))
    } else {
      basecox_coef_result[s,] <- c(sim, as.vector(basecox_coef))
      basecox_sd_result[s,] <- c(sim, as.vector(basecox_sd))
      coxme_coef_result[s,] <- c(sim, as.vector(coxme_coef))
      coxme_sd_result[s,] <- c(sim, as.vector(coxme_sd))
      td_coef_result[s,] <- c(sim, as.vector(td_coef))
      td_sd_result[s, ] <- c(sim, as.vector(td_sd))
      jm_coef_result[s, ] <- c(sim, as.vector(jm_coef))
      jm_sd_result[s, ] <- c(sim, as.vector(jm_sd))
    }
  }
  return(list(basecox_coef_result, basecox_sd_result,
              coxme_coef_result, coxme_sd_result,
              td_coef_result, td_sd_result,
              jm_coef_result, jm_sd_result))
  
}

system.time(simulation.out <- mclapply(1:ncores, simulation.func, 
                                       mc.cores = ncores, mc.silent = F))

simout.df <- as.data.frame(do.call(rbind, simulation.out))


basecox_coef_sim_result <- as.data.frame(do.call(rbind, simout.df$V1))
basecox_sd_sim_result <- as.data.frame(do.call(rbind, simout.df$V2))
coxme_coef_sim_result <- as.data.frame(do.call(rbind, simout.df$V3))
coxme_sd_sim_result <- as.data.frame(do.call(rbind, simout.df$V4))
td_coef_sim_result <- as.data.frame(do.call(rbind, simout.df$V5))
td_sd_sim_result <- as.data.frame(do.call(rbind, simout.df$V6))
jm_coef_sim_result <- as.data.frame(do.call(rbind, simout.df$V7))
jm_sd_sim_result <- as.data.frame(do.call(rbind, simout.df$V8))

colnames(basecox_coef_sim_result) <- c("sim", "alpha", "gamma_x1", "gamma_x2")
colnames(basecox_sd_sim_result) <- c("sim", "alpha", "gamma_x1", "gamma_x2")
colnames(coxme_coef_sim_result) <- c("sim", "alpha", "gamma_x1", "gamma_x2")
colnames(coxme_sd_sim_result) <- c("sim", "alpha", "gamma_x1", "gamma_x2")
colnames(td_coef_sim_result) <- c("sim", "alpha", "gamma_x1", "gamma_x2")
colnames(td_sd_sim_result) <- c("sim", "alpha", "gamma_x1", "gamma_x2")
colnames(jm_coef_sim_result) <- c("sim", "beta0", "beta_obstime", "beta_x1", "beta_x2",
                                  "sigma", "gamma0", "gamma_x1", "gamma_x2",
                                  "alpha", "sigma.t", 
                                  "re_sigma11", "re_sigma12", "re_sigma22")
colnames(jm_sd_sim_result) <- c("sim", "beta0", "beta_obstime", "beta_x1", "beta_x2",
                                "sigma", "gamma0", "gamma_x1", "gamma_x2",
                                "alpha", "sigma.t", 
                                "re_sigma11", "re_sigma12", "re_sigma22")

true_parameter <- c(beta_0, beta_1, beta_2, beta_3, residual_sd, 
                    gamma_0 * frailty_alpha, gamma_1 * frailty_alpha, gamma_2 * frailty_alpha, 
                    alpha * frailty_alpha, delta * frailty_alpha,
                    b_covmat[1,1], b_covmat[1,2], b_covmat[2,2])
result_jm_mean <- colMeans(na.omit(jm_coef_sim_result))[-1]
bias_jm <- (result_jm_mean - true_parameter)
relative_bias_jm <- (result_jm_mean - true_parameter)/true_parameter
mse_jm <- colMeans((sweep(na.omit(jm_coef_sim_result)[-1], 2, true_parameter))^2)
empirical_se_jm <- apply(na.omit(jm_coef_sim_result)[-1], 2, sd)


true_parameter_basecox <- c(alpha * frailty_alpha, gamma_1 * frailty_alpha, gamma_2 * frailty_alpha)
result_basecox_mean <- colMeans(na.omit(basecox_coef_sim_result))[-1]
relative_bias_basecox <- (result_basecox_mean - true_parameter_basecox)/true_parameter_basecox
mse_basecox <- colMeans((sweep(na.omit(basecox_coef_sim_result)[-1], 2, true_parameter_basecox))^2)
empirical_se_basecox <- apply(na.omit(basecox_coef_sim_result)[-1], 2, sd)


true_parameter_coxme <- c(alpha * frailty_alpha, gamma_1 * frailty_alpha, gamma_2 * frailty_alpha)
result_coxme_mean <- colMeans(na.omit(coxme_coef_sim_result))[-1]
relative_bias_coxme <- (result_coxme_mean - true_parameter_coxme)/true_parameter_coxme
mse_coxme <- colMeans((sweep(na.omit(coxme_coef_sim_result)[-1], 2, true_parameter_coxme))^2)
empirical_se_coxme <- apply(na.omit(coxme_coef_sim_result)[-1], 2, sd)


true_parameter_td <- c(alpha * frailty_alpha, gamma_1 * frailty_alpha, gamma_2 * frailty_alpha)
result_td_mean <- colMeans(na.omit(td_coef_sim_result))[-1]
relative_bias_td <- (result_td_mean - true_parameter_td)/result_td_mean
mse_td <- colMeans((sweep(na.omit(td_coef_sim_result)[-1], 2, true_parameter_td))^2)
empirical_se_td <- apply(na.omit(td_coef_sim_result)[-1], 2, sd)



final_jm_result <- round(data.frame(true_parameter, result_jm_mean, relative_bias_jm, 
                                    mse_jm, empirical_se_jm),4)

final_basecox_result <- round(data.frame(true_parameter_basecox, result_basecox_mean, relative_bias_basecox, 
                                         mse_basecox, empirical_se_basecox),4)

final_coxme_result <- round(data.frame(true_parameter_coxme, result_coxme_mean, relative_bias_coxme, 
                                       mse_coxme, empirical_se_coxme),4)

final_td_result <- round(data.frame(true_parameter_td, result_td_mean, relative_bias_td, 
                                    mse_td, empirical_se_td),4)

colnames(final_jm_result) <- c("JM_true_parameter", "mean_est", "relative_bias", "mse", "empirical_se")
colnames(final_basecox_result) <- c("Cox_true_parameter", "mean_est", "relative_bias", "mse", "empirical_se")
colnames(final_coxme_result) <- c("Coxme_true_parameter", "mean_est", "relative_bias", "mse", "empirical_se")
colnames(final_td_result) <- c("TD_true_parameter", "mean_est", "relative_bias", "mse", "empirical_se")

final_jm_result_selected <- final_jm_result[c(9,7,8),]


filename <- paste0("JM.compare.result.sim_",simtotal,"_N_", N, "_frailty_", frailty_alpha, ".RData")
save(file=filename,
     list=c("final_jm_result", "final_basecox_result", "final_coxme_result", "final_td_result",
            "simout.df"))

load(file=filename)


