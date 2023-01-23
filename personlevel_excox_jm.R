source("create_personleveldata.R")


### extended Cox model
td.tl <- coxph(Surv(tstart, tstop, TL) ~ pctpocket5mm + pctabl40 + pctmobil05mm + 
                 baseage + basebmi + basesmoke +college + basenumteeth, data = tddat)
summary(td.tl)  


### joint model

jm.tl <- stan_jm(formulaLong = list(pctpocket5mm ~ year + (1 + year| id), # adding baseline covariates did not change the results
                                    pctmobil05mm ~ year + (1 + year| id),
                                    pctabl40 ~ year + (1 + year| id)),
                 dataLong = tlLong,
                 formulaEvent = Surv(years, TL) ~ baseage + basebmi + basesmoke +college + basenumteeth,
                 dataEvent = tlSurv,
                 time_var = "year",
                 # this next line is only to keep the example small in size!
                 chains = 1, cores = 1, seed = 12345, iter = 1000)
print(jm.tl)
summary(jm.tl, probs = c(0.025, 0.975))