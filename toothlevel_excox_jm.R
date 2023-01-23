source("create_toothleveldata.R")

# base Cox model
basecox <- coxph(Surv(time, tloss) ~  basepock + baseage + basebmi + basesmoke +college,
                 data = vasurv, cluster = id)
summary(basecox)


# extended Cox model
tdcox <- coxph(Surv(tstart, tstop, endpt) ~ maxpock + baseage + basebmi + basesmoke +college,
               data = tdsurv, cluster = id)

summary(tdcox)



# joint model
# first, fit longitudinal submodel
jmlong <- lme(maxpock ~ obstime,
              random = ~ obstime | idtooth,
              data = valong,
              na.action = na.omit)
summary(jmlong)

# next, fit survival submodel
jmcox <- coxph(Surv(time, tloss) ~ baseage + basebmi + basesmoke +college,
               cluster = id,
               data = vasurv,
               x = TRUE)
summary(jmcox)

# finally, fit joint model (current value model)
jmtloss <- jointModel(jmlong, jmcox,
                      timeVar = "obstime",
                      method = "piecewise-PH-aGH")

summary(jmtloss)

