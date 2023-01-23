#### Personal level data with some noise ####

tldat <- read.csv("vadls_jm2.csv", header = TRUE)

tldat <- tldat %>%
  mutate(college = ifelse(edu == 2, 1, 0)) %>%
  dplyr::group_by(id) %>%
  mutate(basepctpocket5mm = pctpocket5mm[1L],
         basepctcal5mm = pctcal5mm[1L],
         basepctmobil05mm = pctmobil05mm[1L],
         basepctabl40 = pctabl40[1L],
         basemeanpocket5mm = meanpocket5mm[1L],
         basemeancal5mm = meancal5mm[1L],
         basemeanmobil05mm = meanmobil05mm[1L],
         basemeanabl40 = meanabl40[1L],
         basesmoke = smoking[1L],
         basemets = atp_mets[1L],
         basetriglyc = triglyc[1L],
         basewtkg = wtkg[1L],
         basebmi = bmi[1L],
         baseage = age[1L],
         basesysbp = sysbp[1L],
         basediabp = diabp[1L],
         basefgluc = fgluc[1L],
         basehdl = hdl[1L]) %>%
  ungroup() %>%
  dplyr::select(id, year,  pctpocket5mm, pctmobil05mm, pctabl40, TL,
                baseage, college, basesmoke, basebmi, basenumteeth) %>% 
  group_by(id) %>%
  mutate(tstart = year, 
    tstop = dplyr::lead(tstart),
         TL = dplyr::lead(TL),
         diff = tstop - tstart) %>%
  dplyr::filter(row_number()!=n()|n()==1) %>% 
  fill(pctpocket5mm, pctmobil05mm, pctabl40, baseage, college, basesmoke, basebmi, .direction = "downup") %>%
  mutate(pctpocket5mm = coalesce(pctpocket5mm, 0),
         pctmobil05mm = coalesce(pctmobil05mm, 0),
         pctabl40 = coalesce(pctabl40, 0))

tldat <- tldat[complete.cases(tldat), ] 

tldat <- tldat %>% 
  group_by(id) %>%
  mutate(visits = n())


tldat_summary <- tldat %>% 
  group_by(id) %>% 
  summarise(baseage = mean(baseage), 
            basebmi = mean(basebmi))
tldat_summary[is.na(tldat_summary)] <- 0



diff_summary <- tldat %>% 
  group_by(id) %>% 
  summarise(diff_mean = mean(diff), 
            diff_sd = sd(diff))
diff_summary[is.na(diff_summary)] <- 0 


pocket_summary <- tldat %>% 
  filter(pctpocket5mm!=0) %>% 
  group_by(id) %>% 
  summarise(pocket_mean = mean(pctpocket5mm), 
            pocket_sd = sd(pctpocket5mm))
pocket_summary[is.na(pocket_summary)] <- 0

mobil_summary <- tldat %>% 
  filter(pctmobil05mm!=0) %>% 
  group_by(id) %>% 
  summarise(mobil_mean = mean(pctmobil05mm), 
            mobil_sd = sd(pctmobil05mm))
mobil_summary[is.na(mobil_summary)] <- 0

abl_summary <- tldat %>%   
  filter(pctabl40!=0) %>% 
  group_by(id) %>% 
  summarise(abl_mean = mean(pctabl40), 
            abl_sd = sd(pctabl40))
abl_summary[is.na(abl_summary)] <- 0


df_list <- list(tldat, pocket_summary, mobil_summary, abl_summary, diff_summary)

tldat_merged <- df_list %>% reduce(full_join, by='id')
tldat_merged[is.na(tldat_merged)] <- 0

set.seed(123456789) 
tldat_new <- tldat_merged %>% 
  group_by(id) %>% 
  mutate(visit = 1:n(),
         id_new = cur_group_id(), 
         baseage_new = baseage + rnorm(1, 0, sd(tldat_summary$baseage)/2), 
         basebmi_new = basebmi + rnorm(1, 0, sd(tldat_summary$basebmi)/2), 
         pctpocket5mm_new = rnorm(visits, pocket_mean, pocket_sd/2),
         pctmobil05mm_new = rnorm(visits, mobil_mean, mobil_sd/2),
         pctabl40_new = rnorm(visits, abl_mean, abl_sd/2), 
         diff_new = rnorm(visits, diff_mean, diff_sd/2), 
         tstop_new = cumsum(diff_new), 
         tstop_new = case_when(tstop_new < 0 ~ 0, 
                               tstop_new >=0 ~ tstop_new),
         tstart_new = dplyr::lag(tstop_new))  %>% 
  mutate(pctpocket5mm_new = case_when(pctpocket5mm==0 ~ 0, 
                                      pctpocket5mm!=0 ~ pctpocket5mm_new), 
         pctmobil05mm_new = case_when(pctmobil05mm==0 ~ 0, 
                                      pctmobil05mm!=0 ~ pctmobil05mm_new), 
         pctabl40_new = case_when(pctabl40==0 ~ 0, 
                                  pctabl40!=0 ~ pctabl40_new)) 
tldat_new[is.na(tldat_new)] <- 0


tldat_new2 <- tldat_new %>% 
  dplyr::select(id, tstop_new, TL, baseage_new, college, basesmoke, basebmi_new, basenumteeth, 
         pctpocket5mm_new, pctmobil05mm_new, pctabl40_new)

tldat_new0 <- tldat_new2 %>% 
  group_by(id) %>% 
  filter(row_number()==1) %>% 
  mutate(tstop_new=0, TL=0)
  
tldat_new2 <- rbind(tldat_new2, tldat_new0) %>% 
  arrange(id, tstop_new)

colnames(tldat_new2) <- c("id", "year", "TL", "baseage", "college", "basesmoke", "basebmi", "basenumteeth", 
                          "pctpocket5mm", "pctmobil05mm", "pctabl40")

write.csv(tldat_new2, "vadls_jm_person_synthetic.csv", row.names=FALSE)


#### Tooth level data with some noise ####

##### original code from create_toothleveldata.R #####
vadls <- read.csv("vadls_jm.csv", header = TRUE)

vadls_toothloss <- vadls %>% 
  mutate(idtooth = as.numeric(paste0(id, tooth))) %>% 
  group_by(idtooth) %>%
  mutate(tloss = ifelse(toothstat == 0, 1, 0), 
         diff = obstime - dplyr::lag(obstime, default = first(obstime))) %>%
  filter(toothstat == 0)

vadls_toothnotloss <- vadls %>%
  mutate(idtooth = as.numeric(paste0(id, tooth))) %>% 
  group_by(idtooth) %>%
  mutate(tloss = ifelse(toothstat == 0, 1, 0), 
         diff = obstime - dplyr::lag(obstime, default = first(obstime))) %>%
  filter(toothstat == 1)

vadls <- rbind(vadls_toothnotloss, vadls_toothloss) %>% 
  arrange(id, tooth) %>% 
  group_by(idtooth) %>% 
  mutate(visits = n(),
         visit = 1:n()) %>% 
  ungroup() 

vadls_drop <- vadls %>% filter(visits==1) 
length(unique(vadls_drop$id))
length(unique(vadls_drop$idtooth))
drop_idtooth <- unique(vadls_drop$idtooth)

vadls <- vadls %>%
  dplyr::filter(!(idtooth %in% drop_idtooth))
length(unique(vadls$id))
length(unique(vadls$idtooth))


valong <- vadls


vasurv <- vadls %>%
  dplyr::group_by(id, tooth) %>%
  mutate(#basecal = maxcal[1L],
    #baseppd = pocket[1L],
    #baseattl = maxattl[1L],
    basepock = maxpock[1L],
    #baseplaq = anyplaque[1L],
    #basebld = bleed[1L],
    basesmoke = smoking[1L],
    basebmi = bmi[1L],
    baseage = age[1L]) %>%
  arrange(obstime) %>%
  slice(n()) %>%
  ungroup() 

vasurv_NAs <- vasurv[is.na(vasurv$basebmi),]
# ID 1052, 1341 and 2134 do not have basebmi, the survival submodel will remove them. 
# Thus, the jmlong data should also remove them to keep the same IDs between jmlong and jmsurv
vasurv <- vasurv %>% filter(! (id %in% c(1052, 1341, 2134)))
table(vasurv$tloss)

valong <- valong %>% filter(!(id %in% c(1052, 1341, 2134)))


tdsurv <- tmerge(vasurv, vasurv, id = idtooth, endpt = event(time, tloss) )
tdsurv <- tmerge(tdsurv, vadls, id = idtooth,
                 #age = tdc(obstime, age),
                 #maxcal = tdc(obstime, maxcal),
                 #pocket = tdc(obstime, pocket),
                 #maxattl = tdc(obstime, maxattl),
                 maxpock = tdc(obstime, maxpock),
                 #anyplaque = tdc(obstime, anyplaque),
                 #bleed = tdc(obstime, bleed),
                 smoking = tdc(obstime, smoking),
                 bmi = tdc(obstime, bmi))

# vasurv, valong, tdsurv have the same IDs and tooth
length(unique(vasurv$id))
length(unique(valong$id))
length(unique(tdsurv$id))
length(unique(vasurv$idtooth))
length(unique(valong$idtooth))
length(unique(tdsurv$idtooth))


# the maxpock of following idtooth are all NAs! 
valong_onlyNA <- valong %>% 
  group_by(idtooth) %>% 
  mutate(NAcounts = sum(is.na(maxpock))) %>%
  filter(NAcounts == visits)
length(unique(valong_onlyNA$idtooth))

# here remove those idtooth
vasurv <- vasurv %>% 
  filter(!(idtooth %in% valong_onlyNA$idtooth))
valong <- valong %>%
  filter(!(idtooth %in% valong_onlyNA$idtooth))
tdsurv <- tdsurv %>% 
  filter(!(idtooth %in% valong_onlyNA$idtooth))

length(unique(vasurv$id))
length(unique(valong$id))
length(unique(tdsurv$id))
length(unique(vasurv$idtooth))
length(unique(valong$idtooth))
length(unique(tdsurv$idtooth))




# restrict patients with at least 26 teeth at baseline
id26teeth <- vasurv %>%
  filter(basenumteeth > 25) %>%
  dplyr::select(id) %>%
  distinct()
vasurv <- vasurv %>% 
  filter(id %in% id26teeth$id) %>%
  group_by(id) %>%
  mutate(newid = cur_group_id())
valong <- valong %>%
  filter(id %in% id26teeth$id) %>%
  group_by(id) %>%
  mutate(newid = cur_group_id())
tdsurv <- tdsurv %>% 
  filter(id %in% id26teeth$id) %>%
  group_by(id) %>%
  mutate(newid = cur_group_id())


##### add some noise to the datasets #####
tldat_summary <- tldat %>% 
  group_by(id) %>% 
  summarise(baseage = mean(baseage), 
            basebmi = mean(basebmi))

tldat_summary[is.na(tldat_summary)] <- 0


diff_summary <- valong %>% 
  group_by(id) %>% 
  filter(diff!=0) %>%
  summarise(diff_mean = mean(diff), 
            diff_sd = sd(diff))
diff_summary[is.na(diff_summary)] <- 0 

diff_new_dat <- tldat_new %>%
  dplyr::select(id, visit, diff_new, baseage_new, basebmi_new)

valong_new <- valong %>% 
  left_join(diff_new_dat, by=c("id", "visit")) %>%
  group_by(idtooth) %>%
  mutate(diff_new = coalesce(diff_new, diff), 
         diff_new = lag(diff_new), 
         diff_new = coalesce(diff_new, 0),
         baseage_new = coalesce(baseage_new, baseage_new[1]),
         basebmi_new = coalesce(basebmi_new, basebmi_new[1]),
         obstime_new = cumsum(diff_new),
         time_new = tail(obstime_new,1))

vasurv_new <- valong_new %>%
  dplyr::group_by(id, tooth) %>%
  mutate(basepock = maxpock[1L],
    basesmoke = smoking[1L],
    basebmi = basebmi_new[1L],
    baseage = baseage_new[1L]) %>%
  arrange(obstime_new) %>%
  slice(n()) %>%
  ungroup() 

tdsurv_new <- tmerge(vasurv_new, vasurv_new, id = idtooth, endpt = event(time_new, tloss) )
tdsurv_new <- tmerge(tdsurv_new, valong_new, id = idtooth,
                 maxpock = tdc(obstime_new, maxpock),
                 smoking = tdc(obstime_new, smoking),
                 bmi = tdc(obstime_new, basebmi_new))

# vasurv, valong, tdsurv have the same IDs and tooth
length(unique(vasurv_new$id))
length(unique(valong_new$id))
length(unique(tdsurv_new$id))
length(unique(vasurv_new$idtooth))
length(unique(valong_new$idtooth))
length(unique(tdsurv_new$idtooth))



valong_new2 <- valong_new %>% 
  dplyr::select(id, tooth, idtooth, time_new, obstime_new, toothstat, tloss, maxpock, baseage_new,
                smoking, college, basebmi_new, basenumteeth)
colnames(valong_new2) <- c("id", "tooth", "idtooth", "time", "obstime", "toothstat", "tloss", "maxpock", "baseage",
                           "smoking", "college", "basebmi", "basenumteeth")

write.csv(valong_new2, "vadls_jm_tooth_synthetic.csv", row.names=FALSE)





