vadls <- read.csv("vadls_jm_tooth_synthetic.csv", header = TRUE)


valong <- vadls


vasurv <- vadls %>%
  dplyr::group_by(id, tooth) %>%
  mutate(basepock = maxpock[1L],
    basesmoke = smoking[1L],
    basebmi = basebmi[1L],
    baseage = baseage[1L]) %>%
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
                 maxpock = tdc(obstime, maxpock),
                 smoking = tdc(obstime, smoking),
                 bmi = tdc(obstime, basebmi))

# vasurv, valong, tdsurv have the same IDs and tooth
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



