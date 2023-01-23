## descriptive analysis -- person level
library(lattice)

tldat <- read.csv("vadls_jm_person_synthetic.csv", header = TRUE)

# number of patients
length(unique(tldat$id))

# number of patients who lost at least one tooth
length(unique(tldat[which(tldat$TL == 1),]$id))

# median length to first tooth loss
median(tldat[which(tldat$TL == 1),]$year)

# maximum number of years of follow-up
max(tldat$year)


va.samp <- subset(tldat, id %in% c(1510,1956,2199,2212)) %>%
  pivot_longer(cols = c(pctpocket5mm, pctabl40, pctmobil05mm),
               names_to = "biomarkers",
               values_to = "percentage") %>%
  mutate(biomarkers=as.factor(biomarkers),
         id=as.factor(id))

levels(va.samp$biomarkers) <- c("Percent ABL >= 40%", "Percent MOB >= 0.5mm", "Percent PPD >= 5mm")
levels(va.samp$id) <- c("Patient 1", "Patient 2", "Patient 3", "Patient 4")


ggplot(va.samp, aes(y=percentage, x=year, color=biomarkers)) + 
  geom_point(size=2) + 
  geom_line() + 
  labs(y = "Percent", x = "Years since baseline") +
  facet_grid(biomarkers ~ id) + 
  scale_color_brewer(palette = "Dark2") + 
  theme_bw(base_size = 12) +
  theme(legend.position = "none")
