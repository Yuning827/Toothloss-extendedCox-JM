## descriptive analysis -- person level
library(ggplot2)
library(dplyr)
library(tidyverse)
library(lattice)

tldat <- read.csv("vadls_jm_personlevel.csv", header = TRUE)


va.samp <- subset(tldat, id %in% c(144,286,364,369)) %>%
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
