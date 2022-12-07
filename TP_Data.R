library(dplyr)

#Read Data
TP <- read.csv("SouthernLM_TP.csv")

mean_TP_byYear = TP %>%
  rename(Year=YEAR) %>%
  group_by(Year) %>%
  dplyr::summarize(Mean = mean(TP_ug_L, na.rm=TRUE), n=n())

