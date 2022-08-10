library(dplyr)

#Read Data
TP <- read.csv("SouthernLM_TP.csv")

mean_TP = TP %>%
  group_by(YEAR) %>%
  dplyr::summarize(Mean = mean(TP_ug_L, na.rm=TRUE), n=n())

