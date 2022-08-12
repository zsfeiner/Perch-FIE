library(tidyverse)

temps <- read_csv("StJoeTemps_1960_2018.csv")
temps


##Summarize a few different ways on annual basis
#Could also do lagging (3 year?) conditions to reflect pre-maturation early life experience?
##Days 0-4C
##Mean winter (Dec-Feb) temps
##Mean summer (Jun-Aug) temps
##GDD5 winter
##GDD5 summer
##Mean annual
##GDD5 annual

#Create winter year variable to subtract 1 from Jan and Feb years
##align to first winter experience
temps$WinterYear <- ifelse(temps$Month %in% c(1:2), temps$Year - 1, temps$Year)

#Create GDD5 tracker
temps$GDD5 <- ifelse(temps$TempC > 5, temps$TempC-5, 0)

winter.index <- temps %>%
  group_by(WinterYear) %>%
  summarize(Days0_4C = sum(TempC < 4 & TempC > 0),
            MeanWint = mean(TempC[which(Month %in% c(12,1,2))]),
            WinterGDD5 = sum(GDD5[which(Month %in% c(12,1,2))]))

other.index <- temps %>%
  group_by(Year) %>%
  summarize(MeanSumm = mean(TempC[which(Month %in% c(6,7,8))]),
            GDD5Summ = sum(TempC[which(Month %in% c(6,7,8))]),
            MeanAnnual = mean(TempC),
            GDD5Annual = sum(GDD5))



#Combine
temps.summary <- left_join(other.index, winter.index, by=c("Year"="WinterYear"))
temps.summary

#write.csv(temps.summary, "StJoeTempSummary.csv", row.names=F)
