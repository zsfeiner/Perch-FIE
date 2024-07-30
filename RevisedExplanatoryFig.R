#TP and temp Summary
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyverse)

###Abiotic variables
#Read Data
TP <- read.csv("SouthernLM_TP.csv")

TP_Year = TP %>%
  rename(Year=YEAR) %>%
  group_by(Year) %>%
  dplyr::summarize(MeanTP = mean(TP_ug_L, na.rm=TRUE), sdTP = sd(TP_ug_L, na.rm=TRUE)) %>%
  filter(Year < 2017)



temps <- read_csv("StJoeTemps_1960_2018.csv")

##Summarize a few different ways on annual basis

#Create winter year variable to subtract 1 from Jan and Feb years
##align to first winter experience
temps$WinterYear <- ifelse(temps$Month %in% c(1:2), temps$Year - 1, temps$Year)

#Create GDD5 tracker
temps$GDD5 <- ifelse(temps$TempC > 5, temps$TempC-5, 0)

GDD <- temps %>%
  group_by(Year) %>%
  summarize(GDD5Annual = sum(GDD5)) %>%
  filter(Year > 1982 & Year < 2017)

abiotic_graph <- left_join(GDD, TP_Year)

abiotic_graph

#Biotic variables
female_yep <- read.csv("female_yep.csv", 
                       colClasses=c("integer","character","integer","integer","character","character","integer","integer","integer","integer"))

#assign cohort year
female_yep$CohortYear=female_yep$YEAR-female_yep$AGE

#assign cohort number
female_yep$Cohort = female_yep$CohortYear - min(female_yep$CohortYear) + 1

#Add relative weight column
female_yep$Ws <- 10^(-5.386 + 3.230 * log10(female_yep$LENGTH))
female_yep$Wr <- female_yep$WEIGHT / female_yep$Ws * 100

#Combine explanatory variables for figure
fig_graph_Data <- data.frame(TL = female_yep$LENGTH,
                             WR = female_yep$Wr,
                             AGE = female_yep$AGE,
                             YEAR = female_yep$YEAR,
                             MAT = female_yep$MAT) %>%
  dplyr::group_by(YEAR) %>%
  dplyr::summarise("Median Age" = median(AGE,na.rm=TRUE),
                   SDAGE = sd(AGE,na.rm=TRUE),
                   "Median Wr" = median(WR,na.rm=TRUE),
                   SDWR = sd(WR,na.rm=TRUE),
                   "Median TL" = median(TL,na.rm=TRUE),
                   SDTL = sd(TL,na.rm=TRUE),
                   "Sample size" = dplyr::n())  %>%
  full_join(abiotic_graph, by=c("YEAR"="Year")) %>%
  pivot_longer(cols=c("Median Age","Median Wr","Median TL","Sample size", "GDD5Annual", "MeanTP"))
fig_graph_Data

fig_graph_Data <- fig_graph_Data %>%
  mutate(name2 = factor(name, levels=c("Sample size","Median Age","Median TL","Median Wr", "GDD5Annual", "MeanTP"),
                        #labels=c("Sample size (N)","Median age (yr)","Median TL (mm)", "Median Wr (%)", expression(Mean*GDD[5]), "Mean TP (ug/L)"))
                        labels=c(expression("`Sample size (N)`"),expression("`Median age (yr)`"),expression("`Median TL (mm)`"), 
                                 expression("`Median Wr (%)`"), expression(Mean~GDD[5]~"("*degree*"C)"), expression(Mean~TP~"("*mu*"g/L)"))))


Predictors_fig <- ggplot(fig_graph_Data, aes(x=YEAR, y=value)) + 
  theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  #geom_line(aes(x =  GLSM_Internal_corr$Date_2, y = Pred_obs$Q5),linetype="dashed") +
  #geom_line(aes(x =  GLSM_Internal_corr$Date_2, y = Pred_obs$Q95),linetype="dashed") + 
  #geom_point(aes(x = Pred_obs$Date, y = Pred_obs$Q50), size=5) +
  geom_line(size=1) +
  #geom_errorbar(mapping=aes(x = COHORT,ymin=value + SD,ymax=Q95))+
  facet_wrap(~name2, scales = "free_y", ncol=2, 
             strip.position="left", labeller=label_parsed) +
  #geom_ribbon(aes(x = GLSM_Internal_corr$Date_2,ymin=Pred_obs$Q10,ymax=Pred_obs$Q10),fill="gray")+
  ylab("")+
  xlab("Year")+
  #ylim(0,6.75) +
  theme(axis.text.x=element_text(size=18,angle=90,hjust=1, vjust=1),
        axis.text.y=element_text(size=18),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        panel.border = element_rect(color="black"),
        axis.line = element_line(colour = "black"),
        strip.text = element_text(size=18),
        # axis.ticks.x = element_blank(),
        # axis.text.x = element_blank(),
        legend.position="none",
        strip.background=element_blank(),
        strip.placement="outside") + 
  scale_x_continuous(breaks=seq(1980,2020,5)) + 
  geom_vline(xintercept=1997, lty=2)

Predictors_fig

ggsave(file="./Figures/RevisedPredictors_Fig.svg", plot=Predictors_fig, width=12, height=10)
ggsave(file="./Figures/RevisedPredictors_Fig.png", plot=Predictors_fig, width=12, height=10)


#Plot GDD and TP history for a cohort
cohortsummary <- female_yep %>%
  group_by(CohortYear, YEAR, AGE) %>%
  summarize(N=n()) %>%
  left_join(abiotic_graph, by=c("YEAR"="Year")) %>%
  mutate(Fishing = ifelse(YEAR <=1996, 1, 0))
  
cohortsummary  

GDDcohortsplot <- ggplot(cohortsummary, aes(x=YEAR, y=GDD5Annual, color=CohortYear)) + 
  geom_line() + geom_point() + theme_bw() +
  facet_wrap(~CohortYear) + xlab("Year") + ylab(expression(Mean~GDD[5]~"("*degree*"C)")) + 
  scale_color_continuous(name="Cohort")
GDDcohortsplot

ggsave(file="./Figures/Revised_Supplement_GDDbyCohort_Fig.png", plot=GDDcohortsplot, width=12, height=10)

TPCohortsplot <- ggplot(cohortsummary, aes(x=YEAR, y=MeanTP, color=CohortYear)) + 
  geom_line() + geom_point() + 
  facet_wrap(~CohortYear) + theme_bw() + xlab("Year") + ylab(expression(Mean~TP~"("*mu*"g/L)")) +
  scale_color_continuous(name="Cohort")
TPCohortsplot
ggsave(file="./Figures/Revised_Supplement_TPbyCohort_Fig.png", plot=TPCohortsplot, width=12, height=10)


