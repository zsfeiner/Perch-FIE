#TP and temp Summary
library(dplyr)

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

abiotic_graph <- cbind(GDD, TP_Year)
abiotic_graph <- abiotic_graph[,-1]

abiotic_plot<-(ggplot(abiotic_graph, aes(x=Year))+
                 theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
                 geom_line(aes(y=GDD5Annual),size=1)+
                 geom_line(aes(y=MeanTP*1000),size=1, linetype="dashed")+
                 # ylab("Growing degree days")+
                 # xlab("Year")+
                 #ylim(0.20,1.00)+
                 scale_y_continuous(name = "Growing degree dadys", 
                                    sec.axis = sec_axis(~./1000, name = "TP \u00B5g/l")) +
                 scale_x_continuous(breaks=c(seq(from = 1984, to = 2016, by = 4)))+
                 theme(axis.text.x=element_text(size=36),
                       axis.text.y=element_text(size=36),
                       axis.title.x=element_text(size=36),
                       axis.title.y=element_text(size=36,angle=90),
                       panel.border = element_blank(),
                       axis.line = element_line(colour = "black"),
                       #axis.ticks.x = element_blank(),
                       #axis.text.x = element_blank()
                       legend.position="none"
                 )
)
abiotic_plot

ggsave(file="./Figures/abiotic.svg", plot=abiotic_plot, width=16, height=10)
ggsave(file="./Figures/abiotic.png", plot=abiotic_plot, width=16, height=10)
