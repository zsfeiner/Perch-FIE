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
                             COHORT = female_yep$CohortYear,
                             MAT = female_yep$MAT) %>%
                  dplyr::group_by(COHORT) %>%
                  dplyr::summarise("Median Age" = median(AGE,na.rm=TRUE),
                                   SDAGE = sd(AGE,na.rm=TRUE),
                                   "Median Wr" = median(WR,na.rm=TRUE),
                                   SDWR = sd(WR,na.rm=TRUE),
                                   "Median TL" = median(TL,na.rm=TRUE),
                                   SDTL = sd(TL,na.rm=TRUE),
                                   "Sample size" = dplyr::n()) %>%
                pivot_longer(cols=c("Median Age","Median Wr","Median TL","Sample size"))


Predictors_fig <- ggplot(fig_graph_Data) + 
                          theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
                          #geom_line(aes(x =  GLSM_Internal_corr$Date_2, y = Pred_obs$Q5),linetype="dashed") +
                          #geom_line(aes(x =  GLSM_Internal_corr$Date_2, y = Pred_obs$Q95),linetype="dashed") + 
                          #geom_point(aes(x = Pred_obs$Date, y = Pred_obs$Q50), size=5) +
                          geom_line(aes(x = COHORT, y = value),size=2) +
                          #geom_errorbar(mapping=aes(x = COHORT,ymin=value + SD,ymax=Q95))+
                          facet_wrap(~name, scales = "free_y", ncol=2) +
                          #geom_ribbon(aes(x = GLSM_Internal_corr$Date_2,ymin=Pred_obs$Q10,ymax=Pred_obs$Q10),fill="gray")+
                          ylab("Response")+
                          xlab("Cohort")+
                          #ylim(0,6.75) +
                          theme(axis.text.x=element_text(size=18,angle=45,hjust=1, vjust=1),
                                axis.text.y=element_text(size=18),
                                axis.title.x=element_text(size=20),
                                axis.title.y=element_text(size=20),
                                panel.border = element_blank(),
                                axis.line = element_line(colour = "black"),
                                strip.text = element_text(size=18),
                                # axis.ticks.x = element_blank(),
                                # axis.text.x = element_blank(),
                                legend.position="none")

Predictors_fig

ggsave(file="./Figures/Predictors_Fig.svg", plot=Predictors_fig, width=16, height=10)
ggsave(file="./Figures/Predictors_Fig.png", plot=Predictors_fig, width=16, height=10)

