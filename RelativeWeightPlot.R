#Relative Weight by year
RW_Cohort <- Cohort+1978
rw_Year = data.frame(cbind(RW,RW_Cohort)) %>%
  dplyr::group_by(RW_Cohort) %>%
  dplyr::summarize(MeanRW = mean(RW, na.rm=TRUE), sdRW = sd(RW, na.rm=TRUE), n = n(), se = sd(RW, na.rm=TRUE)/sqrt(n()))

rw_plot<-(ggplot(rw_Year, aes(x=RW_Cohort,y = MeanRW  ))+
                 theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
                 geom_line(linewidth=1)+
                 # geom_errorbar(aes(ymin=MeanRW-(2*se), ymax=MeanRW+(2*se)), width=1,linewidth=1,
                 #          position=position_dodge(.9)) +
          # geom_line(aes(y=GDD5Annual),size=1)+
                 # geom_line(aes(y=MeanTP*1000),size=1, linetype="dashed")+
                 ylab("Mean relative weight")+
                 xlab("Cohort year")+
                 #ylim(0.20,1.00)+
                 # scale_y_continuous(name = "Growing degree dadys", 
                 #                    sec.axis = sec_axis(~./1000, name = "TP \u00B5g/l")) +
                 scale_x_continuous(breaks=c(seq(from = 1979, to = 2016, by = 4)))+
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
rw_plot

ggsave(file="./Figures/rw_plot.svg", plot=rw_plot, width=16, height=10)
ggsave(file="./Figures/rw_plot.png", plot=rw_plot, width=16, height=10)
