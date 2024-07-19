library(tidyr)
library(ggplot2)
library(matrixStats)

#Coefficient plot

beta <- data.frame(rstan::extract(fit_full_fishing,pars='beta'))


names(beta) <- c("Intercept","Age","Total Length","Relative Weight","Age * Total Length","Age * Relative Weight", "Total Phosphorus","Growing Degree Days","Commercial Fishing")

beta_long <- beta %>%
              pivot_longer(cols=c("Intercept","Age","Total Length","Relative Weight","Age * Total Length","Age * Relative Weight","Total Phosphorus","Growing Degree Days","Commercial Fishing"),
                           names_to = "Coefficient",
                           values_to = "Estimate")

cbPalette=c("#99CC00","#CCFF33","#333300","#999933","#FFFF33","#FFFF99","#DDFFCC","#66CC99","#999933")

level_order <- c("Intercept","Age","Total Length","Relative Weight","Commercial Fishing","Growing Degree Days", "Total Phosphorus","Age * Total Length","Age * Relative Weight")

MargCoefPlot<-(ggplot(data=beta_long,aes(x=factor(Coefficient,level=level_order),y=Estimate,fill=factor(Coefficient,level=level_order)))+
                 theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
                 geom_violin(trim=T,scale="width")+
                 #scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9","#999999","#E69F00"))+
                 #scale_fill_brewer(palette="Dark2")+
                 scale_fill_manual(values=cbPalette)+
                 stat_summary(fun = median,
                              fun.min = function(x) quantile(x, probs = 0.025),
                              fun.max = function(x) quantile(x, probs = 0.975), size = 0.75)+
                 #facet_wrap(~ Param,ncol=3)+
                 ylab("Standardized coefficient")+
                 xlab("Predictor")+
                 geom_hline(yintercept=0)+
                 #ylim(0.20,1.00)+
                 #scale_x_continuous(breaks=c(2014, 2015, 2016, 2017, 2018))+
                 theme(axis.text.x=element_text(size=15,angle=90),
                       axis.text.y=element_text(size=20),
                       axis.title.x=element_text(size=20),
                       axis.title.y=element_text(size=20,angle=90),
                       panel.border = element_blank(),
                       axis.line = element_line(colour = "black"),
                       #axis.ticks.x = element_blank(),
                       #axis.text.x = element_blank()
                       legend.position="none",
                       strip.text = element_text(size=20)
                 )
)
MargCoefPlot

ggsave(file="./Figures/Marg_coef.svg", plot=MargCoefPlot, width=15, height=10)
ggsave(file="./Figures/Marg_coef.png", plot=MargCoefPlot, width=15, height=10)


##Cohort effects
beta_u <- data.frame(rstan::extract(fit_full_fishing,pars='u'))
beta_u_1 <- list()

#extract cohorts for each beta
beta_u_1[[1]] <- beta_u[,seq(from = 1, to = 289, by = 8)]
beta_u_1[[2]] <- beta_u[,seq(from = 2, to = 290, by = 8)]
beta_u_1[[3]] <- beta_u[,seq(from = 3, to = 291, by = 8)]
beta_u_1[[4]] <- beta_u[,seq(from = 4, to = 292, by = 8)]
beta_u_1[[5]] <- beta_u[,seq(from = 5, to = 293, by = 8)]
beta_u_1[[6]] <- beta_u[,seq(from = 6, to = 294, by = 8)]
beta_u_1[[7]] <- beta_u[,seq(from = 7, to = 295, by = 8)]
beta_u_1[[8]] <- beta_u[,seq(from = 8, to = 296, by = 8)]
#beta_u_1[[9]] <- beta_u[,seq(from = 9, to = 333, by = 9)]

#Add mean beta to each cohort effect
cohort_beta_1 <- list()
cohort_beta_1[[1]] <- beta_u_1[[1]] + beta[,1]
cohort_beta_1[[2]] <- beta_u_1[[2]] + beta[,2]
cohort_beta_1[[3]] <- beta_u_1[[3]] + beta[,3]
cohort_beta_1[[4]] <- beta_u_1[[4]] + beta[,4]
cohort_beta_1[[5]] <- beta_u_1[[5]] + beta[,5]
cohort_beta_1[[6]] <- beta_u_1[[6]] + beta[,6]
cohort_beta_1[[7]] <- beta_u_1[[7]] + beta[,7]
cohort_beta_1[[8]] <- beta_u_1[[8]] + beta[,8]
#cohort_beta_1[[9]] <- beta_u_1[[9]] + beta[,9]

#Combine by column
cohort_beta_all <- cbind(cohort_beta_1[[1]],cohort_beta_1[[2]],cohort_beta_1[[3]],cohort_beta_1[[4]],cohort_beta_1[[5]],
                        cohort_beta_1[[6]],cohort_beta_1[[7]],cohort_beta_1[[8]])

#Calcuate quantiles and add chort number
annual_effects <- data.frame((colQuantiles(as.matrix(cohort_beta_all),probs=c(0.025,0.50,0.975))))
annual_effects$cohort <- rep(seq(1979,2015,1),8)

#Add column of beta names
covar_names <- c("Intercept","Age","Total Length","Relative Weight","Age * Total Length","Age * Relative Weight", "Total Phosphorus","Growing Degree Days")
annual_effects$parm <- c(rep(covar_names[1],37),rep(covar_names[2],37),rep(covar_names[3],37),rep(covar_names[4],37),
                        rep(covar_names[5],37),rep(covar_names[6],37),rep(covar_names[7],37),rep(covar_names[8],37))#,
                        #rep(covar_names[9],37))

#label quantiles and specify the order for facets
names(annual_effects) =c("LCL","MED","UCL","Cohort","Param")
level_order <- c("Intercept","Age","Total Length","Relative Weight","Growing Degree Days", "Total Phosphorus","Age * Total Length","Age * Relative Weight")
annual_effects <- annual_effects %>%
                  mutate(across(Param, factor, levels=level_order))

#Facet plot
CoefPlot_subset<-(ggplot(data=annual_effects,aes(x=Cohort,y=MED))+
                    theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
                    geom_point(aes(x=Cohort, y = MED)) +
                    geom_errorbar(aes(x=Cohort, ymin=LCL, ymax=UCL), width=0.05) + 
                    facet_wrap(~ Param,ncol=3)+
                    ylab("Standardized coefficient")+
                    xlab("Cohort")+
                    geom_hline(yintercept=0)+
                    #ylim(0.20,1.00)+
                    #scale_x_continuous(breaks=c(2,4,6,8,10,12))+
                    #scale_x_continuous(breaks=c(2014, 2015, 2016, 2017, 2018))+
                    theme(axis.text.x=element_text(size=20),
                          axis.text.y=element_text(size=20),
                          axis.title.x=element_text(size=20),
                          axis.title.y=element_text(size=20,angle=90),
                          panel.border = element_blank(),
                          axis.line = element_line(colour = "black"),
                          #axis.ticks.x = element_blank(),
                          #axis.text.x = element_blank()
                          legend.position="none",
                          strip.text = element_text(size=20)
                    )
)
CoefPlot_subset

ggsave(file="./Figures/Coef_By_Cohort.svg", plot=CoefPlot_subset, width=14, height=12)
ggsave(file="./Figures/Coef_By_Cohort.png", plot=CoefPlot_subset, width=14, height=12)
