library(rstan)
library(coda)
library(tidyverse)
library(ggplot2)

#Run TP_Data.R to load TP data - will be ignored in model
source("./Scripts/TP_Data.R")
#Run SummarizeWaterTemp.R to load water temp data
source("./Scripts/SummarizeWaterTemps.R")

female_yep <- read.csv("./Data/female_yep.csv", 
                       colClasses=c("integer","character","integer","integer","character","character","integer","integer","integer","integer"))

#assign cohort year
female_yep$CohortYear=female_yep$YEAR-female_yep$AGE


#assign cohort number
female_yep$Cohort = female_yep$CohortYear - min(female_yep$CohortYear) + 1

#Calculate mean for standardizing
female_meanTL = mean(female_yep$LENGTH)

#Calculate SD for standardizing
female_SDTL = sd(female_yep$LENGTH)

#####Do some comparisons of linear vs log-linear model fit#####
plot(female_yep$LENGTH~female_yep$AGE)
plot(LENGTH ~ AGE, data=filter(female_yep, AGE<=6))

mod1 <- lm(LENGTH ~ AGE, data=filter(female_yep, AGE<=6))
mod2 <- lm(LENGTH ~ log(AGE), data=filter(female_yep, AGE<=6))
mod4 <- nls(LENGTH ~ Linf * (1-exp(-k*(AGE-t0))), data=filter(female_yep, AGE <= 6),
            start=list(Linf=300, k=0.2, t0=0))
AIC(mod1, mod2, mod4)
summary(mod2)
summary(mod4)

sqrt(mean((female_yep$LENGTH[female_yep$AGE<=6] - predict(mod1))^2))
sqrt(mean((female_yep$LENGTH[female_yep$AGE<=6] - predict(mod2))^2))
sqrt(mean((female_yep$LENGTH[female_yep$AGE<=6] - predict(mod4))^2))

linecolors <- c("Linear: AIC = 46874, RMSE = 39.51"="blue",
                "Linear-log: AIC = 46541, RMSE = 38.10" = "red",
                "von Bertalanffy: AIC = 46512, RMSE = 37.97" = "green")

growth_model_plot <- ggplot(filter(female_yep, AGE<=6), aes(x=AGE, y=LENGTH)) +
  geom_point() + 
  geom_smooth(method="lm", se=F, lwd=2, aes(color="Linear: AIC = 46874, RMSE = 39.51")) + 
  geom_smooth(method="lm", formula="y~log(x)", se=F, lwd=2, aes(color="Linear-log: AIC = 46541, RMSE = 38.10")) + 
  geom_smooth(method="nls", 
              formula="y~Linf * (1-exp(-k*(x-t0)))", 
              method.args=list(start=c(Linf=300, k=0.2, t0=0)), se=F, lwd=2,
              aes(color="von Bertalanffy: AIC = 46512, RMSE = 37.97")) + 
  theme_classic() + labs(y="Total length (mm)", x="Age", color="Legend") + scale_color_manual(name="Model",values=linecolors) + 
  theme(legend.position=c(0.75,0.1), legend.text=element_text(size=14), legend.title=element_text(size=14),
        axis.title=element_text(size=14), axis.text=element_text(size=12))

#Create supplemental figure for growth model
ggsave(file="./RevisedFigures/GrowthModelPlot.png", plot=growth_model_plot, width=12, height=10)

####Use log(age) model - similar fits, simpler implementation###

###Length/Weight
#Slow way of getting missing weight, only 6 NA's so it doesn't take too long
for (x in 1:nrow(female_yep)){
  if (is.na(female_yep$WEIGHT[x])==TRUE){
    subsetdata=subset(female_yep,female_yep$YEAR==female_yep$YEAR[x])
    lm1=lm(log10(subsetdata$WEIGHT)~log10(subsetdata$LENGTH))
    female_yep$WEIGHT[x] = 10^(coef(lm1)[1] + coef(lm1)[2]* log10(female_yep$LENGTH[x]) )
  }
}
plot(female_yep$WEIGHT~female_yep$LENGTH)

#Confirm no NA's
table(is.na(female_yep$WEIGHT))

#Add relative weight column using Willis et al standard weight
female_yep$Ws <- 10^(-5.386 + 3.230 * log10(female_yep$LENGTH))
female_yep$Wr <- female_yep$WEIGHT / female_yep$Ws * 100


#####Add environmental variables - use mean TP and annual GDD5#####
mean_TP_byYear
temps.summary

#Combine TP, temp
GDD_byYear <- temps.summary %>%
  select("GDD5Annual","Year")

#Create abiotics data frame with TP, Fishing, and GDD by year and cohort year
abiotics <- GDD_byYear %>%
  inner_join(mean_TP_byYear) %>%
  mutate(CohortYear = Year - min(female_yep$CohortYear) + 1) %>%
  filter(Year >= min(female_yep$YEAR) - 1) %>%
  rename(TP=Mean) %>%
  mutate(Fishing = ifelse(Year <= 1996, 1, 0)) %>%
  select(Year, CohortYear, GDD5Annual, TP, Fishing)

#Join perch data up to temps, TP, and binary fishing
female_yep <- left_join(female_yep, select(temps.summary, Year, GDD5Annual), by=c("YEAR"="Year")) %>%
  rename("GDD5"="GDD5Annual") %>%
  left_join(select(mean_TP_byYear, Year, Mean), by=c("YEAR"="Year")) %>%
  rename("TP"="Mean") %>%
  mutate(Fishing = ifelse(YEAR <= 1996, 1, 0))


####PMRN model set up######
##Create dataset
#Female YEP
N = nrow(female_yep)

#Explanatory variables
Age <- (female_yep$AGE)
TL.mm <- female_yep$LENGTH
mean_TL <- mean(female_yep$LENGTH) #mean for scaling
sd_TL <- sd(female_yep$LENGTH) #sd for scaling
RW <- female_yep$Wr #body condition, relative weight
mean_RW <- mean(female_yep$Wr) #mean for scaling
sd_RW <- sd(female_yep$Wr) #sd for scaling
TP <- female_yep$TP #total phosphorus
mean_TP <- mean(abiotics$TP) #mean for scaling
sd_TP <- sd(abiotics$TP) #sd for scaling
Fished <- female_yep$Fishing #Binary fishing
GDD5 <- female_yep$GDD5 #GDD base 5C
mean_GDD5 <- mean(abiotics$GDD5Annual) #mean for scaling
sd_GDD5 <- sd(abiotics$GDD5Annual) #sd for scaling
Cohort <- female_yep$Cohort
Year <- female_yep$YEAR - min(female_yep$YEAR) + 1  #Year marker starting from 1 to get prev year's GDD from abiotics
nCohorts <- length(unique(female_yep$Cohort))
X <- cbind(Age,TL.mm,RW,GDD5,TP,Fished)
K <- ncol(X)
Mat <- ifelse(female_yep$MAT=="I",0,1)
nAbiotics <- nrow(abiotics)

#Create datalist
dat <- list('N'=nrow(female_yep), 'K'=K, 'Mat'=Mat, 'Age'=Age, 'TL'=TL.mm, 
            'Cohort'=Cohort, 'nCohorts'=nCohorts, 'mean_TL'=mean_TL, 'sd_TL'=sd_TL,
            'RW'=RW, 'mean_RW'=mean_RW, 'sd_RW'=sd_RW,
            'GDD5'=GDD5, 'mean_GDD5'=mean_GDD5, 'sd_GDD5'=sd_GDD5,
            'TP'=TP, 'mean_TP'=mean_TP, 'sd_TP'=sd_TP,
            'abiotics'=abiotics, 'Year'=Year, 'nAbiotics'=nAbiotics,
            'Fished'=Fished)

#Function to initialize Stan run
inits_full_fishing <- function() {
  beta <- c(rnorm(1,-3,0.05),rnorm(1,1,0.05),rnorm(1,2,0.05),rnorm(1,0,0.05),rnorm(1,0,0.05),rnorm(1,0,0.05),rnorm(1,0,0.05),rnorm(1,0,0.05),rnorm(1,0,0.05))
  sigma_u <- c(runif(1,1,5),runif(1,0.1,1),runif(1,0.5,2),runif(1,0.1,1),runif(1,0.1,1),runif(1,0.1,1),runif(1,0.1,1),runif(1,0.1,1))
  phi_mu <- rnorm(1,100,0.05)
  phi_sigma <- runif(1,0.5,1)
  gamma_mu <- rnorm(1,100,0.05)
  gamma_sigma <- runif(1,0.5,2)
  omega1 <- rnorm(nCohorts,0,1)
  omega2 <- rnorm(nCohorts,0,1)
  omega3 <- rnorm(nCohorts,0,1)
  omega4 <- rnorm(nCohorts,0,1)
  zeta <- runif(1,0.5,2)
  tau <- runif(1,1,10)
  z_u <- matrix(rnorm(8*nCohorts,0,0.5),nrow=8,ncol=nCohorts)
  init <- list("beta"=beta, "sigma_u"=sigma_u, 
               "phi_mu"=phi_mu, "phi_sigma"=phi_sigma, "gamma_mu"=gamma_mu, "gamma_sigma"=gamma_sigma,
               "omega1"=omega1, "omega2"=omega2,"omega3"=omega3,"omega4"=omega4,"zeta"=zeta,
               "tau"=tau, "z_u"=z_u)
  return(init)
}

#Test
#inits_full_fishing()

stanmatcode_full_fishing = stan_model(file = './Scripts/yep_fie_pmrn_revised.stan')
fit_full_fishing = sampling(stanmatcode_full_fishing, data=dat, init=inits_full_fishing, 
                    iter=4000, warmup=2000, thin=1, chains=3, cores=3, #was 4000 and 2000
                    control=list(adapt_delta=0.90,max_treedepth=10) )

#saveRDS(fit_full_fishing,"./Data/YEPFIE_covar_enviro_revised_12.2.2024.RDS")

#Examine convergence and results
print(fit_full_fishing, pars=c('beta','sigma_u','phi_mu','gamma_mu','tau'), digits=3, prob=c(0.025,0.5,0.975))
print(fit_full_fishing, pars=c('beta'), digits=3, prob=c(0.025,0.5,0.975))

print(fit_full_fishing, pars=c('m'), digits=3, prob=c(0.025,0.5,0.975))
print(fit_full_fishing, pars=c('s'), digits=3, prob=c(0.025,0.5,0.975))

stan_trace(fit_full_fishing,pars=c('omega1'))
stan_trace(fit_full_fishing,pars=c('omega2'))
stan_trace(fit_full_fishing,pars=c('omega3'))
stan_trace(fit_full_fishing,pars=c('omega4'))
stan_trace(fit_full_fishing,pars=c('zeta'))
stan_trace(fit_full_fishing,pars=c('p[1]','p[2]','p[3]','p[4]','p[5]','p[6]'))
stan_trace(fit_full_fishing,pars=c('prev_p[1]','prev_p[2]','prev_p[3]','prev_p[4]','prev_p[5]','prev_p[6]'))
stan_trace(fit_full_fishing,pars=c('s[1]','s[2]','s[3]','s[4]','s[5]','s[6]'))
stan_trace(fit_full_fishing,pars=c('prev_sc_wr[1]','prev_sc_wr[2]','prev_sc_wr[3]','prev_sc_wr[4]','prev_sc_wr[5]','prev_sc_wr[6]'))
stan_trace(fit_full_fishing,pars=c('wr_inc[1]','wr_inc[2]','wr_inc[3]','wr_inc[4]','wr_inc[5]','wr_inc[6]'))
stan_trace(fit_full_fishing,pars=c('m[1]','m[2]','m[3]','m[4]','m[5]','m[6]'))
stan_trace(fit_full_fishing,pars=c('beta','gamma_mu','phi_mu'))
stan_trace(fit_full_fishing,pars=c('L_u'))
pairs(fit_full_fishing,pars=c('beta','gamma_mu','phi_mu'))

#Create supplemental trace plots of parameters
#betas and sigmas
mat_trace<-stan_trace(fit_full_fishing,pars=c('beta','sigma_u'))
mat_trace
ggsave(file="./RevisedFigures/MatModelTraceplot.png", plot=mat_trace, width=16, height=10)

#phi growth 
trace_growth_phi <- stan_trace(fit_full_fishing, pars=c('phi'))
trace_growth_phi
ggsave(file="./RevisedFigures/Growth_Phi_Traceplot.png", plot=trace_growth_phi, width=12, height=10)

#gamma growth
trace_growth_gamma <- stan_trace(fit_full_fishing, pars=c('gamma'))
trace_growth_gamma
ggsave(file="./RevisedFigures/Growth_Gamma_Traceplot.png", plot=trace_growth_gamma, width=12, height=10)

#Growth hyperparameters
trace_growth_hyper <- stan_trace(fit_full_fishing, pars=c("phi_mu","phi_sigma","gamma_mu","gamma_sigma","tau"), ncol=2)
trace_growth_hyper
ggsave(file="./RevisedFigures/Growth_hyper_Traceplot.png", plot=trace_growth_hyper, width=12, height=10)

#Condition omega1 and zeta
trace_cond_omega1 <- stan_trace(fit_full_fishing, pars=c("omega1", "zeta"))
trace_cond_omega1
ggsave(file="./RevisedFigures/Cond_omega1_Traceplot.png", plot=trace_cond_omega1, width=12, height=10)

#Condition omega2
trace_cond_omega2 <- stan_trace(fit_full_fishing, pars=c("omega2"))
trace_cond_omega2
ggsave(file="./RevisedFigures/Cond_omega2_Traceplot.png", plot=trace_cond_omega1, width=12, height=10)

#Condition omega3
trace_cond_omega3 <- stan_trace(fit_full_fishing, pars=c("omega3"))
trace_cond_omega3
ggsave(file="./RevisedFigures/Cond_omega3_Traceplot.png", plot=trace_cond_omega1, width=12, height=10)

#Condition omega4
trace_cond_omega4 <- stan_trace(fit_full_fishing, pars=c("omega4"))
trace_cond_omega4
ggsave(file="./RevisedFigures/Cond_omega4_Traceplot.png", plot=trace_cond_omega1, width=12, height=10)


########
######Extract posteriors to calculate PMRN midpoints########
##Posterior draws of probability of first maturation m
m <- rstan::extract(fit_full_fishing,pars='m')$m
dim(m)

#Create dataset with length, age and cohort for each m
w <- bind_cols(robustbase::colMedians(m),TL.mm,Age, Cohort)
head(w)
names(w) <- c('m','TL.mm','Age','Cohort')
head(w)

#Examine instances of negative m - fish that were already mature, poor model fits
table(filter(w, m<0)$Cohort)
neg.cohorts <- names(table(filter(w, m<0)$Cohort))

#View cohorts with individuals with median mat probabilities that are negative
ggplot(filter(w, Cohort %in% neg.cohorts), aes(x=TL.mm, y=m, color=as.factor(Age))) + 
  geom_point() + facet_wrap(~Cohort, scales="free")

#Calculate PMRN midpoints, 25%, and 75% using logistic regression
##NA if < 2 datapoints, <2 unique lengths, or slope of logistic is negative
#Function to run logistic regression and return 25, 50, and 75% or NAs
midpoints <- function(m, TL.mm) {
  m[m<0] <- 0
  m <- as.numeric(m)
  
  if (length(TL.mm)<2 | n_distinct(TL.mm)==1) {
    Lp50 <- NA
    Lp75 <- NA
    Lp25 <- NA
  } else {
    
    reg <- suppressWarnings(glm(m ~ TL.mm, family=binomial))
    if (coef(reg)[2] < 0) { 
      Lp50 <- NA
      Lp75 <- NA 
      Lp25 <- NA 
      
    } else {
      Lp50 <- -coef(reg)[1]/coef(reg)[2]
      Lp75 <- (log(1/.75-1) + coef(reg)[1])/-coef(reg)[2]
      Lp25 <- (log(1/.25-1) + coef(reg)[1])/-coef(reg)[2]
    }
  }
  #pb$tick()  #theoretically able to use this to track progress in midpoints pipe
  #but kept throwing error so commented out for now
  return(c(Lp50, Lp25, Lp75))
}

#Create data.frame of posterior draws of m for each fish, length, age, and cohort
mats <- data.frame(bind_cols(t(m), TL.mm, Age, Cohort))
dim(mats)
#head(mats)

ms <- paste0("m",1:dim(m)[1]) #external vector of names for obs of m
names(mats) <- c(paste0("m",1:dim(m)[1]), "TL.mm","Age","Cohort") #rename dataframe

timestart=Sys.time() #check time to run
#New and improved pipe to estimate and summarize midpoints using logistic regression
#Takes approx 9 minutes for 3000 posterior draws
Lps <- mats %>%
  filter(Age %in% c(2:5)) %>%
  group_by(Cohort, Age) %>%
  summarize(across(all_of(ms), ~midpoints(m=.x, TL.mm=TL.mm))) %>%
  mutate(Lp = rep(c("Lp50","Lp25",'Lp75'), n_distinct(Cohort, Age))) %>%
  ungroup(.) %>%
  mutate(mean=rowMeans(x=select(., all_of(ms)), na.rm=T),
         median = apply(select(.,all_of(ms)), 1, median, na.rm=T),
         CI2.5 = apply(select(., all_of(ms)), 1, quantile, probs=0.025, na.rm=T),
         CI97.5 = apply(select(., all_of(ms)), 1, quantile, probs=0.975, na.rm=T),
         NAs = rowSums(is.na(select(., all_of(ms))))) %>%
  select(Cohort, Age, Lp, mean, median, CI2.5, CI97.5, NAs)
Sys.time() - timestart

#Match up years and combine with Lp estimates
yearnames <- unique(select(female_yep, CohortYear, Cohort))
yearnames
Lps <- left_join(Lps, yearnames)
Lps

#Filter out uncertain Lps (more than 5% NA) or biologically unrealistic values (median > 600 mm TL)
filtdat <- filter(Lps, Lp=="Lp50",NAs < length(ms)*0.05, median < 600) %>%
  group_by(Age) %>%
  mutate(xjitter = sort(rnorm(n=n(), mean=0, sd=0.1)))

#Create a full dataset for plotting and visualization
all.age.cohort <- filtdat %>% expand(Age, Cohort) %>% left_join(yearnames) %>% filter(Age %in% c(2:5))
filtdat2 <- left_join(all.age.cohort, filtdat) %>% arrange(Cohort, Age) %>%
  group_by(Age) %>%
  mutate(xjitter = sort(rnorm(n=n(), mean=0, sd=0.1)))

filtdat2

#Plot Lp50s with no more than 5% NAs
timeplot <- ggplot(filtdat2, aes(x=CohortYear, y=median, color=factor(Age))) + 
  geom_point(size=2) + geom_line(size=1.25) + 
  scale_y_continuous(limits=c(0,500), oob=scales::squish) + 
  geom_errorbar(aes(x=CohortYear, ymin=CI2.5, ymax=CI97.5), width=0.0) + 
  theme_classic() + scale_color_viridis_d() + 
  labs(x="Cohort year",y=bquote(Lp[50]), color="Age") + scale_x_continuous(breaks=c(seq(1980,2020,5)))

timeplot

ageplot <- ggplot(filtdat2, aes(x=Age+xjitter, y=median, group=CohortYear, color=CohortYear)) + 
  geom_line(lwd=1) + geom_point() +
  geom_errorbar(aes(x=Age+xjitter, ymin=CI2.5, ymax=CI97.5, group=CohortYear), width=0.00) +
  scale_color_viridis_c() + theme_classic() + scale_y_continuous(limits=c(0,500), oob=scales::squish)+
  labs(x="Age", y=bquote(Lp[50]), color="Cohort year") 
ageplot

library(gridExtra)
library(ggpubr)
finalplot <- ggarrange(ageplot, timeplot, ncol=1, labels=c("a)","b)"), font.label=list(face="plain"),
                       hjust=-3.8)
finalplot
ggsave("~/External Projects/Lake Michigan YEP FIE/Perch-FIE/Figures/PMRN_plot.png", 
       plot=finalplot, width=3.5, height=5, units="in", dpi=500,
       scale=1.5)

#save.image("YEP_PMRNrun_Revised_7.26.2024.Rdata")

########Summary table information#########
###########Sample size tables
female_yep

Ntable <- female_yep %>%
  group_by(CohortYear, AGE) %>%
  summarize(N=n(), n_imm = sum(MAT =="I"), n_mat=sum(MAT != "I"),
            meanTL = mean(LENGTH), meanWr = mean(Wr, na.rm=T))
Ntable
print(Ntable, n=Inf)


NCohorttable <- female_yep %>%
  group_by(CohortYear) %>%
  summarize(N=n(), n_imm = sum(MAT =="I"), n_mat=sum(MAT != "I"),
            meanTL = mean(LENGTH), meanWr = mean(Wr, na.rm=T))
NCohorttable
print(NCohorttable, n=Inf)

#Plot proportion mature by cohort
cohortplot <- NCohorttable %>%
  pivot_longer(cols=c(n_imm, n_mat), names_to="Maturity")
cohortplot
cohortmatplot <- ggplot(cohortplot, aes(x=CohortYear, y=value, fill=Maturity)) + 
  geom_bar(stat="identity") + theme_minimal() + theme(axis.line=element_line(color="black")) + 
  scale_fill_discrete(label=c("Immature","Mature")) + scale_y_continuous(breaks=seq(0,1000,100)) + 
  ylab("Number") + xlab("Cohort")

#Create supplemental figure
ggsave(file="./RevisedFigures/CohortMaturityPlot.png", plot=cohortmatplot, width=8, height=6)






