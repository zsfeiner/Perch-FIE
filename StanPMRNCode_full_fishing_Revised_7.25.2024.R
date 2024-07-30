library(rstan)
library(coda)
library(tidyverse)
library(ggplot2)

#Run TP_Data.R to load TP data - will be ignored in model
source("TP_Data.R")
#Run SummarizeWaterTemp.R to load water temp data
source("SummarizeWaterTemps.R")

female_yep <- read.csv("female_yep.csv", 
                       colClasses=c("integer","character","integer","integer","character","character","integer","integer","integer","integer"))

#assign cohort year
female_yep$CohortYear=female_yep$YEAR-female_yep$AGE


#assign cohort number
female_yep$Cohort = female_yep$CohortYear - min(female_yep$CohortYear) + 1

female_meanTL = mean(female_yep$LENGTH)

female_SDTL = sd(female_yep$LENGTH)

#####Do some comparisons of linear vs log-linear model fit#####
plot(female_yep$LENGTH~female_yep$AGE)
plot(LENGTH ~ AGE, data=filter(female_yep, AGE<=6))

mod1 <- lm(LENGTH ~ AGE, data=filter(female_yep, AGE<=6))
mod2 <- lm(LENGTH ~ log(AGE), data=filter(female_yep, AGE<=6))
mod3 <- nls(LENGTH ~ exp(phi + gamma * log(AGE)), data=filter(female_yep, AGE <= 6),
            start=list(phi=1, gamma=0))
mod4 <- nls(LENGTH ~ Linf * (1-exp(-k*(AGE-t0))), data=filter(female_yep, AGE <= 6),
            start=list(Linf=300, k=0.2, t0=0))
AIC(mod1, mod2, mod3, mod4)

plot(LENGTH ~ AGE, data=filter(female_yep, AGE<=6))
lines(predict(mod1, newdata=list(AGE=c(0,1,2,3,4,5,6))) ~ c(0:6), col="red")
lines(predict(mod2, newdata=list(AGE=c(0,1,2,3,4,5,6))) ~ c(0:6), col="green")
lines(predict(mod3, newdata=list(AGE=c(0,1,2,3,4,5,6))) ~ c(0:6), col="blue")
lines(predict(mod4, newdata=list(AGE=c(0,1,2,3,4,5,6))) ~ c(0:6), col="purple")

#Use log(age) model


#Length/Weight
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

#Add relative weight column
female_yep$Ws <- 10^(-5.386 + 3.230 * log10(female_yep$LENGTH))
female_yep$Wr <- female_yep$WEIGHT / female_yep$Ws * 100

#Add environmental variables - use mean TP and annual GDD5
mean_TP_byYear
temps.summary

#Combine TP, temp
GDD_byYear <- temps.summary %>%
  select("GDD5Annual","Year")

#Not used
#TP_byYear <- mean_TP_byYear %>%
# select("Year","Mean") %>%
# filter(Year > 1982 & Year < 2016)

#Combine abiotics and fishing indicator
#I think GDD should be matched up by year of capture, not cohort year
#Scaling GDD and abiotics GDD is also probably not the right thing to do - 
##it has them on different scales

#Going to rethink this and have GDD by year instead
#Start abiotics one year before last cohortyear

abiotics <- GDD_byYear %>%
  inner_join(mean_TP_byYear) %>%
  mutate(CohortYear = Year - min(female_yep$CohortYear) + 1) %>%
  filter(Year >= min(female_yep$YEAR) - 1) %>%
  rename(TP=Mean) %>%
  mutate(Fishing = ifelse(Year <= 1996, 1, 0)) %>%
  select(Year, CohortYear, GDD5Annual, TP, Fishing)

#likely going to have to feed fishing and GDD in separately as GDD should index on year
##and fishing should index on cohort - or not, should it be on year?
#Fishing <- tibble(CohortYear=seq(min(female_yep$CohortYear)-1, max(female_yep$CohortYear))) %>%
#  mutate(Cohort=seq(0,n()-1), Fishing = ifelse(CohortYear <= 1996, 1, 0))


female_yep <- left_join(female_yep, select(temps.summary, Year, GDD5Annual), by=c("YEAR"="Year")) %>%
  rename("GDD5"="GDD5Annual") %>%
  left_join(select(mean_TP_byYear, Year, Mean), by=c("YEAR"="Year")) %>%
  rename("TP"="Mean") %>%
  mutate(Fishing = ifelse(YEAR <= 1996, 1, 0))

######FOR NOW FILTER TO ENVIRONMENTAL DATA######
#Not needed - no GDD NAs
#female_yep <- filter(female_yep, !is.na(mean_TP), !is.na(GDD5))
#female_yep

#######Reassign cohorts#####Not Needed
#assign cohort year
#female_yep$CohortYear=female_yep$YEAR-female_yep$AGE

#assign cohort number
#female_yep$Cohort = female_yep$CohortYear - min(female_yep$CohortYear) +1
##########################





#Female YEP
N = nrow(female_yep)

#Explanatory variables
Age <- (female_yep$AGE)
TL.mm <- female_yep$LENGTH
mean_TL <- mean(female_yep$LENGTH)
sd_TL <- sd(female_yep$LENGTH)
RW <- female_yep$Wr
mean_RW <- mean(female_yep$Wr)
sd_RW <- sd(female_yep$Wr)
TP <- female_yep$TP
mean_TP <- mean(abiotics$TP)
sd_TP <- sd(abiotics$TP)
Fished <- female_yep$Fishing
GDD5 <- female_yep$GDD5
mean_GDD5 <- mean(abiotics$GDD5Annual) #I think this is how this should be scaled
sd_GDD5 <- sd(abiotics$GDD5Annual) #I think this is how this should be scaled
Cohort <- female_yep$Cohort
Year <- female_yep$YEAR - min(female_yep$YEAR) + 1  #Year marker starting from 1 to get prev year's GDD from abiotics
nCohorts <- length(unique(female_yep$Cohort))
X <- cbind(Age,TL.mm,RW,GDD5,TP,Fished)
K <- ncol(X)
Mat <- ifelse(female_yep$MAT=="I",0,1)
nAbiotics <- nrow(abiotics)
#nFishing <- nrow(Fishing)


#Create datalist
dat <- list('N'=nrow(female_yep), 'K'=K, 'Mat'=Mat, 'Age'=Age, 'TL'=TL.mm, 
            'Cohort'=Cohort, 'nCohorts'=nCohorts, 'mean_TL'=mean_TL, 'sd_TL'=sd_TL,
            'RW'=RW, 'mean_RW'=mean_RW, 'sd_RW'=sd_RW,
            'GDD5'=GDD5, 'mean_GDD5'=mean_GDD5, 'sd_GDD5'=sd_GDD5,
            'TP'=TP, 'mean_TP'=mean_TP, 'sd_TP'=sd_TP,
            'abiotics'=abiotics, 'Year'=Year, 'nAbiotics'=nAbiotics,
            'Fished'=Fished)


inits_full_fishing <- function() {
  beta <- c(rnorm(1,-3,0.05),rnorm(1,1,0.05),rnorm(1,2,0.05),rnorm(1,0,0.05),rnorm(1,0,0.05),rnorm(1,0,0.05),rnorm(1,0,0.05),rnorm(1,0,0.05),rnorm(1,0,0.05))
  sigma_u <- c(runif(1,1,5),runif(1,0.1,1),runif(1,0.5,2),runif(1,0.1,1),runif(1,0.1,1),runif(1,0.1,1),runif(1,0.1,1),runif(1,0.1,1))
  phi_mu <- rnorm(1,100,0.05)
  phi_sigma <- runif(1,0.5,1)
  gamma_mu <- rnorm(1,100,0.05)
  gamma_sigma <- runif(1,0.5,2)
  sigma <- runif(1,1,10)
  z_u <- matrix(rnorm(8*nCohorts,0,0.5),nrow=8,ncol=nCohorts)
  init <- list("beta"=beta, "sigma_u"=sigma_u, "phi_mu"=phi_mu, "phi_sigma"=phi_sigma, "gamma_mu"=gamma_mu, "gamma_sigma"=gamma_sigma, "sigma"=sigma, "z_u"=z_u)
  return(init)
}


stanmatcode_full_fishing = stan_model(file = 'yep_fie_covar_revised_7.25.2024.stan')
fit_full_fishing = sampling(stanmatcode_full_fishing, data=dat, init=inits_full_fishing, 
                    iter=4000, warmup=2000, thin=1, chains=3, cores=3, #was 4000 and 2000
                    control=list(adapt_delta=0.90,max_treedepth=10) )
saveRDS(fit_full_fishing,"YEPFIE_covar_enviro_revised_7.25.2024.RDS")

print(fit_full_fishing, pars=c('beta','sigma_u','phi_mu','gamma_mu','sigma'), digits=3, prob=c(0.025,0.5,0.975))
print(fit_full_fishing, pars=c('beta'), digits=3, prob=c(0.025,0.5,0.975))

print(fit_full_fishing, pars=c('m'), digits=3, prob=c(0.025,0.5,0.975))
print(fit_full_fishing, pars=c('s'), digits=3, prob=c(0.025,0.5,0.975))

stan_trace(fit_full_fishing,pars=c('p[1]','p[2]','p[3]','p[4]','p[5]','p[6]'))
stan_trace(fit_full_fishing,pars=c('prev_p[1]','prev_p[2]','prev_p[3]','prev_p[4]','prev_p[5]','prev_p[6]'))
stan_trace(fit_full_fishing,pars=c('s[1]','s[2]','s[3]','s[4]','s[5]','s[6]'))
stan_trace(fit_full_fishing,pars=c('m[1]','m[2]','m[3]','m[4]','m[5]','m[6]'))
stan_trace(fit_full_fishing,pars=c('beta','gamma_mu','phi_mu'))
stan_trace(fit_full_fishing,pars=c('L_u'))
pairs(fit_full_fishing,pars=c('beta','gamma_mu','phi_mu'))

print(fit_full_fishing,pars=c('p[1251]','prev_p[1251]'))
print(fit_full_fishing,pars=c('m[1251]'))


trace_1<-stan_trace(fit_full_fishing,pars=c('beta','sigma_u','phi_mu','phi_sigma','gamma_mu','gamma_sigma','sigma'))
ggsave(file="./Figures/App_traceplot.svg", plot=trace_1, width=16, height=10)
ggsave(file="./Figures/App_traceplot.png", plot=trace_1, width=16, height=10)



m <- rstan::extract(fit_full_fishing,pars='m')$m
dim(m)

w <- bind_cols(robustbase::colMedians(m),TL.mm,Age, Cohort)
head(w)
names(w) <- c('m','TL.mm','Age','Cohort')
head(w)

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

#library(progress) #attempted progress bar, couldn't figure it out
#pb <- progress_bar$new(format = "[:bar] :current/:total (:percent) elapsed :elapsed eta :eta",
#                            total = n_distinct(mats$Cohort, mats$Age))

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

#Match up years
yearnames <- unique(select(female_yep, CohortYear, Cohort))
yearnames
Lps <- left_join(Lps, yearnames)
Lps

filtdat <- filter(Lps, Lp=="Lp50",NAs < length(ms)*0.05, median < 600) %>%
  group_by(Age) %>%
  mutate(xjitter = sort(rnorm(n=n(), mean=0, sd=0.1)))

all.age.cohort <- filtdat %>% expand(Age, Cohort) %>% left_join(yearnames) %>% filter(Age %in% c(2:5))
filtdat2 <- left_join(all.age.cohort, filtdat) %>% arrange(Cohort, Age) %>%
  group_by(Age) %>%
  mutate(xjitter = sort(rnorm(n=n(), mean=0, sd=0.1)))

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


save.image("YEP_PMRNrun_Revised_7.26.2024.Rdata")


#Examine variance of Lp estimates over time
mats
head(mats)
dim(mats)

CV <- function(x) { return(sd(x, na.rm=T)/mean(x, na.rm=T))}

makezero <- function(x) {x[x<0] <- 0; return(x)}

matCVs <- mats %>%
  filter(Age %in% c(2:5)) %>%
  mutate_all(makezero) %>%
  group_by(Cohort, Age) %>%
  summarize_at(vars(ms), var) %>%
  ungroup(.) %>%
  mutate(mean=rowMeans(x=select(., all_of(ms)), na.rm=T),
         median = apply(select(.,all_of(ms)), 1, median, na.rm=T),
         CI2.5 = apply(select(., all_of(ms)), 1, quantile, probs=0.025, na.rm=T),
         CI97.5 = apply(select(., all_of(ms)), 1, quantile, probs=0.975, na.rm=T),
         NAs = rowSums(is.na(select(., all_of(ms))))) %>%
  select(Cohort, Age, mean, median, CI2.5, CI97.5, NAs)

filtCVs <- filter(matCVs, NAs < length(ms)*0.05) %>%
  group_by(Age) %>%
  mutate(xjitter = sort(rnorm(n=n(), mean=0, sd=0.1)))

filtCVs <- left_join(filtCVs, yearnames)

all.age.cohort <- filtCVs %>% expand(Age, Cohort) %>% left_join(yearnames) %>% filter(Age %in% c(2:5))
filtCVs2 <- left_join(all.age.cohort, filtCVs) %>% arrange(Cohort, Age) %>%
  group_by(Age) %>%
  mutate(xjitter = sort(rnorm(n=n(), mean=0, sd=0.1)))

#Plot CVs with no more than 5% NAs
timeCVs <- ggplot(filtCVs2, aes(x=CohortYear, y=median, color=factor(Age))) + 
  geom_point(size=2) + geom_line(size=1.25) + 
  #scale_y_continuous(lim, oob=scales::squish) + 
  geom_errorbar(aes(x=CohortYear, ymin=CI2.5, ymax=CI97.5), width=0.0) + 
  theme_classic() + scale_color_viridis_d() + geom_smooth(se=F) + 
  labs(x="Cohort year",y="Among-individual maturation probability CV", color="Age") + scale_x_continuous(breaks=c(seq(1980,2020,5)))
timeCVs

ageCVs <- ggplot(filtCVs2, aes(x=Age+xjitter, y=median, group=CohortYear, color=CohortYear)) + 
  #geom_line(lwd=1) + 
  geom_point() +
  geom_errorbar(aes(x=Age+xjitter, ymin=CI2.5, ymax=CI97.5, group=CohortYear), width=0.00) +
  scale_color_viridis_c() + theme_classic() + #scale_y_continuous(limits=c(0,500), oob=scales::squish)+
  labs(x="Age", y="Among-individual maturation probability CV", color="Cohort year") 
ageCVs

CVplot <- ggarrange(ageCVs, timeCVs, ncol=1, labels=c("a)","b)"), font.label=list(face="plain"),
                       hjust=-3.8)
CVplot
ggsave("~/External Projects/Lake Michigan YEP FIE/Perch-FIE/Figures/PMRNSD_plot.png", 
       plot=CVplot, width=3.5, height=5, units="in", dpi=500,
       scale=1.5)



#Sample size tables
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

cohortplot <- NCohorttable %>%
  pivot_longer(cols=c(n_imm, n_mat), names_to="Maturity")
cohortplot
ggplot(cohortplot, aes(x=CohortYear, y=value, fill=Maturity)) + 
  geom_bar(stat="identity") + theme_minimal() + theme(axis.line=element_line(color="black")) + 
  scale_fill_discrete(label=c("Immature","Mature")) + scale_y_continuous(breaks=seq(0,1000,100))


dim(mats)
robustbase::colMedians(mats)
new.w <- w %>%
  filter(Age %in% c(2)) %>%
  mutate(m=ifelse(m<0,0,m))
ggplot(data=new.w, aes(x=m, group=Cohort, color=Cohort, after_stat(scaled))) + 
  geom_density() + facet_wrap(~Age * Cohort)






