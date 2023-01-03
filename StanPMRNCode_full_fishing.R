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

plot(female_yep$LENGTH~female_yep$AGE)

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
  select(Year, CohortYear, GDD5Annual, TP)

#likely going to have to feed fishing and GDD in separately as GDD should index on year
##and fishing should index on cohort - or not, should it be on year?

Fishing <- tibble(CohortYear=seq(min(female_yep$CohortYear)-1, max(female_yep$CohortYear))) %>%
  mutate(Cohort=seq(0,n()-1), Fishing = ifelse(CohortYear <= 1996, 1, 0))


female_yep <- left_join(female_yep, select(temps.summary, Year, GDD5Annual), by=c("YEAR"="Year")) %>%
  rename("GDD5"="GDD5Annual") %>%
  left_join(select(mean_TP_byYear, Year, Mean), by=c("YEAR"="Year")) %>%
  rename("TP"="Mean") %>%
  mutate(Fishing = ifelse(female_yep$CohortYear <= 1996, 1, 0))

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
nFishing <- nrow(Fishing)


#Create datalist
dat <- list('N'=nrow(female_yep), 'K'=K, 'Mat'=Mat, 'Age'=Age, 'TL'=TL.mm, 
            'Cohort'=Cohort, 'nCohorts'=nCohorts, 'mean_TL'=mean_TL, 'sd_TL'=sd_TL,
            'RW'=RW, 'mean_RW'=mean_RW, 'sd_RW'=sd_RW,
            'GDD5'=GDD5, 'mean_GDD5'=mean_GDD5, 'sd_GDD5'=sd_GDD5,
            'TP'=TP, 'mean_TP'=mean_TP, 'sd_TP'=sd_TP,
            'abiotics'=abiotics, 'Year'=Year, 'nAbiotics'=nAbiotics,
            'fishing'=Fishing, 'nFishing'=nFishing, 'Fished'=Fished)


inits_full_fishing <- function() {
  beta <- c(rnorm(1,-3,0.05),rnorm(1,1,0.05),rnorm(1,2,0.05),rnorm(1,0,0.05),rnorm(1,0,0.05),rnorm(1,0,0.05),rnorm(1,0,0.05),rnorm(1,0,0.05),rnorm(1,0,0.05))
  sigma_u <- c(runif(1,1,5),runif(1,0.1,1),runif(1,0.5,2),runif(1,0.1,1),runif(1,0.1,1),runif(1,0.1,1),runif(1,0.1,1),runif(1,0.1,1),runif(1,0.1,1))
  phi_mu <- rnorm(1,200,0.05)
  phi_sigma <- runif(1,0.5,1)
  gamma_mu <- rnorm(1,40,0.05)
  gamma_sigma <- runif(1,0.5,2)
  sigma <- runif(1,1,10)
  z_u <- matrix(rnorm(9*nCohorts,0,0.5),nrow=9,ncol=nCohorts)
  init <- list("beta"=beta, "sigma_u"=sigma_u, "phi_mu"=phi_mu, "phi_sigma"=phi_sigma, "gamma_mu"=gamma_mu, "gamma_sigma"=gamma_sigma, "sigma"=sigma, "z_u"=z_u)
  return(init)
}


stanmatcode_full_fishing = stan_model(file = 'yep_fie_covar_full_fishing.stan')
fit_full_fishing = sampling(stanmatcode_full_fishing, data=dat, init=inits_full_fishing, 
                    iter=4000, warmup=2000, thin=1, chains=3, cores=3, #was 4000 and 2000
                    control=list(adapt_delta=0.90,max_treedepth=10) )
saveRDS(fit_full_fishing,"YEPFIE_covar_enviro_full_fishing.RDS")

print(fit_full_fishing, pars=c('beta','sigma_u','phi_mu','gamma_mu','sigma'), digits=3, prob=c(0.025,0.5,0.975))
print(fit_full_fishing, pars=c('m'), digits=3, prob=c(0.025,0.5,0.975))
stan_trace(fit_full_fishing,pars=c('p[1]','p[2]','p[3]','p[4]','p[5]','p[6]'))
stan_trace(fit_full_fishing,pars=c('m[100]','m[200]','m[300]','m[400]','m[500]','m[600]'))
stan_trace(fit_full_fishing,pars=c('beta','gamma_mu','phi_mu'))
stan_trace(fit_full_fishing,pars=c('L_u'))
pairs(fit_full_fishing,pars=c('beta','gamma_mu','phi_mu'))

print(fit_full_fishing,pars=c('p[1251]','prev_p[1251]'))
print(fit_full_fishing,pars=c('m[1251]'))



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

#Plot Lp50s with no more than 5% NAs
ggplot(filter(Lps, Lp=="Lp50", NAs < length(ms)*0.05), aes(x=CohortYear, y=median, color=factor(Age))) + 
  geom_point(size=2) + geom_line(size=1.25) + ylim(0,500) +
  geom_errorbar(aes(x=CohortYear, ymin=CI2.5, ymax=CI97.5), width=0.05) + 
  theme_classic()

save.image("YEP_PMRNrun_full.Rdata")
