library(rstan)
library(coda)

#Run TP_Data.R to load TP data - will be ignored in model
#source("TP_Data.R")
#Run SummarizeWaterTemp.R to load water temp data
source("SummarizeWaterTemps.R")

#Read Data
#female_yep <- read.csv("~/Research/LMYEP_Genetics/Perch-FIE/female_yep.csv",
#                       colClasses=c("integer","character","integer","integer","character","character","integer","integer","integer","integer"))
#male_yep <- read.csv("~/Research/LMYEP_Genetics/Perch-FIE/male_yep.csv")

female_yep <- read.csv("female_yep.csv", 
                       colClasses=c("integer","character","integer","integer","character","character","integer","integer","integer","integer"))

#assign cohort year
female_yep$CohortYear=female_yep$YEAR-female_yep$AGE
#male_yep$CohortYear=male_yep$YEAR-male_yep$AGE

#assign cohort number
female_yep$Cohort = female_yep$CohortYear - min(female_yep$CohortYear) + 1
#male_yep$Cohort = male_yep$CohortYear - min(male_yep$CohortYear) +1

female_meanTL = mean(female_yep$LENGTH)
#male_meanTL = mean(male_yep$LENGTH)

female_SDTL = sd(female_yep$LENGTH)
#male_SDTL = sd(male_yep$LENGTH)

plot(female_yep$LENGTH~female_yep$AGE)
#plot(male_yep$LENGTH~male_yep$AGE)


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
abiotics <- GDD_byYear
abiotics$CohortYear = abiotics$Year-1982
abiotics$Fishing = ifelse(abiotics$Year <= 1996, 1, 0)
abiotics <- select(abiotics, -Fishing) #Drop Fishing


female_yep <- left_join(female_yep, select(temps.summary, Year, GDD5Annual), by=c("CohortYear"="Year")) %>%
  rename("GDD5"="GDD5Annual")
female_yep

#Add commercial fishing indicator
female_yep$Fishing = ifelse(female_yep$CohortYear <= 1996, 1, 0)

######FOR NOW FILTER TO ENVIRONMENTAL DATA######
#Not needed - no GDD NAs
#female_yep <- filter(female_yep, !is.na(mean_TP), !is.na(GDD5))
#female_yep

#######Reassign cohorts#####
#assign cohort year
female_yep$CohortYear=female_yep$YEAR-female_yep$AGE

#assign cohort number
female_yep$Cohort = female_yep$CohortYear - min(female_yep$CohortYear) +1
#########################


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
#TP <- female_yep$mean_TP
#mean_TP <- mean(female_yep$mean_TP)
#sd_TP <- sd(female_yep$mean_TP)
#Fishing <- female_yep$Fishing
GDD5 <- female_yep$GDD5
mean_GDD5 <- mean(female_yep$GDD5)
sd_GDD5 <- sd(female_yep$GDD5)
Cohort <- female_yep$Cohort
nCohorts <- length(unique(female_yep$Cohort))
X <- cbind(Age,TL.mm,RW,GDD5)
K <- ncol(X)
Mat <- ifelse(female_yep$MAT=="I",0,1)


#Create datalist
dat <- list('N'=nrow(female_yep), 'K'=K, 'Mat'=Mat, 'Age'=Age, 'TL'=TL.mm, 
            'Cohort'=Cohort, 'nCohorts'=nCohorts, 'mean_TL'=mean_TL, 'sd_TL'=sd_TL,
            'RW'=RW, 'mean_RW'=mean_RW, 'sd_RW'=sd_RW,
            'GDD5'=GDD5, 'mean_GDD5'=mean_GDD5, 'sd_GDD5'=sd_GDD5,
            'abiotics'=abiotics)


inits_nofish <- function() {
  beta <- c(rnorm(1,-3,0.05),rnorm(1,1,0.05),rnorm(1,2,0.05),rnorm(1,0,0.05),rnorm(1,0,0.05),rnorm(1,0,0.05),rnorm(1,0,0.05))
  sigma_u <- c(runif(1,1,5),runif(1,0.1,1),runif(1,0.5,2),runif(1,0.1,1),runif(1,0.1,1),runif(1,0.1,1),runif(1,0.1,1))
  phi_mu <- rnorm(1,200,0.05)
  phi_sigma <- runif(1,0.5,1)
  gamma_mu <- rnorm(1,40,0.05)
  gamma_sigma <- runif(1,0.5,2)
  sigma <- runif(1,1,10)
  z_u <- matrix(rnorm(7*nCohorts,0,0.5),nrow=7,ncol=nCohorts)
  init <- list("beta"=beta, "sigma_u"=sigma_u, "phi_mu"=phi_mu, "phi_sigma"=phi_sigma, "gamma_mu"=gamma_mu, "gamma_sigma"=gamma_sigma, "sigma"=sigma, "z_u"=z_u)
  return(init)
}


stanmatcode_nofish = stan_model(file = 'yep_fie_covar_NoFishing.stan')
fit_nofish = sampling(stanmatcode_nofish, data=dat, init=inits_nofish, 
               iter=4000, warmup=2000, thin=1, chains=3, cores=3, #was 4000 and 2000
               control=list(adapt_delta=0.90,max_treedepth=10) )
saveRDS(fit_nofish,"YEPFIE_covar_enviro_nofish.RDS")

print(fit_nofish, pars=c('beta','sigma_u','phi_mu','gamma_mu','sigma'), digits=3, prob=c(0.025,0.5,0.975))
print(fit_nofish, pars=c('m'), digits=3, prob=c(0.025,0.5,0.975))
stan_trace(fit_nofish,pars=c('p[1]','p[2]','p[3]','p[4]','p[5]','p[6]'))
stan_trace(fit_nofish,pars=c('m[100]','m[200]','m[300]','m[400]','m[500]','m[600]'))
stan_trace(fit_nofish,pars=c('beta','gamma_mu','phi_mu'))
stan_trace(fit_nofish,pars=c('L_u'))
pairs(fit_nofish,pars=c('beta','gamma_mu','phi_mu'))

print(fit_nofish,pars=c('p[1251]','prev_p[1251]'))
print(fit_nofish,pars=c('m[1251]'))


#plot(female_yep$TL.mm[female_yep$Cohort==9],female_yep$Mat[female_yep$Cohort==9],col=female_yep$Age[female_yep$Cohort==9])
## Use the L matrices to compute the correlation matrices
# L matrices
#L_u <- extract(fit, pars = "L_u")$L_u

# correlation parameters
#cor_u <- apply(L_u, 1, function(x) tcrossprod(x)[4, 2])
#print(signif(quantile(cor_u, probs = c(0.025, 0.5, 0.975)), 2))
#print(mean(cor_u))

#<<<<<<< HEAD
m <- rstan::extract(fit_nofish, pars='m')$m
#=======
m <- data.frame(rstan::extract(fit_nofish,pars='m'))
#>>>>>>> ba2c4ca1ef35fef90db7f003f5199a4fcadcbc18
m

#Create dataframe of maturation prob, length, age, and cohort for plotting
w <- bind_cols(as.numeric(m[1,]),TL.mm,Age, Cohort)
names(w) <- c("m", "TL.mm","Age","Cohort")

#Check cohorts with negative maturation probability
table(w[w$m<0,4])

#Plot maturation data for cohorts with negative maturation probability
#In some cohorts very large fish have negative maturation probabilities because
#some are identified as immature - senescence, mis-ID of spent fish?  Worth removing immature fish > ~275 mm as outliers?
ggplot(filter(female_yep, Cohort %in% c(8,11,16,19,20,21,22,24,26,28)), aes(x=LENGTH, y=MAT, color=as.factor(AGE), group=as.factor(AGE))) + 
  geom_point() + 
  facet_wrap(~Cohort, scales="free")

ggplot(filter(w, Cohort %in% c(8,11,16,19,20,21,22,24,26,28)), aes(x=TL.mm, y=m, group=as.factor(Age), color=as.factor(Age))) + 
  geom_point() + facet_wrap(~Cohort, scales='free')


#Use spline interpolation to determine PMRN midpoints for ages 2-5 (older ages declining in m, already mature)
w
nAges <- length(unique(w$Age[w$Age %in% c(2:5)]))
nCohorts <- length(unique(w$Cohort[w$Age %in% c(2:5)]))
nReps <- dim(m)[1]
Ages <- sort(unique(w$Age[w$Age %in% c(2:5)]))
Cohorts <- sort(unique(w$Cohort[w$Age %in% c(2:5)]))

mat.all <- as.data.frame(cbind(t(m),TL.mm,Age,Cohort))
head(mat.all)
dim(mat.all)

plot(mat.all[mat.all$Cohort==6,1] ~ TL.mm, data=mat.all[mat.all$Cohort==6,],col=mat.all$Age[mat.all$Cohort==6],ylim=c(-0.5,1),xlim=c(50,350))
# Commented out - takes a long time
# for (i in 2:1000) {
#   points(thing[thing$Cohort==6,i] ~ TL.mm, data=thing[thing$Cohort==6,],col=thing$Age[thing$Cohort==6])
# }

#Estimate PMRN midpoints and quantiles for each age and cohort using all posterior draws
#Make array of i ages, j cohorts, and k replicates
Lp50 <- array(0, dim = c(nAges,nCohorts,nReps))
dim(Lp50)
Lp75 <- array(0, dim = c(nAges,nCohorts,nReps))
dim(Lp75)
Lp25 <- array(0, dim = c(nAges,nCohorts,nReps))
dim(Lp25)
Lp50.all<- array(0, dim = c(nAges,nCohorts,nReps))
dim(Lp50.all)

#Determine slope to make sure positive
####First option - Across all replicates, cohorts, and ages use spline interpolation to determine Lp50, Lp25, and Lp75
###Lp50.all extrapolates midpoint, while Lp50 only uses individuals with maturation prob that crosses 50%
###Could probably make this a function at some point
#Note drop fish with m < 0, i.e., negative maturation probability
#Throws NA if - logistic regression of m ~ TL returns negative slope
#Throws NA if - only 1 fish in age class or only 1 unique length in age class
for (k in 1:nReps) {
  
  d <- as.data.frame(cbind(m[k,],TL.mm,Age,Cohort))
  names(d)[1] <- "m"
  
  for (i in 1:nAges) {
    for (j in 1:nCohorts) {
      sub.d <- subset(d, Age == Ages[i] & Cohort == Cohorts[j] & m > 0)
      
      if (nrow(sub.d)<2 | length(unique(sub.d$TL.mm))==1) {
        Lp50[i,j,k] <- NA
        Lp75[i,j,k] <- NA
        Lp25[i,j,k] <- NA
        Lp50.all[i,j,k] <- NA
      } else {
        
        reg <- glm(m ~ TL.mm, family=binomial, data=sub.d)
        if (coef(reg)[2] < 0 | min(sub.d$m, na.rm=TRUE) > 0.5 | max(sub.d$m, na.rm=TRUE) < 0.5) { 
          Lp50[i,j,k] <- NA } else {
            midpoint <- spline(sub.d$m,sub.d$TL.mm,xout=0.5)
            Lp50[i,j,k] <-  midpoint$y
          }
        if (coef(reg)[2] < 0 | min(sub.d$m, na.rm=TRUE) > 0.75 | max(sub.d$m, na.rm=TRUE) < 0.75) { 
          Lp75[i,j,k] <- NA } else {
            midpoint <- spline(sub.d$m,sub.d$TL.mm,xout=0.75)
            Lp75[i,j,k] <-  midpoint$y
          }
        if (coef(reg)[2] < 0 | min(sub.d$m, na.rm=TRUE) > 0.25 | max(sub.d$m, na.rm=TRUE) < 0.25) { 
          Lp25[i,j,k] <- NA } else {
            midpoint <- spline(sub.d$m,sub.d$TL.mm,xout=0.25)
            Lp25[i,j,k] <-  midpoint$y
          }
        if (coef(reg)[2] < 0) { 
          Lp50.all[i,j,k] <- NA } else {
            midpoint <- spline(sub.d$m,sub.d$TL.mm,xout=0.5)
            Lp50.all[i,j,k] <-  midpoint$y
          }
      }
    }
  }
}


#### ######Second option - Use logistic regression####
# for (k in 1:nReps) {
#   
#   d <- as.data.frame(cbind(m[k,],TL.mm,Age,Cohort))
#   names(d)[1] <- "m"
#   
#   for (i in 1:nAges) {
#     for (j in 1:nCohorts) {
#       sub.d <- subset(d, Age == Ages[i] & Cohort == Cohorts[j])
#       sub.d$m[sub.d$m<0] <- 0
#       
#       if (nrow(sub.d)<2 | length(unique(sub.d$TL.mm))==1) {
#         Lp50[i,j,k] <- NA
#         Lp75[i,j,k] <- NA
#         Lp25[i,j,k] <- NA
#         Lp50.all[i,j,k] <- NA
#       } else {
#         
#         reg <- suppressWarnings(glm(m ~ TL.mm, family=binomial, data=sub.d))
#         if (coef(reg)[2] < 0) { 
#           Lp50[i,j,k] <- NA
#           Lp75[i,j,k] <- NA 
#           Lp25[i,j,k] <- NA 
#           
#           } else {
#             Lp50[i,j,k] <- -coef(reg)[1]/coef(reg)[2]
#             Lp75[i,j,k] <- (log(1/.75-1) + coef(reg)[1])/-coef(reg)[2]
#             Lp25[i,j,k]  <- (log(1/.25-1) + coef(reg)[1])/-coef(reg)[2]
#             }
#       }
#     }
#   }
# }
##### 

rowSums(is.na(Lp50[1,,]))
rowSums(is.na(Lp50[2,,]))
rowSums(is.na(Lp50[3,,]))

rowSums(is.na(Lp50.all[1,,]))
rowSums(is.na(Lp50.all[2,,]))
rowSums(is.na(Lp50.all[3,,]))

credible.intervals <- function(x) { 
  creds <- quantile(x, probs=c(0.025,0.5,0.975),na.rm=TRUE)
  return(creds)}

credible.intervals(Lp50.all[3,1,])

apply(Lp50.all[1,,],1,FUN=credible.intervals)
apply(Lp50.all[2,,],1,FUN=credible.intervals)
apply(Lp50.all[3,,],1,FUN=credible.intervals)

PMRNs.Lp50 <- PMRNs.Lp25 <- PMRNs.Lp75 <- PMRNs.Lp50.all <- matrix(0,nrow=0,ncol=7)
colnames(PMRNs.Lp50) <- colnames(PMRNs.Lp25) <- colnames(PMRNs.Lp75) <- colnames(PMRNs.Lp50.all) <- c("Cohort","Age","NAs","Mean","CI2.5","Median","CI97.5")
head(PMRNs.Lp50)
for (i in 1:nCohorts) {
  #Summarize Lp50s
  temp.Lp50 <- matrix(0,nrow=4,ncol=7)
  colnames(temp.Lp50) <- colnames(PMRNs.Lp50)
  z.Lp50 <- Lp50[,i,]
  temp.Lp50[,1] <- i
  temp.Lp50[,2] <- c(2,3,4,5)
  temp.Lp50[,3] <- rowSums(is.na(z.Lp50))
  temp.Lp50[,4] <- rowMeans(z.Lp50,na.rm=TRUE)
  temp.Lp50[,5:7] <- t(apply(z.Lp50,1,FUN=credible.intervals))
  PMRNs.Lp50 <- rbind(PMRNs.Lp50,temp.Lp50)
  
  #Summarize Lp25s
  temp.Lp25 <- matrix(0,nrow=4,ncol=7)
  colnames(temp.Lp25) <- colnames(PMRNs.Lp25)
  z.Lp25 <- Lp25[,i,]
  temp.Lp25[,1] <- i
  temp.Lp25[,2] <- c(2,3,4,5)
  temp.Lp25[,3] <- rowSums(is.na(z.Lp25))
  temp.Lp25[,4] <- rowMeans(z.Lp25,na.rm=TRUE)
  temp.Lp25[,5:7] <- t(apply(z.Lp25,1,FUN=credible.intervals))
  PMRNs.Lp25 <- rbind(PMRNs.Lp25,temp.Lp25)
  
  #Summarize Lp75s
  temp.Lp75 <- matrix(0,nrow=4,ncol=7)
  colnames(temp.Lp75) <- colnames(PMRNs.Lp75)
  z.Lp75 <- Lp75[,i,]
  temp.Lp75[,1] <- i
  temp.Lp75[,2] <- c(2,3,4,5)
  temp.Lp75[,3] <- rowSums(is.na(z.Lp75))
  temp.Lp75[,4] <- rowMeans(z.Lp75,na.rm=TRUE)
  temp.Lp75[,5:7] <- t(apply(z.Lp75,1,FUN=credible.intervals))
  PMRNs.Lp75 <- rbind(PMRNs.Lp75,temp.Lp75)
  
  #Summarize Lp50.alls
  temp.Lp50.all <- matrix(0,nrow=4,ncol=7)
  colnames(temp.Lp50.all) <- colnames(PMRNs.Lp50.all)
  z.Lp50.all <- Lp50.all[,i,]
  temp.Lp50.all[,1] <- i
  temp.Lp50.all[,2] <- c(2,3,4,5)
  temp.Lp50.all[,3] <- rowSums(is.na(z.Lp50.all))
  temp.Lp50.all[,4] <- rowMeans(z.Lp50.all,na.rm=TRUE)
  temp.Lp50.all[,5:7] <- t(apply(z.Lp50.all,1,FUN=credible.intervals))
  PMRNs.Lp50.all <- rbind(PMRNs.Lp50.all,temp.Lp50.all)
  
}

PMRNs.Lp50
PMRNs.Lp75
PMRNs.Lp25
PMRNs.Lp50.all

#Filter out unrealistic or uncertain estimates with too many NAs
PMRNs.Lp50 <- as.data.frame(PMRNs.Lp50)
sub.Lp50 <- subset(PMRNs.Lp50, NAs < 250 & CI2.5 >= -100 & CI97.5 <= 2000)
sub.Lp50

PMRNs.Lp75 <- as.data.frame(PMRNs.Lp75)
sub.Lp75 <- subset(PMRNs.Lp75, NAs < 250 & CI2.5 >= -100 & CI97.5 <= 2000)
sub.Lp75

PMRNs.Lp25 <- as.data.frame(PMRNs.Lp25)
sub.Lp25 <- subset(PMRNs.Lp25, NAs < 250 & CI2.5 >= -100 & CI97.5 <= 2000)
sub.Lp25

PMRNs.Lp50.all <- as.data.frame(PMRNs.Lp50.all)
sub.Lp50.all <- subset(PMRNs.Lp50.all, NAs < 500 & Median<1000 & CI2.5 > 0 & CI97.5 < 1000)
sub.Lp50.all


#Match up years
yearnames <- cbind(Cohort, female_yep$YEAR-female_yep$AGE)
colnames(yearnames) <-c("Cohort","Year")
yearnames <- as.data.frame(yearnames)
sub.Lp50$Year <- yearnames$Year[match(sub.Lp50$Cohort, yearnames$Cohort)]
sub.Lp75$Year <- yearnames$Year[match(sub.Lp75$Cohort, yearnames$Cohort)]
sub.Lp25$Year <- yearnames$Year[match(sub.Lp25$Cohort, yearnames$Cohort)]
sub.Lp50.all$Year <- yearnames$Year[match(sub.Lp50.all$Cohort, yearnames$Cohort)]

sub.Lp50

plot(0,0,ylim=c(100,400),xlim=c(1.5,5.5))
for (i in 1:48) {
  points(sub.Lp50$Age[sub.Lp50$Cohort==i], sub.Lp50$Median[sub.Lp50$Cohort==i],col=i,type="b",pch=20,cex=2)
  
  #lines(x=rep(sub.Lp50$Age[sub.Lp50$Cohort==i],2),y=c(sub.Lp50$CI2.5[sub.Lp50$Cohort==i],sub.Lp50$CI97.5[sub.Lp50$Cohort==i]),col=i)
}

library(ggplot2)

####Everything past this is just plotting the PMRNs in different ways that can be cleaned up later but works for now####
ggplot(data=sub.Lp50,aes(x=Age,y=Median,group=as.factor(Year)))+
  geom_line(aes(color=as.factor(Year)),size=1.10) + theme_classic()

par(mfrow=c(1,3))
plot(Median ~ Year, data=sub.Lp50[sub.Lp50$Age==3,],type="b",col=1,pch=20,ylim=c(100,400),xlim=c(1979,2015),main="Age 3",cex=2,lwd=2,ylab="Median Lp50")
lines(CI2.5 ~ Year, data=sub.Lp50[sub.Lp50$Age==3,],lty=2,col=1)
lines(CI97.5 ~ Year, data=sub.Lp50[sub.Lp50$Age==3,],lty=2,col=1)
abline(reg=lm(Median ~ Year, data=sub.Lp50[sub.Lp50$Age==3,]),col=1)

plot(Median ~ Year, data=sub.Lp50[sub.Lp50$Age==4,],type="b",col=2,pch=20,ylim=c(100,400),xlim=c(1979,2015),main="Age 4",cex=2,lwd=2,ylab="Median Lp50")
lines(CI2.5 ~ Year, data=sub.Lp50[sub.Lp50$Age==4,],lty=2,col=2)
lines(CI97.5 ~ Year, data=sub.Lp50[sub.Lp50$Age==4,],lty=2,col=2)
abline(reg=lm(Median ~ Year, data=sub.Lp50[sub.Lp50$Age==4,]),col=2)

plot(Median ~ Year, data=sub.Lp50[sub.Lp50$Age==5,],type="b",col=3,pch=20,ylim=c(100,400),xlim=c(1979,2015), main="Age 5",cex=2,lwd=2,ylab="Median Lp50")
lines(CI2.5 ~ Year, data=sub.Lp50[sub.Lp50$Age==5,],lty=2,col=3)
lines(CI97.5 ~ Year, data=sub.Lp50[sub.Lp50$Age==5,],lty=2,col=3)
abline(reg=lm(Median ~ Year, data=sub.Lp50[sub.Lp50$Age==5,]),col=3)

par(mfrow=c(1,1))
plot(Median ~ Year, data=sub.Lp50[sub.Lp50$Age==3,],type="b",col=1,pch=20,ylim=c(100,400),xlim=c(1979,2015),main="Age 3",cex=2,lwd=2,ylab="Median Lp50")
points(Median ~ Year, data=sub.Lp50[sub.Lp50$Age==4,],type="b",col=2,pch=20,cex=2,lwd=2)
points(Median ~ Year, data=sub.Lp50[sub.Lp50$Age==5,],type="b",col=3,pch=20,cex=2,lwd=2)
legend(legend=c("Age 3","Age 4",'Age 5'), col=c(1:3),x="topleft",bty="n",lty=1,lwd=2,pch=20,cex=2)
abline(reg=lm(Median ~ Year, data=sub.Lp50[sub.Lp50$Age==3,]),col=1)
abline(reg=lm(Median ~ Year, data=sub.Lp50[sub.Lp50$Age==4,]),col=2)
abline(reg=lm(Median ~ Year, data=sub.Lp50[sub.Lp50$Age==5,]),col=3)

summary(lm(Median ~ Year, data=sub.Lp50[sub.Lp50$Age==3,]))
summary(lm(Median ~ Year, data=sub.Lp50[sub.Lp50$Age==4,]))
summary(lm(Median ~ Year, data=sub.Lp50[sub.Lp50$Age==5,]))


##
par(mfrow=c(1,1))
plot(Median ~ Year, data=sub.Lp75[sub.Lp75$Age==3,],type="b",col=1,pch=20,ylim=c(100,400),xlim=c(1979,2015),main="Age 3",cex=2,lwd=2,ylab="Median Lp75")
points(Median ~ Year, data=sub.Lp75[sub.Lp75$Age==4,],type="b",col=2,pch=20,cex=2,lwd=2)
points(Median ~ Year, data=sub.Lp75[sub.Lp75$Age==5,],type="b",col=3,pch=20,cex=2,lwd=2)
legend(legend=c("Age 3","Age 4",'Age 5'), col=c(1:3),x="topleft",bty="n",lty=1,lwd=2,pch=20,cex=2)
abline(reg=lm(Median ~ Year, data=sub.Lp75[sub.Lp75$Age==3,]),col=1)
abline(reg=lm(Median ~ Year, data=sub.Lp75[sub.Lp75$Age==4,]),col=2)
abline(reg=lm(Median ~ Year, data=sub.Lp75[sub.Lp75$Age==5,]),col=3)

summary(lm(Median ~ Year, data=sub.Lp75[sub.Lp75$Age==3,]))
summary(lm(Median ~ Year, data=sub.Lp75[sub.Lp75$Age==4,]))
summary(lm(Median ~ Year, data=sub.Lp75[sub.Lp75$Age==5,]))


##
par(mfrow=c(1,1))
plot(Median ~ Year, data=sub.Lp25[sub.Lp25$Age==3,],type="b",col=1,pch=20,ylim=c(100,400),xlim=c(1979,2015),main="Age 3",cex=2,lwd=2,ylab="Median Lp25")
points(Median ~ Year, data=sub.Lp25[sub.Lp25$Age==4,],type="b",col=2,pch=20,cex=2,lwd=2)
points(Median ~ Year, data=sub.Lp25[sub.Lp25$Age==5,],type="b",col=3,pch=20,cex=2,lwd=2)
legend(legend=c("Age 3","Age 4",'Age 5'), col=c(1:3),x="topleft",bty="n",lty=1,lwd=2,pch=20,cex=2)
abline(reg=lm(Median ~ Year, data=sub.Lp25[sub.Lp25$Age==3,]),col=1)
abline(reg=lm(Median ~ Year, data=sub.Lp25[sub.Lp25$Age==4,]),col=2)
abline(reg=lm(Median ~ Year, data=sub.Lp25[sub.Lp25$Age==5,]),col=3)

summary(lm(Median ~ Year, data=sub.Lp25[sub.Lp25$Age==3,]))
summary(lm(Median ~ Year, data=sub.Lp25[sub.Lp25$Age==4,]))
summary(lm(Median ~ Year, data=sub.Lp25[sub.Lp25$Age==5,]))

##
par(mfrow=c(1,1))
plot(Median ~ Year, data=sub.Lp50[sub.Lp50$Age==4,],type="b",col=1,pch=20,ylim=c(100,400),xlim=c(1979,2015),main="Age 3",cex=2,lwd=2,ylab="Median Lp25")
points(Median ~ Year, data=sub.Lp25[sub.Lp25$Age==4,],type="b",col=2,pch=20,cex=2,lwd=2)
points(Median ~ Year, data=sub.Lp75[sub.Lp75$Age==4,],type="b",col=3,pch=20,cex=2,lwd=2)

##
par(mfrow=c(1,1))
plot(Median ~ Year, data=sub.Lp50[sub.Lp50$Age==5,],type="b",col=1,pch=20,ylim=c(100,400),xlim=c(1979,2015),main="Age 3",cex=2,lwd=2,ylab="Median Lp25")
points(Median ~ Year, data=sub.Lp25[sub.Lp25$Age==5,],type="b",col=2,pch=20,cex=2,lwd=2)
points(Median ~ Year, data=sub.Lp75[sub.Lp75$Age==5,],type="b",col=3,pch=20,cex=2,lwd=2)


##
par(mfrow=c(1,1))
plot(Median ~ Year, data=sub.Lp50[sub.Lp50$Age==3,],type="b",col=1,pch=20,ylim=c(100,400),xlim=c(1979,2015),main="",cex=2,lwd=2,ylab="Median Lp50")
points(Median ~ Year, data=sub.Lp50[sub.Lp50$Age==4,],type="b",col=2,pch=20,cex=2,lwd=2)
points(Median ~ Year, data=sub.Lp50[sub.Lp50$Age==5,],type="b",col=3,pch=20,cex=2,lwd=2)
legend(legend=c("Age 3","Age 4",'Age 5'), col=c(1:3),x="topleft",bty="n",lty=1,lwd=2,pch=20,cex=2)

for (i in 1:nrow(sub.Lp50[sub.Lp50$Age==3,])){
  lines(x=rep(sub.Lp50$Year[sub.Lp50$Age==3][i],2),y=c(sub.Lp50$CI2.5[sub.Lp50$Age==3][i],sub.Lp50$CI97.5[sub.Lp50$Age==3][i]),col=1,lwd=2)
}

for (i in 1:nrow(sub.Lp50[sub.Lp50$Age==4,])){
  lines(x=rep(sub.Lp50$Year[sub.Lp50$Age==4][i],2),y=c(sub.Lp50$CI2.5[sub.Lp50$Age==4][i],sub.Lp50$CI97.5[sub.Lp50$Age==4][i]),col=2,lwd=2)
}
for (i in 1:nrow(sub.Lp50[sub.Lp50$Age==5,])){
  lines(x=rep(sub.Lp50$Year[sub.Lp50$Age==5][i],2),y=c(sub.Lp50$CI2.5[sub.Lp50$Age==5][i],sub.Lp50$CI97.5[sub.Lp50$Age==5][i]),col=3,lwd=2)
}

##
par(mfrow=c(1,1))
plot(Median ~ Year, data=sub.Lp50.all[sub.Lp50.all$Age==3,],type="b",col=1,pch=20,ylim=c(100,400),xlim=c(1979,2015),main="",cex=2,lwd=2,ylab="Median Lp50.all")
points(Median ~ Year, data=sub.Lp50.all[sub.Lp50.all$Age==4,],type="b",col=2,pch=20,cex=2,lwd=2)
points(Median ~ Year, data=sub.Lp50.all[sub.Lp50.all$Age==5,],type="b",col=3,pch=20,cex=2,lwd=2)
legend(legend=c("Age 3","Age 4",'Age 5'), col=c(1:3),x="topleft",bty="n",lty=1,lwd=2,pch=20,cex=2)

for (i in 1:nrow(sub.Lp50.all[sub.Lp50.all$Age==3,])){
  lines(x=rep(sub.Lp50.all$Year[sub.Lp50.all$Age==3][i],2),y=c(sub.Lp50.all$CI2.5[sub.Lp50.all$Age==3][i],sub.Lp50.all$CI97.5[sub.Lp50.all$Age==3][i]),col=1,lwd=2)
}

for (i in 1:nrow(sub.Lp50.all[sub.Lp50.all$Age==4,])){
  lines(x=rep(sub.Lp50.all$Year[sub.Lp50.all$Age==4][i],2),y=c(sub.Lp50.all$CI2.5[sub.Lp50.all$Age==4][i],sub.Lp50.all$CI97.5[sub.Lp50.all$Age==4][i]),col=2,lwd=2)
}
for (i in 1:nrow(sub.Lp50.all[sub.Lp50.all$Age==5,])){
  lines(x=rep(sub.Lp50.all$Year[sub.Lp50.all$Age==5][i],2),y=c(sub.Lp50.all$CI2.5[sub.Lp50.all$Age==5][i],sub.Lp50.all$CI97.5[sub.Lp50.all$Age==5][i]),col=3,lwd=2)
}








abline(reg=lm(Median ~ Year, data=sub.Lp50.all[sub.Lp50.all$Age==3,]),col=1)
abline(reg=lm(Median ~ Year, data=sub.Lp50.all[sub.Lp50.all$Age==4,]),col=2)
abline(reg=lm(Median ~ Year, data=sub.Lp50.all[sub.Lp50.all$Age==5,]),col=3)

summary(lm(Median ~ Year, data=sub.Lp50.all[sub.Lp50.all$Age==3,]))
summary(lm(Median ~ Year, data=sub.Lp50.all[sub.Lp50.all$Age==4,]))
summary(lm(Median ~ Year, data=sub.Lp50.all[sub.Lp50.all$Age==5,]))

save.image("YEP_PMRNrun_NoFishing_12.8.2022.Rdata")
