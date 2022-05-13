#Add in and fix data

dat <- read.table("EscanabaWAE.txt", header=TRUE)
head(dat)

#Subset angling data (only data with maturity status)
rec <- subset(dat, Gear==1)

#Subset to females
fem <- subset(rec, OtherSex=="F")
head(fem)

#Check for errors
plot(Wt.g ~ TL.mm, data=fem)
plot(TL.mm ~ Age, data=fem)

#Add in ages using year-specific length-age key
table(fem$Year)

library(measurements)
library(FSA)
library(ggplot2)
years <- unique(fem$Year)
complete.ages <- fem[0,]

for (i in 1:length(years)) {
  
  lake <- subset(fem, Year==years[i])
  
  if (sum(is.na(lake$Age)) < 1 | sum(!is.na(lake$Age)) < 1) {
    complete.ages <- rbind(complete.ages, lake)
    
  } else {
    
    ages <- lake[!is.na(lake$Age),]
    lens <- lake[is.na(lake$Age),]
    
    startcat <- round(min(ages$TL.mm,na.rm=TRUE),0)-5
    len.cats <- lencat(~TL.mm, data=ages, w=10, startcat=startcat)
    
    raw.key <- with(len.cats, table(LCat, Age))
    prop.key <- prop.table(raw.key, margin=1)
    
    lens.short <- lens[lens$TL.mm < min(len.cats$LCat,na.rm=TRUE),] 
    lens.long <- lens[lens$TL.mm >= min(len.cats$LCat,na.rm=TRUE),]
    
    lens1 <- alkIndivAge(prop.key, Age ~ TL.mm, data=lens.long)
    lens2 <- rbind(lens.short,lens1)
    
    comb <- as.data.frame(rbind(ages, lens2))
    complete.ages <- rbind(complete.ages, comb)
  }
}

#Clean out any unaged or unknown maturity fish in old data
data <- subset(complete.ages, !is.na(Age))
data <- subset(data, !is.na(Maturity) & Maturity != 3 & Year < 2005)
summary(data)

ggplot(data, aes(x=TL.mm, y=Maturity)) + 
  geom_point() + facet_wrap(~Year)

#Create cohorts and set maturity as 0 or 1
data$Cohort <- data$Year - data$Age
data$Mat <- data$Maturity - 1


#########Random effect intercept and slopes in stan with interaction and correlation matrix###########################################
table(data$Mat, data$Age)

#Visualize maturity by age
ggplot(data=data, aes(x=TL.mm, y=Maturity)) + geom_point() + facet_wrap(~Age)

#Just use fish with decent number of immatures
set <- data[data$Age %in% c(2:7),]
nrow(set)

set$Cohort <- as.numeric(factor(set$Cohort))

str(set)
#Logistic regression in stan
# create a N x k matrix of covariates
N = nrow(set)

#Explanatory variables
Age <- (set$Age)
TL.mm <- set$TL.mm
mean_TL <- mean(TL.mm)
sd_TL <- sd(TL.mm)
Cohort <- set$BinCohort
nCohorts <- length(unique(set$Cohort))
X <- cbind(Age,TL.mm)
K <- ncol(X)
Mat <- set$Mat

# # Run lm for later comparison; but go ahead and examine now if desired
# modlm = glmer(Mat ~ (1+Age*scale(TL.mm)|as.factor(Cohort))+Age*scale(TL.mm), family="binomial",control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
# summary(modlm)
# 
# modlm.2 = glmer(Mat ~ (1+Age*scale(TL.mm)|as.factor(Cohort))+Age+scale(TL.mm), family="binomial",control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
# summary(modlm.2)
# 
# modlm.3 = glmer(Mat ~ (1+Age+scale(TL.mm)|as.factor(Cohort))+Age+scale(TL.mm), family="binomial",control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
# summary(modlm.3)
# 
# AIC(modlm, modlm.2,modlm.3)
# 
# #thing <- nls(TL.mm ~ exp(phi + gamma*log(Age)), start=list('phi'=5,'gamma'=1))
# #summary(thing)
#Use linear model instead

#Create datalist
dat <- list('N'=N, 'K'=K, 'Mat'=Mat, 'Age'=Age, 'TL'=TL.mm, 'Cohort'=Cohort, 'nCohorts'=nCohorts, 'mean_TL'=mean_TL, 'sd_TL'=sd_TL)
library(rstan)
library(coda)

stanmatcode <- "
//Data block
data {
int<lower=1> N;   //Sample size
int<lower=1> K;   //Number of continuous parameters
vector[N] Age;   //Age, explanatory variable
vector[N] TL;    //Length, response variable
int<lower=1> nCohorts;   //Number of cohorts
int<lower=1,upper=nCohorts> Cohort[N];  //Cohort ID
int<lower=0,upper=1> Mat[N];   // Maturation status
real mean_TL;    //mean total length for scaling
real<lower=0> sd_TL;  //sd total length for scaling
}


transformed data {
vector[N] sc_TL;
sc_TL = (TL - mean_TL)/sd_TL;
}


//Parameter block
parameters {
//Maturation parameters
vector[K+2] beta;  //Mean coefficient vector of length K+2 (intercept + age + TL + age:TL)
matrix[4,nCohorts] z_u;  //Random effects (intercept, age effect, TL effect, interaction effect)
vector<lower=0>[4] sigma_u;  //Variance for random effects (intercept, age effect, TL effect, interaction effect)
cholesky_factor_corr[4] L_u;

//Growth parameters
vector[nCohorts] phi;  //Coefficient vector of length nCohorts
vector[nCohorts] gamma;  //Coefficient vector of length nCohorts
real phi_mu;
real<lower=0> phi_sigma;
real gamma_mu;
real<lower=0> gamma_sigma;
real<lower=0> sigma;    //Error, lower limit zero
}

transformed parameters {
matrix[4,nCohorts] u;
u = diag_pre_multiply(sigma_u, L_u) * z_u;  //Cohort random effects
}

//Model block
//Maturation model
model {
real mu[N];  //Maturation
vector[N] y;  //Growth

for (i in 1:N) {
//Maturation
mu[i] = (beta[1] + u[1,Cohort[i]]) + 
(beta[2] + u[2,Cohort[i]]) * Age[i] +
(beta[3] + u[3,Cohort[i]]) * sc_TL[i] +
(beta[4] + u[4,Cohort[i]]) * Age[i] * sc_TL[i];

//Growth
y[i] = exp(phi[Cohort[i]] + gamma[Cohort[i]] * log(Age[i]));

}

Mat ~ bernoulli_logit(mu);  //Maturation
TL ~ normal(y, sigma);  //Growth

//Priors
//Maturation
beta ~ normal(0,10);
L_u ~ lkj_corr_cholesky(2.0);
to_vector(z_u) ~ normal(0,1);

//Growth
phi ~ normal(phi_mu, phi_sigma);
phi_mu ~ normal(0, 10);
phi_sigma ~ cauchy(0,5);

gamma ~ normal(gamma_mu, gamma_sigma);
gamma_mu ~ normal(0,10);
gamma_sigma ~ cauchy(0,5);

sigma ~ cauchy(0,5);
}


generated quantities {
vector[N] s;
vector<lower=0,upper=1>[N] p;
vector<lower=0,upper=1>[N] prev_p;
vector[N] m;

for (i in 1:N) {
s[i] = ((TL[i] - (1 - pow((Age[i] - 1)/Age[i],gamma[Cohort[i]])) * TL[i]) - mean_TL) / sd_TL;

p[i] = inv_logit((beta[1] + u[1,Cohort[i]]) + 
(beta[2] + u[2,Cohort[i]]) * Age[i] +
(beta[3] + u[3,Cohort[i]]) * sc_TL[i] +
(beta[4] + u[4,Cohort[i]]) * Age[i] * sc_TL[i]);

prev_p[i] = inv_logit((beta[1] + u[1,Cohort[i]]) + 
(beta[2] + u[2,Cohort[i]]) * (Age[i] - 1) +
(beta[3] + u[3,Cohort[i]]) * s[i] +
(beta[4] + u[4,Cohort[i]]) * (Age[i] - 1) * s[i]);

m[i] = (p[i] - prev_p[i]) / (1 - prev_p[i]);
}
} 

"

inits <- function() {
  beta <- c(rnorm(1,-3,0.05),rnorm(1,1,0.05),rnorm(1,2,0.05),rnorm(1,0,0.05))
  sigma_u <- c(runif(1,1,5),runif(1,0.1,1),runif(1,0.5,2),runif(1,0.1,1))
  phi_mu <- rnorm(1,5.5,0.05)
  phi_sigma <- runif(1,0.5,1)
  gamma_mu <- rnorm(1,0.5,0.05)
  gamma_sigma <- runif(1,0.5,2)
  sigma <- runif(1,1,10)
  z_u <- matrix(rnorm(4*nCohorts,0,0.5),nrow=4,ncol=nCohorts)
  init <- list("beta"=beta, "sigma_u"=sigma_u, "phi_mu"=phi_mu, "phi_sigma"=phi_sigma, "gamma_mu"=gamma_mu, "gamma_sigma"=gamma_sigma, "sigma"=sigma, "z_u"=z_u)
  return(init)
  }



fit = stan(model_code = stanmatcode, data=dat, init=inits, 
           iter=4000, warmup=2000, thin=8, chains=4, control=list(adapt_delta=0.90,max_treedepth=15))
print(fit, pars=c('beta','sigma_u','phi_mu','gamma_mu','sigma'), digits=3, prob=c(0.025,0.5,0.975))
print(fit, pars=c('m'), digits=3, prob=c(0.025,0.5,0.975))
stan_trace(fit,pars=c('p[1]','p[2]','p[3]','p[4]','p[5]','p[6]'))
stan_trace(fit,pars=c('m[100]','m[200]','m[300]','m[400]','m[500]','m[600]'))
stan_trace(fit,pars=c('beta','gamma_mu','phi_mu'))
stan_trace(fit,pars=c('L_u'))
pairs(fit,pars=c('beta','gamma_mu','phi_mu'))

plot(set$TL.mm[set$Cohort==9],set$Mat[set$Cohort==9],col=set$Age[set$Cohort==9])
## Use the L matrices to compute the correlation matrices
# L matrices
#L_u <- extract(fit, pars = "L_u")$L_u

# correlation parameters
#cor_u <- apply(L_u, 1, function(x) tcrossprod(x)[4, 2])
#print(signif(quantile(cor_u, probs = c(0.025, 0.5, 0.975)), 2))
#print(mean(cor_u))


#Mess around with plotting maturation ogives
m <- extract(fit,pars='m')$m
m

w <- cbind(m[1,],TL.mm,Age, Cohort)
table(w[w[,1]<0,4])

plot(set$TL.mm[set$Cohort==2],set$Mat[set$Cohort==2],col=set$Age[set$Cohort==2])
plot(set$TL.mm[set$Cohort==5],set$Mat[set$Cohort==5],col=set$Age[set$Cohort==5])
plot(set$TL.mm[set$Cohort==7],set$Mat[set$Cohort==7],col=set$Age[set$Cohort==7])
plot(set$TL.mm[set$Cohort==9],set$Mat[set$Cohort==9],col=set$Age[set$Cohort==9])
plot(set$TL.mm[set$Cohort==12],set$Mat[set$Cohort==12],col=set$Age[set$Cohort==12])
plot(set$TL.mm[set$Cohort==14],set$Mat[set$Cohort==14],col=set$Age[set$Cohort==14])
plot(set$TL.mm[set$Cohort==8],set$Mat[set$Cohort==8],col=set$Age[set$Cohort==8])
plot(set$TL.mm[set$Cohort==20],set$Mat[set$Cohort==20],col=set$Age[set$Cohort==20])
plot(set$TL.mm[set$Cohort==22],set$Mat[set$Cohort==22],col=set$Age[set$Cohort==22])


plot(w[w[,4]==5 & w[,3]==4,1] ~ w[w[,4]==5 & w[,3]==4,2],col=w[w[,4]==5 & w[,3]==4,3])

points(spline(x=w[w[,4]==5 & w[,3]==4,2], y=w[w[,4]==5 & w[,3]==4,1]))

spline(y=w[w[,4]==5 & w[,3]==4,2], x=w[w[,4]==5 & w[,3]==4,1], xout=0.5)
palette(c("black","red","green","blue","gray","magenta","orange","purple"))
plot(w[w[,4]==1,1] ~ w[w[,4]==1,2],col=w[w[,4]==1,3])
plot(w[w[,4]==6,1] ~ w[w[,4]==6,2],col=w[w[,4]==6,3])
plot(w[w[,4]==8,1] ~ w[w[,4]==8,2],col=w[w[,4]==8,3])
plot(w[w[,4]==9,1] ~ w[w[,4]==9,2],col=w[w[,4]==9,3])
plot(w[w[,4]==14,1] ~ w[w[,4]==14,2],col=w[w[,4]==14,3])
plot(w[w[,4]==18,1] ~ w[w[,4]==18,2],col=w[w[,4]==18,3])
plot(w[w[,4]==23,1] ~ w[w[,4]==23,2],col=w[w[,4]==23,3])
plot(w[w[,4]==25,1] ~ w[w[,4]==25,2],col=w[w[,4]==25,3])

plot(w[w[,3]==3,1] ~ w[w[,3]==3,2],col=w[w[,3]==3,4],ylim=c(0,1))
plot(w[w[,3]==4,1] ~ w[w[,3]==4,2],col=w[w[,3]==4,4])
plot(w[w[,3]==5,1] ~ w[w[,3]==5,2],col=w[w[,3]==5,4])
plot(w[w[,3]==6,1] ~ w[w[,3]==6,2],col=w[w[,3]==6,4])
plot(w[w[,3]==7,1] ~ w[w[,3]==7,2],col=w[w[,3]==7,4])
table(data$BinCohort, data$Age)

#Use spline interpolation to determine PMRN midpoints for ages 3-5 (older ages declining in m, already mature)
head(w)
w <- as.data.frame(w)
head(w)
names(w) <- c('m','TL.mm','Age','BinCohort')
head(w)

nAges <- length(unique(w$Age[w$Age %in% c(3:5)]))
nCohorts <- length(unique(w$BinCohort[w$Age %in% c(3:5)]))
nReps <- dim(m)[1]
Ages <- sort(unique(w$Age[w$Age %in% c(3:5)]))
Cohorts <- sort(unique(w$BinCohort[w$Age %in% c(3:5)]))

thing <- as.data.frame(cbind(t(m),TL.mm,Age,Cohort))
head(thing)
dim(thing)

plot(thing[thing$Cohort==1,1] ~ TL.mm, data=thing[thing$Cohort==1,],col=thing$Age[thing$Cohort==1],ylim=c(-2,1),xlim=c(300,800))
for (i in 2:1000) {
  points(thing[thing$Cohort==1,i] ~ TL.mm, data=thing[thing$Cohort==1,],col=thing$Age[thing$Cohort==1])
}


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
#for (k in 1:nReps) {
#      
#  d <- as.data.frame(cbind(m[k,],TL.mm,Age,Cohort))
#  names(d)[1] <- "m"
#  
#  for (i in 1:nAges) {
#    for (j in 1:nCohorts) {
#      sub.d <- subset(d, Age == Ages[i] & Cohort == Cohorts[j] & m > 0)
#      
#      if (nrow(sub.d)<2 | length(unique(sub.d$TL.mm))==1) {
#        Lp50[i,j,k] <- NA
#        Lp75[i,j,k] <- NA
#        Lp25[i,j,k] <- NA
#        Lp50.all[i,j,k] <- NA
#      } else {
#
    #  reg <- glm(m ~ TL.mm, family=binomial, data=sub.d)
    #  if (coef(reg)[2] < 0 | min(sub.d$m, na.rm=TRUE) > 0.5 | max(sub.d$m, na.rm=TRUE) < 0.5) { 
   #     Lp50[i,j,k] <- NA } else {
  #        midpoint <- spline(sub.d$m,sub.d$TL.mm,xout=0.5)
   #       Lp50[i,j,k] <-  midpoint$y
  #      }
  #    if (coef(reg)[2] < 0 | min(sub.d$m, na.rm=TRUE) > 0.75 | max(sub.d$m, na.rm=TRUE) < 0.75) { 
  #      Lp75[i,j,k] <- NA } else {
 #         midpoint <- spline(sub.d$m,sub.d$TL.mm,xout=0.75)
#          Lp75[i,j,k] <-  midpoint$y
#        }
#      if (coef(reg)[2] < 0 | min(sub.d$m, na.rm=TRUE) > 0.25 | max(sub.d$m, na.rm=TRUE) < 0.25) { 
#        Lp25[i,j,k] <- NA } else {
#          midpoint <- spline(sub.d$m,sub.d$TL.mm,xout=0.25)
#          Lp25[i,j,k] <-  midpoint$y
#        }
#      if (coef(reg)[2] < 0) { 
#        Lp50.all[i,j,k] <- NA } else {
#          midpoint <- spline(sub.d$m,sub.d$TL.mm,xout=0.5)
#          Lp50.all[i,j,k] <-  midpoint$y
#        }
#    }
#  }
#}
#}
#Use logistic regression
for (k in 1:nReps) {
  
  d <- as.data.frame(cbind(m[k,],TL.mm,Age,Cohort))
  names(d)[1] <- "m"
  
  for (i in 1:nAges) {
    for (j in 1:nCohorts) {
      sub.d <- subset(d, Age == Ages[i] & Cohort == Cohorts[j])
      sub.d$m[sub.d$m<0] <- 0
      
      if (nrow(sub.d)<2 | length(unique(sub.d$TL.mm))==1) {
        Lp50[i,j,k] <- NA
        Lp75[i,j,k] <- NA
        Lp25[i,j,k] <- NA
        Lp50.all[i,j,k] <- NA
      } else {
        
        reg <- suppressWarnings(glm(m ~ TL.mm, family=binomial, data=sub.d))
        if (coef(reg)[2] < 0) { 
          Lp50[i,j,k] <- NA
          Lp75[i,j,k] <- NA 
          Lp25[i,j,k] <- NA 
          
          } else {
            Lp50[i,j,k] <- -coef(reg)[1]/coef(reg)[2]
            Lp75[i,j,k] <- (log(1/.75-1) + coef(reg)[1])/-coef(reg)[2]
            Lp25[i,j,k]  <- (log(1/.25-1) + coef(reg)[1])/-coef(reg)[2]
            }
      }
    }
  }
}

colSums(is.na(Lp50[1,,]))
colSums(is.na(Lp50[2,,]))
colSums(is.na(Lp50[3,,]))
dim(Lp50)
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
  temp.Lp50 <- matrix(0,nrow=3,ncol=7)
  colnames(temp.Lp50) <- colnames(PMRNs.Lp50)
  z.Lp50 <- Lp50[,i,]
  temp.Lp50[,1] <- i
  temp.Lp50[,2] <- c(3,4,5)
  temp.Lp50[,3] <- rowSums(is.na(z.Lp50))
  temp.Lp50[,4] <- rowMeans(z.Lp50,na.rm=TRUE)
  temp.Lp50[,5:7] <- t(apply(z.Lp50,1,FUN=credible.intervals))
  PMRNs.Lp50 <- rbind(PMRNs.Lp50,temp.Lp50)
  
  #Summarize Lp25s
  temp.Lp25 <- matrix(0,nrow=3,ncol=7)
  colnames(temp.Lp25) <- colnames(PMRNs.Lp25)
  z.Lp25 <- Lp25[,i,]
  temp.Lp25[,1] <- i
  temp.Lp25[,2] <- c(3,4,5)
  temp.Lp25[,3] <- rowSums(is.na(z.Lp25))
  temp.Lp25[,4] <- rowMeans(z.Lp25,na.rm=TRUE)
  temp.Lp25[,5:7] <- t(apply(z.Lp25,1,FUN=credible.intervals))
  PMRNs.Lp25 <- rbind(PMRNs.Lp25,temp.Lp25)
  
  #Summarize Lp75s
  temp.Lp75 <- matrix(0,nrow=3,ncol=7)
  colnames(temp.Lp75) <- colnames(PMRNs.Lp75)
  z.Lp75 <- Lp75[,i,]
  temp.Lp75[,1] <- i
  temp.Lp75[,2] <- c(3,4,5)
  temp.Lp75[,3] <- rowSums(is.na(z.Lp75))
  temp.Lp75[,4] <- rowMeans(z.Lp75,na.rm=TRUE)
  temp.Lp75[,5:7] <- t(apply(z.Lp75,1,FUN=credible.intervals))
  PMRNs.Lp75 <- rbind(PMRNs.Lp75,temp.Lp75)
  
  #Summarize Lp50.alls
  temp.Lp50.all <- matrix(0,nrow=3,ncol=7)
  colnames(temp.Lp50.all) <- colnames(PMRNs.Lp50.all)
  z.Lp50.all <- Lp50.all[,i,]
  temp.Lp50.all[,1] <- i
  temp.Lp50.all[,2] <- c(3,4,5)
  temp.Lp50.all[,3] <- rowSums(is.na(z.Lp50.all))
  temp.Lp50.all[,4] <- rowMeans(z.Lp50.all,na.rm=TRUE)
  temp.Lp50.all[,5:7] <- t(apply(z.Lp50.all,1,FUN=credible.intervals))
  PMRNs.Lp50.all <- rbind(PMRNs.Lp50.all,temp.Lp50.all)
  
}

PMRNs.Lp50
PMRNs.Lp75
PMRNs.Lp25
PMRNs.Lp50.all

PMRNs.Lp50 <- as.data.frame(PMRNs.Lp50)
sub.Lp50 <- subset(PMRNs.Lp50, NAs < 500)
sub.Lp50

PMRNs.Lp75 <- as.data.frame(PMRNs.Lp75)
sub.Lp75 <- subset(PMRNs.Lp75, NAs < 500)
sub.Lp75

PMRNs.Lp25 <- as.data.frame(PMRNs.Lp25)
sub.Lp25 <- subset(PMRNs.Lp25, NAs < 500)
sub.Lp25

PMRNs.Lp50.all <- as.data.frame(PMRNs.Lp50.all)
sub.Lp50.all <- subset(PMRNs.Lp50.all, NAs < 500 & Median<1000 & CI2.5 > 0 & CI97.5 < 1000)
sub.Lp50.all


#Match up years
yearnames <- cbind(Cohort, set$Year-set$ObsAge)
colnames(yearnames) <-c("Cohort","Year")
yearnames <- as.data.frame(yearnames)
sub.Lp50$Year <- yearnames$Year[match(sub.Lp50$Cohort, yearnames$Cohort)]
sub.Lp75$Year <- yearnames$Year[match(sub.Lp75$Cohort, yearnames$Cohort)]
sub.Lp25$Year <- yearnames$Year[match(sub.Lp25$Cohort, yearnames$Cohort)]
sub.Lp50.all$Year <- yearnames$Year[match(sub.Lp50.all$Cohort, yearnames$Cohort)]

sub.Lp50

plot(0,0,ylim=c(100,800),xlim=c(2.5,5.5))
for (i in 1:48) {
  points(sub.Lp50$Age[sub.Lp50$Cohort==i], sub.Lp50$Median[sub.Lp50$Cohort==i],col=i,type="b",pch=20,cex=2)
  #lines(x=rep(sub.Lp50$Age[sub.Lp50$Cohort==i],2),y=c(sub.Lp50$CI2.5[sub.Lp50$Cohort==i],sub.Lp50$CI97.5[sub.Lp50$Cohort==i]),col=i)
}

par(mfrow=c(1,3))
plot(Median ~ Year, data=sub.Lp50[sub.Lp50$Age==3,],type="b",col=1,pch=20,ylim=c(250,500),xlim=c(1950,2000),main="Age 3",cex=2,lwd=2,ylab="Median Lp50")
lines(CI2.5 ~ Year, data=sub.Lp50[sub.Lp50$Age==3,],lty=2,col=1)
lines(CI97.5 ~ Year, data=sub.Lp50[sub.Lp50$Age==3,],lty=2,col=1)
abline(reg=lm(Median ~ Year, data=sub.Lp50[sub.Lp50$Age==3,]),col=1)

plot(Median ~ Year, data=sub.Lp50[sub.Lp50$Age==4,],type="b",col=2,pch=20,ylim=c(250,500),xlim=c(1950,2000),main="Age 4",cex=2,lwd=2,ylab="Median Lp50")
lines(CI2.5 ~ Year, data=sub.Lp50[sub.Lp50$Age==4,],lty=2,col=2)
lines(CI97.5 ~ Year, data=sub.Lp50[sub.Lp50$Age==4,],lty=2,col=2)
abline(reg=lm(Median ~ Year, data=sub.Lp50[sub.Lp50$Age==4,]),col=2)

plot(Median ~ Year, data=sub.Lp50[sub.Lp50$Age==5,],type="b",col=3,pch=20,ylim=c(250,500),xlim=c(1950,2000), main="Age 5",cex=2,lwd=2,ylab="Median Lp50")
lines(CI2.5 ~ Year, data=sub.Lp50[sub.Lp50$Age==5,],lty=2,col=3)
lines(CI97.5 ~ Year, data=sub.Lp50[sub.Lp50$Age==5,],lty=2,col=3)
abline(reg=lm(Median ~ Year, data=sub.Lp50[sub.Lp50$Age==5,]),col=3)

par(mfrow=c(1,1))
plot(Median ~ Year, data=sub.Lp50[sub.Lp50$Age==3,],type="b",col=1,pch=20,ylim=c(250,500),xlim=c(1950,2000),main="Age 3",cex=2,lwd=2,ylab="Median Lp50")
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
plot(Median ~ Year, data=sub.Lp75[sub.Lp75$Age==3,],type="b",col=1,pch=20,ylim=c(250,500),xlim=c(1950,2000),main="Age 3",cex=2,lwd=2,ylab="Median Lp75")
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
plot(Median ~ Year, data=sub.Lp25[sub.Lp25$Age==3,],type="b",col=1,pch=20,ylim=c(250,500),xlim=c(1950,2000),main="Age 3",cex=2,lwd=2,ylab="Median Lp25")
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
plot(Median ~ Year, data=sub.Lp50[sub.Lp50$Age==4,],type="b",col=1,pch=20,ylim=c(250,500),xlim=c(1950,2000),main="Age 3",cex=2,lwd=2,ylab="Median Lp25")
points(Median ~ Year, data=sub.Lp25[sub.Lp25$Age==4,],type="b",col=2,pch=20,cex=2,lwd=2)
points(Median ~ Year, data=sub.Lp75[sub.Lp75$Age==4,],type="b",col=3,pch=20,cex=2,lwd=2)

##
par(mfrow=c(1,1))
plot(Median ~ Year, data=sub.Lp50[sub.Lp50$Age==5,],type="b",col=1,pch=20,ylim=c(250,500),xlim=c(1950,2000),main="Age 3",cex=2,lwd=2,ylab="Median Lp25")
points(Median ~ Year, data=sub.Lp25[sub.Lp25$Age==5,],type="b",col=2,pch=20,cex=2,lwd=2)
points(Median ~ Year, data=sub.Lp75[sub.Lp75$Age==5,],type="b",col=3,pch=20,cex=2,lwd=2)


##
par(mfrow=c(1,1))
plot(Median ~ Year, data=sub.Lp50[sub.Lp50$Age==3,],type="b",col=1,pch=20,ylim=c(100,700),xlim=c(1950,2000),main="",cex=2,lwd=2,ylab="Median Lp50")
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


abline(reg=lm(Median ~ Year, data=sub.Lp50.all[sub.Lp50.all$Age==3,]),col=1)
abline(reg=lm(Median ~ Year, data=sub.Lp50.all[sub.Lp50.all$Age==4,]),col=2)
abline(reg=lm(Median ~ Year, data=sub.Lp50.all[sub.Lp50.all$Age==5,]),col=3)

summary(lm(Median ~ Year, data=sub.Lp50.all[sub.Lp50.all$Age==3,]))
summary(lm(Median ~ Year, data=sub.Lp50.all[sub.Lp50.all$Age==4,]))
summary(lm(Median ~ Year, data=sub.Lp50.all[sub.Lp50.all$Age==5,]))

save.image("InitialPMRNrun_12.12.2018.Rdata")





