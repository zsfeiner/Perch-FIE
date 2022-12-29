###GRAVEYARD SCRIPT FOR OLD WAY OF SUMMARIZING PMRNS

#Test against old way - works, send this to graveyard
nAges <- length(unique(w$Age[w$Age %in% c(3:5)]))
nCohorts <- length(unique(w$Cohort[w$Age %in% c(3:5)]))
nReps <- dim(m)[1]
Ages <- sort(unique(w$Age[w$Age %in% c(3:5)]))
Cohorts <- sort(unique(w$Cohort[w$Age %in% c(3:5)]))

thing <- as.data.frame(cbind(t(m),TL.mm,Age,Cohort))
head(thing)
dim(thing)

plot(thing[thing$Cohort==6,1] ~ TL.mm, data=thing[thing$Cohort==6,],col=thing$Age[thing$Cohort==6],ylim=c(-0.5,1),xlim=c(50,350))

nAges=length(2:5)
nCohorts=length(1:3)
nReps=3000
#Make array of i ages, j cohorts, and k replicates
Lp50 <- array(0, dim = c(nAges,nCohorts,nReps))
dim(Lp50)
Lp75 <- array(0, dim = c(nAges,nCohorts,nReps))
dim(Lp75)
Lp25 <- array(0, dim = c(nAges,nCohorts,nReps))
dim(Lp25)
Lp50.all<- array(0, dim = c(nAges,nCohorts,nReps))
dim(Lp50.all)


ms <- t(Lps)

plot(ms[-c(1,2),4] ~ Lp50[1,3,])

######SPLINE INTERPOLATION OPTION DO NOT USE#####
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
#######

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
sub.Lp50 <- subset(PMRNs.Lp50, NAs < 150)
sub.Lp50

test <- ggplot(data=sub.Lp50, aes(x=Cohort, y=Median, color=factor(Age))) + 
  geom_point() + geom_line() + ylim(0, 500) + theme_classic()

test
test +  geom_line(data = filter(Lps, Lp=="Lp50", NAs < 3000*0.05), aes(x=Cohort, y=median, group=factor(Age)), color="black", inherit.aes=F)

ggplot(filter(Lps, Lp=="Lp50", NAs < 3000*0.05), aes(x=Cohort, y=median, color=factor(Age))) + 
  geom_point(size=2) + geom_line(size=1.25) + ylim(0,500) +
  geom_errorbar(aes(x=Cohort, ymin=CI2.5, ymax=CI97.5), width=0.05) + 
  theme_classic()


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

plot(0,0,ylim=c(100,400),xlim=c(2.5,5.5))
for (i in 1:48) {
  points(sub.Lp50$Age[sub.Lp50$Cohort==i], sub.Lp50$Median[sub.Lp50$Cohort==i],col=i,type="b",pch=20,cex=2)
  
  #lines(x=rep(sub.Lp50$Age[sub.Lp50$Cohort==i],2),y=c(sub.Lp50$CI2.5[sub.Lp50$Cohort==i],sub.Lp50$CI97.5[sub.Lp50$Cohort==i]),col=i)
}

library(ggplot2)
ggplot(data=sub.Lp50,aes(x=Age,y=Median,group=as.factor(Year)))+
  geom_line(aes(color=as.factor(Year),size=1.10)) + geom_point()

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

save.image("InitialPMRNrun_1.7.2019.Rdata")