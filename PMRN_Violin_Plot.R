####Make new PMRN plot as age separated violin of posteriors
#loadRDS

fit_full_fishing <- readRDS("YEPFIE_covar_enviro_full_fishing_removefishingslopes.RDS") 

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
  reframe(across(all_of(ms), ~midpoints(m=.x, TL.mm=TL.mm))) %>%
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


#New age plot
dim(mats)

Lps.posteriors <- mats %>%
  filter(Age %in% c(2:5)) %>%
  group_by(Cohort, Age) %>%
  summarize(across(all_of(ms), ~midpoints(m=.x, TL.mm=TL.mm))) #%>%
  #mutate(Lp = rep(c("Lp50","Lp25",'Lp75'), n_distinct(Cohort, Age))) %>%
  #ungroup(.) %>%
Lps.posteriors

lps.long <- Lps.posteriors %>%
  pivot_longer(cols=starts_with("m"),
               names_to="iter",
               values_to="estimate")
lps.long <- lps.long %>%
  left_join(select(filtdat, Cohort, Age, CohortYear))
lps.long <- filter(lps.long, !is.na(CohortYear))
lps.long

ages <- c("2"="Age-2","3"="Age-3","4"="Age-4","5"="Age-5")

ggplot(lps.long, aes(x=CohortYear, y=estimate, group=CohortYear, fill=factor(Age))) + 
  geom_violin(trim=T, scale="width", alpha=0.6) + 
  facet_wrap(~Age, labeller=as_labeller(ages)) + 
  ylim(c(0,500)) + 
  stat_summary(fun = median,
               fun.min = function(x) quantile(x, probs = 0.025),
               fun.max = function(x) quantile(x, probs = 0.975), size = 0.5) + 
  theme_bw() + #theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                  #   legend.position="none") + 
  scale_fill_manual(values=viridis::viridis(n=4)) + 
  labs(x="Cohort", y="Length (mm)") + 
  scale_x_continuous(breaks=seq(1900,2020,5))

lps.long
max(filter(lps.long, Age==2)$CohortYear)
sum(sort(filter(lps.long, Age==2, CohortYear==1983)$estimate) < sort(filter(lps.long, Age==2, CohortYear==2014)$estimate), na.rm=T)


#Look at overlap between posteriors of earliest, pre-fishing, and post-fishing last years
#Create function to calculate posterior overlap and plot density plots by age and cohort year


post_overlap  <- function(age, cohortyear) {
  plotdat <- filter(lps.long, Age==age, CohortYear %in% cohortyear) %>%
    mutate(fCohortYear = paste0("y",CohortYear)) %>%
    group_by(fCohortYear) %>%
    mutate(step = rep(1:18000))
  
  compdat <- plotdat %>%
    pivot_wider(id_cols = c("Age","step"), names_from="fCohortYear", values_from="estimate")
  
  postplot <- (ggplot(plotdat, aes(x=estimate, color=factor(CohortYear))) +
    geom_density(stat="density", na.rm=T, lwd=1) + xlim(c(0,500)) + theme_bw() + 
    scale_color_viridis_d(guide=guide_legend(title="Cohort")) + xlab("Total length (mm)") + ylab("Density") +
    ggtitle(paste0("    Age-",age)) + theme(plot.title=element_text(margin=margin(t=30,b=-30))))
  
comptable <- tibble(expand.grid(cohortyear, cohortyear)) %>%
  filter(Var1 < Var2) %>%
  mutate(overlap = NA)

for (i in 1:nrow(comptable)) {
  comp <- tibble(x = filter(plotdat, CohortYear==comptable$Var1[i])$estimate, 
                     y = filter(plotdat, CohortYear==comptable$Var2[i])$estimate) %>%
    filter(complete.cases(.)) %>%
    mutate(x=sort(x), y=rev(sort(y)))
  
  comptable$overlap[i] <- sum(comp$x > comp$y)/length(comp$x)

}

comptable <- comptable[order(comptable$Var1),]

complist <- list(comptable, postplot)

return(complist)

}

age2overlaps <- post_overlap(age=2, cohortyear=c(1983,1996,2006,2014))
age2overlaps

age3overlaps <- post_overlap(age=3, cohortyear=c(1982,1996,2004,2013))
age3overlaps

age4overlaps <- post_overlap(age=4, cohortyear=c(1983,1998,2005,2012))
age4overlaps

age5overlaps <- post_overlap(age=5, cohortyear=c(1984,2001,2006,2011))
age5overlaps

test <- ggarrange(age2overlaps[[2]], age3overlaps[[2]], age4overlaps[[2]], age5overlaps[[2]], common.legend=F)
test
ggsave("./Figures/Supplement_OverlapPlots.png", test, dpi=500)
