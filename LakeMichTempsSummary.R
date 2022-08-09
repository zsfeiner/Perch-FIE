setwd(paste0(getwd(),"/45007hTemps_1981.2021"))
library(tidyverse)
library(readr)
library(data.table)

filelist <- list.files(getwd())

filelist


#This gets kind of dumb because different years are in slightly different formats so need to splice them together
alltempsto98 <- filelist[1:18] %>%
  set_names(.) %>%
  map_df(fread, .id="FileName", colClasses="numeric") %>%
  tibble(.) %>%
  rename(YYYY=YY) %>%
  mutate(YYYY=YYYY+1900)
alltempsto98

alltemps99 <- fread("~/LakeMichTemp/45007hTemps_1981.2021/45007h1999.txt", colClasses="numeric") %>%
  tibble(.) %>%
  mutate(FileName="45007h1999.txt", .before=YYYY)
alltemps99

alltemps00.04 <- filelist[20:24] %>%
  set_names(.) %>%
  map_df(fread, .id="FileName", colClasses="numeric", drop="TIDE", fill=T) %>%
  tibble(.) 
alltemps00.04

alltemps05.08 <- filelist[25:28] %>%
  set_names(.) %>%
  map_df(fread, .id="FileName", colClasses="numeric", drop=c("TIDE", "mm")) %>%
  tibble(.) %>%
  select(-"WDIR","PRES")
alltemps05.08

alltemps09.21 <- filelist[29:41] %>%
  set_names(.) %>%
  map_df(fread, .id="FileName", colClasses="numeric", skip=2, 
         col.names=c("YYYY", "MM", "DD", "hh", "mm", "WDIR", "WSPD", "GST",  "WVHT",   
                     "DPD",   "APD", "MWD",   "PRES",  "ATMP",  "WTMP",  "DEWP",  "VIS",  "TIDE")) %>%
  tibble(.) %>%
  select(-"mm", -"TIDE")
alltemps09.21

alltemps <- bind_rows(alltempsto98,alltemps99, alltemps00.04,alltemps05.08,alltemps09.21)
alltemps
alltemps <- alltemps %>%
  mutate(DATE = as.Date(paste0(MM,"/",DD,"/",YYYY), format="%m/%d/%Y"))
alltemps

temps <- alltemps %>%
  filter(WTMP<40, !is.na(WTMP)) %>%
  group_by(DATE, YYYY, MM, DD) %>%
  summarize(Temp = mean(WTMP))
plot(Temp ~ DATE, temps)

#Have temperature measurements from July, Aug, Sept, Oct except for July 2015
plot(MM ~ YYYY, temps)



###Read data from NOAA GLERL FTP site to get St Joes Water Plant Data

files <- paste0("TM096D", seq(60,92),".EDT")
url = "https://www.glerl.noaa.gov/ftp/publications/tech_reports/glerl-096/"

library(threadr)
download_ftp_file(paste0(url, files), file_local=paste0("~/External Projects/Lake Michigan YEP FIE/Perch-FIE/StJoeTemps/",files), curl=F)

StJoe_1960.1992 <- paste0("~/External Projects/Lake Michigan YEP FIE/Perch-FIE/StJoeTemps/",files) %>%
  set_names(.) %>%
  map_df(fread, .id="FileName", colClasses="numeric", col.names=c("Month","Day","Year","Temp"), fill=T) %>%
  tibble(.) %>%
  filter(!is.na(Month)) %>%
  mutate(Date = as.Date(paste(Month,Day,Year, sep="/"), format="%m/%d/%Y"))
StJoe_1960.1992

plot(Temp ~ Date, data=StJoe_1960.1992)
lines(Temp ~ DATE, data=temps, col="red")
colSums(is.na(StJoe_1960.1992))


##Read excel files from St Joe Plant Manager
library(readxl)
#Loop to read in data, clean, and organize it
Years <- 2005:2018
fileloc <- "~/External Projects/Lake Michigan YEP FIE/Perch-FIE/StJoeTemps/"
StJoes_2005.2018 <-tibble()
for (i in 1:length(Years)) {
  set <- read_excel(paste0("~/External Projects/Lake Michigan YEP FIE/Perch-FIE/StJoeTemps/", Years[i], ".xlsx"), sheet=1, skip=1)
  names(set) <- make.names(names(set))
  set <- set %>%
    select(seq(1:ncol(.))[seq(1:ncol(.)) %% 3 != 0]) %>%
    pivot_longer(cols=everything(),
                 names_to=".value",
                 names_pattern = '([A-Za-z]+)\\d?') %>%
    mutate(Year=Years[i], Month = match(substr(Date, 1, 3), month.abb), Day = parse_number(substr(Date, 5, str_length(Date)))) %>%
    mutate(NewDate = as.Date(paste(Month, Day, Year, sep="/"), format="%m/%d/%Y"))
  
  set <- set[order(set$NewDate),]
  StJoes_2005.2018 <- bind_rows(StJoes_2005.2018, set)
}

StJoes_2005.2018$TempC <- measurements::conv_unit(StJoes_2005.2018$Temperature, from="F", to="C")

plot(TempC ~ NewDate, data=StJoes_2005.2018)

#Combine St Joes Data

StJoe_1960.1992
StJoes_2005.2018

SJ_05.18 <- StJoes_2005.2018 %>%
  select(NewDate, Month, Day, Year, TempC) %>%
  rename("Date"="NewDate","Temp"="TempC")
SJ_05.18

SJ_60.92 <- StJoe_1960.1992 %>%
  select(Date, Month, Day, Year, Temp)

StJoeTemp <- bind_rows(SJ_60.92, SJ_05.18) %>%
  mutate(Temp = ifelse(Temp<0, 0, Temp), GDD5 = ifelse(Temp>5, Temp-5, 0))

GDD5 <- aggregate(GDD5 ~ Year, data=StJoeTemp, FUN="sum")
sumtemp <- aggregate(Temp ~ Year, data=filter(StJoeTemp, Month %in% c(6,7,8,9)), FUN="mean")

plot(GDD5 ~ Year, data=GDD5)
plot(Temp ~ Year, data=sumtemp)

plot(Temp ~ Date, data=StJoeTemp, type="l")


##Graveyard - compare nearshore buoy to offshore buoy
nearfiles <- list.files("C:/Users/feinezs/Documents/LakeMichTemp/45026hTemps_2011.2021")
nearfiles

setwd("~/LakeMichTemp/45026hTemps_2011.2021/")
neartemps <- nearfiles %>%
  set_names(.) %>%
  map_df(fread, colClasses="numeric", skip=2,
         col.names=c("YYYY", "MM", "DD", "hh", "mm", "WDIR", "WSPD", "GST",  "WVHT",   
                     "DPD",   "APD", "MWD",   "PRES",  "ATMP",  "WTMP",  "DEWP",  "VIS",  "TIDE")) %>%
  tibble(.) %>%
  select(-"mm", -"TIDE")
neartemps


neartemps <- neartemps %>%
  filter(WTMP<40, !is.na(WTMP)) %>%
  mutate(DATE=as.Date(paste0(MM,"/",DD,"/",YYYY), format="%m/%d/%Y")) %>%
  group_by(DATE, YYYY, MM, DD) %>%
  summarize(NearTemp=mean(WTMP))
neartemps
plot(NearTemp ~ DATE, data=neartemps)



bothtemps <- left_join(neartemps, temps)
bothtemps

ggplot(bothtemps, aes(x=DATE, y=NearTemp)) + 
  geom_point(color="blue") + 
  geom_point(aes(x=DATE, y=Temp), color="red") + 
  theme_classic()
