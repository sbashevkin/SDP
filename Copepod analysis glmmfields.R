require(dplyr)
require(zooper)
require(lubridate)
require(hms)
require(tidyr)
require(glmmfields)


# Reduce data to monthly resolution by picking data points closest to mid-month
PF<-zoopComb%>%
  filter(Taxname=="Pseudodiaptomus forbesi")%>%
  left_join(zoopEnvComb%>%
              select(-Source),
            by="SampleID")%>%
  filter(Year>1995)%>% # Introduced in 1987, give them 8 years to get established
  mutate(Date_num=as.numeric(Date))%>%
  mutate(SalSurf_l = log(SalSurf),
         Julian_day = yday(Date))%>%
  drop_na(Date, SalSurf, Temperature, Latitude, Longitude, CPUE)%>%
  mutate(Station=if_else(Station%in%unique(zooper::stationsEMPEZ$Station), paste(Source, Station, year(Date), month(Date), day(Date)), paste(Source, Station)))%>%
  group_by_at(vars(-CPUE, -Volume, -SampleID))%>%
  summarize(CPUE=mean(CPUE), Volume=mean(Volume))%>%
  ungroup()%>%
  mutate(Count=round(CPUE*Volume),
         Year_fac=factor(Year),
         Date=if_else(is.na(Datetime), parse_date_time(paste0(year(Date), "-", month(Date), "-", day(Date), " 12:00"), "%Y-%m-%d %H:%M", tz="America/Los_Angeles"), Datetime),
         Month=month(Date, label = T),
         Week=week(Date))%>%
  mutate(Noon_diff=abs(hms(hours=12)-as_hms(Datetime)),
         mday_15_diff=abs(mday(Date)-15))%>% # Find how far each date is from the 15th of the month
  group_by(Station, Month, Year)%>%
  filter(mday_15_diff==min(mday_15_diff))%>%
  filter(Date==min(Date))%>%
  filter(Noon_diff==min(Noon_diff))%>%
  ungroup()%>%
  mutate(Date_int=(Year-1996)*12+month(Date),
         CPUE2=CPUE+1)%>%
  mutate_at(vars(Date_num, SalSurf, SalSurf_l, Temperature, Latitude, Longitude, Julian_day), list(s=~(.-mean(., na.rm=T))/sd(., na.rm=T)))

m<-glmmfields(CPUE2 ~ Month + poly(SalSurf_l_s, 2) + Source, 
              time="Date_int", lat="Latitude", lon="Longitude", 
              data=PF, iter=1e3, chains=2, cores=2, family=lognormal())
