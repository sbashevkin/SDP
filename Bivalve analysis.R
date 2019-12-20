require(leaflet)
require(leaflet.minicharts)
require(readxl)
require(lubridate)
require(tidyverse)
require(brms)

Fieldfiles <- list.files(path = "Data/Water quality", full.names = T, pattern="Field")

Labfiles <- list.files(path = "Data/Water quality", full.names = T, pattern="Lab")

WQ<-sapply(Fieldfiles, function(x) read_excel(x, guess_max = 5e4))%>%
  bind_rows()%>%
  select(Date=SampleDate, Station=StationCode, Parameter=AnalyteName, Value=Result, Notes=TextResult)%>%
  filter(Parameter%in%c("Temperature", "Secchi Depth", "Conductance (EC)", "Oxygen", "Depth"))%>%
  group_by(Date, Station, Parameter, Notes)%>%
  summarise(Value=mean(Value, na.rm=T))%>%
  ungroup()%>%
  bind_rows(sapply(Labfiles, function(x) read_excel(x, guess_max = 5e4))%>%
              bind_rows()%>%
              select(Station=StationCode, Date=SampleDate, Parameter=ConstituentName, Value=Result, Notes=LabAnalysisRemarks)%>%
              filter(Parameter=="Chlorophyll a")%>%
              group_by(Date, Station, Parameter, Notes)%>%
              summarise(Value=mean(Value, na.rm=T))%>%
              ungroup())%>%
  spread(key=Parameter, value=Value)%>%
  rename(Chlorophyll=`Chlorophyll a`, Secchi_depth=`Secchi Depth`, Conductivity=`Conductance (EC)`)%>%
  bind_rows(read_excel("Data/Water quality/EMP WQ Combined_2000-2018.xlsx", na=c("N/A", "<R.L.", "Too dark"), col_types = c(rep("text", 3), "date", rep("text", 37)))%>%
              select(Station=`Station Name`, Date, Chlorophyll=starts_with("Chlorophyll"), Latitude=`North Latitude Degrees (d.dd)`, Longitude=`West Longitude Degrees (d.dd)`, Secchi_depth=`Secchi Depth Centimeters`, Temperature=starts_with("Water Temperature"), Conductivity=starts_with("Specific Conductance"),
                     Oxygen=starts_with("Dissolved Oxygen"), Depth=starts_with("Water Depth"), Oxygen_Bottom=starts_with("(Bottom) Dissolved"), Conductivity_Bottom=starts_with("(Bottom) Specific"), Temperature_Bottom=starts_with("(Bottom) Water Temperature"), Turbidity_Bottom=starts_with("(Bottom) Turbidity"), Fluorescence_Bottom=starts_with("(Bottom) Fluorescence"))%>%
              mutate(Chlorophyll=parse_double(ifelse(Chlorophyll%in%c("<0.05", "<0.5"), 0, Chlorophyll)),
                     Turbidity_Bottom=parse_double(ifelse(Turbidity_Bottom=="<0.10000000000000001", 0, Turbidity_Bottom)))%>%
              mutate_at(c("Secchi_depth", "Temperature", "Conductivity", "Oxygen", "Depth", "Temperature_Bottom", "Conductivity_Bottom", "Oxygen_Bottom", "Fluorescence_Bottom"), parse_double))%>%
  mutate(MonthYear=floor_date(Date, unit = "month"),
         Year=year(Date),
         Salinity=((0.36966/(((Conductivity*0.001)^(-1.07))-0.00074))*1.28156),
         Salinity_Bottom=((0.36966/(((Conductivity_Bottom*0.001)^(-1.07))-0.00074))*1.28156))%>%
  select(-Conductivity, -Conductivity_Bottom)%>%
  group_by(MonthYear, Year, Station)%>%
  summarise_at(vars(c("Chlorophyll", "Secchi_depth", "Temperature", "Salinity", "Oxygen", "Depth", "Temperature_Bottom", "Salinity_Bottom", "Oxygen_Bottom", "Fluorescence_Bottom", "Turbidity_Bottom")), ~mean(., na.rm=T))%>%
  ungroup()


#Stations<-read_excel("Bivalves app/1975-18 CPUE bivalves only, 2019Sept9.xlsx",
#                     sheet = "75-17 station locations", skip=1)%>%
#  select(BivStation=Site_Code, Latitude, Longitude)

Stations<-read_csv("Data/Water quality/Master station key.csv",
                   col_types = "ccddc")%>%
  select(-StationID)%>%
  filter(Source=="EMP")%>%
  drop_na()

Biv<-read_excel("Bivalves app/1975-18 CPUE bivalves only, 2019Sept9.xlsx",
                sheet = "75-18 CPUE per m2", skip=1)%>%
  select(Date, BivStation=StationCode, `Potamocorbula amurensis`, `Corbicula fluminea`)%>%
  gather(key="Taxa", value="CPUE", -BivStation, -Date)%>%
  mutate(Year=year(Date),
         MonthYear=floor_date(Date, unit = "month"))%>%
  separate(BivStation, into=c("Station", "Position"), sep="-", remove=F)%>%
  group_by(Year, MonthYear, Taxa, Station)%>%
  summarise(CPUE=mean(CPUE, na.rm=T))%>%
  ungroup()%>%
  left_join(WQ, by=c("Year", "MonthYear", "Station"))%>%
  #left_join(Stations, by="BivStation")%>%
  left_join(Stations, by="Station")%>%
  filter(!is.na(Salinity) & !is.na(Oxygen))%>%
  mutate(Date_num=as.numeric(MonthYear))%>%
  mutate_at(vars(c("Chlorophyll", "Secchi_depth", "Temperature", "Salinity", "Salinity_Bottom", "Oxygen", "Depth", "Temperature_Bottom", "Oxygen_Bottom", "Fluorescence_Bottom", "Turbidity_Bottom", "Date_num")), list(s=~(.-mean(., na.rm=T))/sd(., na.rm=T)))%>%
  mutate(Month=month(MonthYear, label=T))

PA<-Biv%>%
  filter(Taxa=="Potamocorbula amurensis" & Year>1990)

model<-brm(CPUE~Salinity_s, data=PA, family=hurdle_gamma(),
           prior=prior(normal(0,10), class="Intercept")+
             prior(normal(0,5), class="b"),
           chains=1, cores=1,
           iter = 1e4, warmup = 2.5e3)
model<-add_criterion(model, c("waic", "loo"))

model2<-brm(CPUE~s(Salinity_s), data=PA, family=hurdle_gamma(),
           prior=prior(normal(0,10), class="Intercept")+
             prior(normal(0,5), class="b"),
           chains=1, cores=1, control=list(adapt_delta=0.99),
           iter = 1e4, warmup = 2.5e3)
model2<-add_criterion(model2, c("waic", "loo"), reloo=T)

modelln<-brm(CPUE~Salinity_s, data=PA, family=hurdle_lognormal(),
           prior=prior(normal(0,10), class="Intercept")+
             prior(normal(0,5), class="b"),
           chains=1, cores=1,
           iter = 1e4, warmup = 2.5e3)
modelln<-add_criterion(modelln, c("waic", "loo"))

modelln2<-brm(CPUE~s(Salinity_s), data=PA, family=hurdle_lognormal(),
            prior=prior(normal(0,10), class="Intercept")+
              prior(normal(0,5), class="b"),
            chains=1, cores=1, control=list(adapt_delta=0.99),
            iter = 1e4, warmup = 2.5e3)
modelln2<-add_criterion(modelln2, c("waic", "loo"), reloo=T)

#hurdle_gamma is better than lognormal by loo
#elpd_diff se_diff
#model2      0.0       0.0 
#modelln2  -42.5      24.8 
#model     -94.2      24.8 
#modelln  -227.2      28.7 

model3<-brm(CPUE~poly(Salinity_s,3), data=PA, family=hurdle_gamma(),
            prior=prior(normal(0,10), class="Intercept")+
              prior(normal(0,5), class="b"),
            chains=1, cores=1, control=list(adapt_delta=0.99),
            iter = 1e4, warmup = 2.5e3)
model3<-add_criterion(model3, c("waic", "loo"), reloo=T)


model4<-brm(CPUE~poly(Salinity_s,2), data=PA, family=hurdle_gamma(),
            prior=prior(normal(0,10), class="Intercept")+
              prior(normal(0,5), class="b"),
            chains=1, cores=1,
            iter = 1e4, warmup = 2.5e3)
model4<-add_criterion(model4, c("waic", "loo"), reloo=T)

#Model with smoother is slightly better than poly(2) model, which is slightly better than poly(3) model. Smoother probably not worth the computational time though.
#elpd_diff se_diff
#model2   0.0       0.0  
#model4  -5.6       6.2  
#model3 -10.9       7.0  
#model  -94.2      24.8

model5<-brm(CPUE~poly(Salinity_s,2)+Month, data=PA, family=hurdle_gamma(),
            prior=prior(normal(0,10), class="Intercept")+
              prior(normal(0,5), class="b"),
            chains=1, cores=1,
            iter = 1e4, warmup = 2.5e3)
model5<-add_criterion(model5, c("waic", "loo"), reloo=T)
#Month has no effect and is a worse model

model6<-brm(CPUE~poly(Salinity_s,2)+poly(Oxygen_s,2), data=PA, family=hurdle_gamma(),
            prior=prior(normal(0,10), class="Intercept")+
              prior(normal(0,5), class="b"),
            chains=1, cores=1,
            iter = 1e4, warmup = 2.5e3)
model6<-add_criterion(model6, c("waic", "loo"), reloo=T)

model7<-brm(CPUE~poly(Salinity_s,2)+Oxygen_s, data=PA, family=hurdle_gamma(),
            prior=prior(normal(0,10), class="Intercept")+
              prior(normal(0,5), class="b"),
            chains=1, cores=1,
            iter = 1e4, warmup = 2.5e3)
model7<-add_criterion(model7, c("waic", "loo"), reloo=T)

model8<-brm(CPUE~s(Salinity_s)+s(Oxygen_s), data=PA, family=hurdle_gamma(),
            prior=prior(normal(0,10), class="Intercept")+
              prior(normal(0,5), class="b"),
            chains=1, cores=1,
            iter = 1e4, warmup = 2.5e3)
model8<-add_criterion(model8, c("waic", "loo"), reloo=T)

model9<-brm(CPUE~poly(Salinity_s,2) + (1|Station), data=PA, family=hurdle_gamma(),
            prior=prior(normal(0,10), class="Intercept")+
              prior(normal(0,5), class="b"),
            chains=1, cores=1, control=list(adapt_delta=0.99), 
            iter = 1e4, warmup = 2.5e3)
model9<-add_criterion(model9, c("waic", "loo"), reloo=T)

model10<-brm(CPUE~poly(Salinity_s,2) + t2(Latitude, Longitude, Date_num_s), data=PA, family=hurdle_gamma(),
            prior=prior(normal(0,10), class="Intercept")+
              prior(normal(0,5), class="b"),
            chains=1, cores=1, control=list(adapt_delta=0.99), 
            iter = 1e4, warmup = 2.5e3)
model10<-add_criterion(model10, c("waic", "loo"), reloo=T)
