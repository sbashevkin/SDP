require(leaflet)
require(leaflet.minicharts)
require(readxl)
require(lubridate)
require(tidyverse)

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
                     Oxygen=starts_with("Dissolved Oxygen"), Depth=starts_with("Water Depth"))%>%
              mutate(Chlorophyll=parse_double(ifelse(Chlorophyll%in%c("<0.05", "<0.5"), 0, Chlorophyll)),
                     Secchi_depth=parse_double(Secchi_depth),
                     Temperature=parse_double(Temperature),
                     Conductivity=parse_double(Conductivity),
                     Oxygen=parse_double(Oxygen),
                     Depth=parse_double(Depth)))%>%
  mutate(MonthYear=floor_date(Date, unit = "month"),
         Year=year(Date),
         Salinity=((0.36966/(((Conductivity*0.001)^(-1.07))-0.00074))*1.28156))%>%
  select(-Conductivity)%>%
  group_by(MonthYear, Year, Station)%>%
  summarise_at(vars(c("Chlorophyll", "Secchi_depth", "Temperature", "Salinity", "Oxygen", "Depth")), ~mean(., na.rm=T))%>%
  ungroup()

Biv<-read_excel("Bivalves app/1975-18 CPUE bivalves only, 2019Sept9.xlsx",
                sheet = "75-18 CPUE per m2", skip=1)%>%
  select(Date, Station=StationCode, `Potamocorbula amurensis`, `Corbicula fluminea`)%>%
  gather(key="Taxa", value="CPUE", -Station, -Date)%>%
  mutate(Year=year(Date),
         MonthYear=floor_date(Date, unit = "month"))%>%
  separate(Station, into=c("Station", "Position"), sep="-")%>%
  group_by(Year, MonthYear, Taxa, Station)%>%
  summarise(CPUE=mean(CPUE, na.rm=T))%>%
  ungroup()%>%
  left_join(WQ, by=c("Year", "MonthYear", "Station"))


Stations<-read_excel("Bivalves app/1975-18 CPUE bivalves only, 2019Sept9.xlsx",
                     sheet = "75-17 station locations", skip=1)%>%
  select(Station=Site_Code, Latitude, Longitude)