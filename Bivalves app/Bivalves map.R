require(leaflet)
require(leaflet.minicharts)
require(readxl)
require(lubridate)
require(tidyverse)

Biv<-read_excel("Bivalves app/1975-18 CPUE bivalves only, 2019Sept9.xlsx",
                sheet = "75-18 CPUE per m2", skip=1)%>%
  select(Date, Station=StationCode, `Potamocorbula amurensis`, `Corbicula fluminea`)%>%
  gather(key="Taxa", value="CPUE", -Station, -Date)%>%
  mutate(Year=year(Date),
         MonthYear=floor_date(Date, unit = "month"))%>%
  separate(Station, into=c("Station", "Position"), sep="-")

Stations<-read_excel("Bivalves app/1975-18 CPUE bivalves only, 2019Sept9.xlsx",
                     sheet = "75-17 station locations", skip=1)%>%
  select(Station=Site_Code, Latitude, Longitude)

Bivsum<-Biv%>%
  left_join(Stations, by="Station")%>%
  group_by(Year, Taxa, Latitude, Longitude)%>%
  summarise(CPUE=mean(CPUE, na.rm=T))%>%
  ungroup()

mapdata<-Bivsum%>%
  spread(key=Taxa, value=CPUE)%>%
  mutate(Total=`Corbicula fluminea`+`Potamocorbula amurensis`)

colors <- c("#1b9e77", "#7570b3")

leaflet(data = mapdata)%>%
  addProviderTiles("Esri.WorldGrayCanvas")%>%
  addMinicharts(lng = mapdata$Longitude, lat = mapdata$Latitude,
                type = "pie",
                chartdata = mapdata%>%select_at(vars(unique(Bivsum$Taxa)))%>%as.matrix(), 
                colorPalette = colors, transitionTime = 0, opacity=0.5, width=60*(sqrt(mapdata$Total)/sqrt(max(mapdata$Total))), time=mapdata$Year)
