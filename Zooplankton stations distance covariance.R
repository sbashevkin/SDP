require(tidyverse)
require(zooper)
require(sp)
require(rgdal)
require(spacetools)
require(sf)
require(leaflet)
require(geoR)

# Clustering distance (m)
d <- 1000

Stations<-zooper::stations%>%
  mutate(ID=paste(Source, Station))%>%
  filter(Source!="YBFMP")%>%
  drop_na(Latitude, Longitude)

xy <- Stations%>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs="+proj=longlat +datum=NAD83")%>%
  st_transform(crs="+init=epsg:3488 +datum=NAD83")

chc <- hclust(dist(st_coordinates(xy)), method="complete")

# Compute distance
chc.d40 <- cutree(chc, h=d)

# Join results to meuse sp points
xy <- xy%>%
  mutate(Clust = chc.d40)%>%
  st_drop_geometry()%>%
  select(ID, Clust)

Stations_clust<-Stations%>%
  left_join(xy, by="ID")%>%
  group_by(Clust)%>%
  summarise(Latitude = mean(Latitude), Longitude = mean(Longitude), Station = list(ID))

#Get in-water distance between stations
map<-read_sf("~/Hydro_CH2MHill_Merged/hydro_merged_dslv.shp")
distance <- Waterdist(map, Stations_clust, Latitude, Longitude, Clust)

map2 <- sf::st_union(map)%>%
  sf::st_as_sf()%>%
  dplyr::mutate(Inside=TRUE)%>%
  dplyr::rename(geometry = .data$x)%>%
  st_transform(crs="+proj=utm +zone=10 ellps=WGS84")

Stations <- Stations_clust%>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326)%>%
  st_transform(crs="+proj=utm +zone=10 ellps=WGS84")%>%
  st_join(map2, join=st_intersects)%>%
  Pointmover(., Inside, map2)

#Transfrom to crs preferred by leaflet

Stations <- st_transform(Stations, crs=4326)
map2 <- st_transform(map2, crs=4326)

#Which are showing infinite distances?

badpoints<-as.integer(names(which(!is.finite(distance[1,]))))

leaflet()%>%
  addProviderTiles("Esri.WorldGrayCanvas")%>%
  addPolygons(data=map2)%>%
  addCircleMarkers(data=filter(Stations, !(Clust%in%badpoints)), weight = 1, color = "black", 
                   radius = 4, label = ~Clust)%>%
  addCircleMarkers(data=filter(Stations, Clust%in%badpoints), weight = 1, color = "red", 
                   radius = 4, label = ~Clust, labelOptions = list(permanent=TRUE))

distance<-distance/max(distance)

distance_cov <-cov.spatial(distance, cov.pars=c(1,1))
