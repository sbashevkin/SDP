require(gdistance)
require(rgdal)
require(raster)
require(fasterize)
require(maptools)
require(sf)
require(tibble)
require(dplyr)
require(rlang)

# read in the shape file
map<-read_sf("Map/DeltaMap.shp")%>%
  st_transform(crs="+proj=utm +zone=10 +datum=NAD83")

# read in the stations table (just a csv of it, to avoid DB stuff)

Stations <- read.csv("20mm_stations_table.csv")

# restrict to stations included

stations_in <- c("405", "411", "418", "501", "504", "508", "513", "519", 
                 "520", "602", "606", "609", "610", "703", "704", "705", 
                 "706", "707", "711", "716", "801", "804", "809", "812", 
                 "815", "901", "902", "906", "910", "912", "914", "915", 
                 "918", "919")

Stations <- Stations[which(Stations$Station %in% stations_in), ]

# convert DMS to UTM for stations

Lat <- rep(NA, nrow(Stations))
Lon <- rep(NA, nrow(Stations))
    
for(i in 1:nrow(Stations)){
  lat <- Stations$LatD[i] + Stations$LatM[i] / 60 + Stations$LatS[i] / 3600
  lon <- Stations$LonD[i] + Stations$LonM[i] / 60 + Stations$LonS[i] / 3600
  Lat[i] <- lat
  Lon[i] <- -1 * lon
}

xy <- cbind(Lon, Lat)
stationlocs <- project(xy, "+proj=utm +zone=10 ellps=WGS84")
Easting <- stationlocs[,1]
Northing <- stationlocs[,2]

Points <- as_tibble(stationlocs)%>%
  mutate(Station = Stations$Station)%>%
  st_as_sf(coords = c("Lon", "Lat"),
                   crs = "+proj=utm +zone=10 ellps=WGS84",
                   remove=T)

# Are all points in the water polygon? (st_intersects returns a null object when there is no intersection)

length(unlist(st_intersects(Points, map)))==nrow(Points)

#Join the points to the polygon to identify the land-based point(s)

Points_joined<-st_join(Points, map, join = st_intersects)%>%
  arrange(Station)

# Replace point outside polygon with closest point within polygon. Hopefully a generalizable approach

pointreplacer <- function(data, attributecol, shapefile){
  attributecol <- enquo(attributecol)
  mappoints<-st_cast(shapefile, "POINT", warn=F)
  badpoints<-filter(data, is.na(!!attributecol))
  new<-mappoints[st_nearest_feature(badpoints, mappoints),]
  
  badpoints$geometry<-new$geometry
  out<-rbind(badpoints, filter(data, !is.na(!!attributecol)))
  return(out)
  
}

Points_joined<-pointreplacer(Points_joined, TYPE, map)%>%
  arrange(Station)

# Are all points in the water polygon now?

length(unlist(st_intersects(Points_joined, map)))==nrow(Points)

# rasterize the polygon and designate water = 1, land = 0
# using 75 m x 75 m grid squares and rasterizing the extent of the map

mapextent<-st_bbox(map)

(mapextent["ymax"] - mapextent["ymin"]) / 75 # 1500 rows
(mapextent["xmax"] - mapextent["xmin"]) / 75 # 1394 columns

r <- raster(ncol = 1394, nrow = 1500)
extent(r) <- extent(c(mapextent["xmin"], mapextent["xmax"], mapextent["ymin"], mapextent["ymax"]))
rp <- fasterize(map, r)
rp[is.na(rp)] <- 0

# measure distances between points within the polygon (i.e. water distances)
# using the mean function in 16 directions (knight and one-cell queen moves) 
# first need to measure transitions between grid squares
# requires geographic correction 

tp <- transition(rp, mean, 16)
tpc <- geoCorrection(tp, "c", scl = FALSE)
waterDist <- costDistance(tpc, st_coordinates(Points_joined))
