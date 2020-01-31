require(gdistance)
require(rgdal)
require(raster)
require(fasterize)
require(maptools)
require(sf)
require(tibble)
require(dplyr)
require(rlang)

pointmover <- function(Data, Attributecol, Shapefile){
  Attributecol <- sym(Attributecol)
  Attributecol <- enquo(Attributecol)
  mappoints<-st_cast(Shapefile, "POINT", warn=F)
  badpoints<-filter(Data, is.na(!!Attributecol))
  new<-mappoints[st_nearest_feature(badpoints, mappoints),]
  
  badpoints$geometry<-new$geometry
  out<-rbind(badpoints, filter(Data, !is.na(!!Attributecol)))
  return(out)
  
}

# read in the stations table (just a csv of it, to avoid DB stuff)

Stations <- read.csv("distance/20mm_stations_table.csv")

# restrict to stations included

stations_in <- c("405", "411", "418", "501", "504", "508", "513", "519", 
                 "520", "602", "606", "609", "610", "703", "704", "705", 
                 "706", "707", "711", "716", "801", "804", "809", "812", 
                 "815", "901", "902", "906", "910", "912", "914", "915", 
                 "918", "919")

Stations <- Stations%>%
  filter(Station %in% stations_in)%>%
  mutate(Lat = LatD + LatM / 60 + LatS / 3600,
         Lon = -1 * (LonD + LonM / 60 + LonS / 3600))

Waterdist <- function(Shapefile_path = "~/ZoopSynth/Old data and code/DeltaShapefile", 
                      Shapefile_crs = "+proj=utm +zone=10 +datum=NAD83",
                      Points = Stations,
                      Attributecol = "Shape_Area"){
  
  Attributecol2<-sym(Attributecol)
  Attributecol2<-enquo(Attributecol2)
  
  # read in the shape file
  map<-read_sf(Shapefile_path)%>%
    st_transform(crs=Shapefile_crs)
  
  stationlocs <- project(as.matrix(select(Points, Lon, Lat)), "+proj=utm +zone=10 ellps=WGS84")
  
  Points <- as_tibble(stationlocs)%>%
    mutate(Station = Points$Station)%>%
    st_as_sf(coords = c("Lon", "Lat"),
             crs = "+proj=utm +zone=10 ellps=WGS84",
             remove=T)
  
  # Are all points in the water polygon? (st_intersects returns a null object when there is no intersection)
  
  #Join the points to the polygon to identify the land-based point(s)
  
  Points_joined<-st_join(Points, map, join = st_intersects)%>%
    arrange(Station)
  
  
  
  # If all points are not within polygon, replace point outside polygon with closest point within polygon.
  if(!all(!is.na(pull(Points_joined, !!Attributecol2)))){
    Points_joined<-pointmover(Points_joined, Attributecol, map)%>%
      arrange(Station)
  }
  
  # Are all points in the water polygon now?
  
  if(!(length(unlist(st_intersects(Points_joined, map)))==nrow(Points))){
    stop("Points could not be moved within shapefile.")
  }
  
  # rasterize the polygon and designate water = 1, land = 0
  # using 75 m x 75 m grid squares and rasterizing the extent of the map
  
  mapextent<-st_bbox(map)
  
  cols <- round((mapextent["xmax"] - mapextent["xmin"]) / 75) # 1394 columns
  rows <- round((mapextent["ymax"] - mapextent["ymin"]) / 75) # 1500 rows
  
  r <- raster(ncol = cols, nrow = rows)
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
  
  return(waterDist)
}


