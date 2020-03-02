require(sf)
require(nhdplusTools)

lon <- -122.396703
lat <- 38.046869

start_point <- st_sfc(st_point(c(lon, lat)),
                          crs = 4326)
start_comid <- discover_nhdplus_id(start_point)

discover_nldi_sources()$source

nldi_feature <- list(featureSource = "comid", featureID = start_comid)

discover_nldi_navigation(nldi_feature)

flowline_nldi <- navigate_nldi(nldi_feature, 
                               mode = "upstreamMain", 
                               data_source = "")

plot(st_geometry(flowline_nldi))

output_file_download <- file.path(tempdir(), "subset_download.gpkg")

output_file_download <-subset_nhdplus(comids = flowline_nldi$nhdplus_comid,
                                      output_file = output_file_download,
                                      nhdplus_data = "download", return_data = TRUE)

sf::st_layers(output_file_download)

flowline_download <- sf::read_sf(file.path(tempdir(), "subset_download.gpkg"), 
                                 "NHDFlowline_Network")

plot(sf::st_geometry(flowline_download), 
     lwd = 3, col = "red")

#Using HR

hr_urls <- download_nhdplushr("~/nhdplusHR", c("1805", "1804", "1802"), download_files = TRUE)

hr<-get_nhdplushr(hr_urls, layers=NULL)
hr<-get_nhdplushr(hr_urls, "~/hr_data.gpkg", layers=NULL)
