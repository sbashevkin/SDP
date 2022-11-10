require(tidyverse)
require(ggspatial)
require(sf)
require(maps)

Stations<-zooper::stations%>%
  filter(Source!="YBFMP")%>%
  mutate(Source=recode(Source, twentymm="20mm"),
         Type="Fixed")%>%
  bind_rows(zooper::stationsEMPEZ%>%
              mutate(Type="EZ",
                     Source="EMP")%>%
              select(-Date))%>%
  drop_na()%>%
  mutate(Source=factor(Source, levels=c("EMP", "20mm", "FMWT", "STN", "FRP")))%>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=F)%>%
  st_transform(crs=26910)

base<-deltamapr::WW_Delta%>%
  st_transform(crs=26910)

labels<-tibble(label=c("San Francisco Bay", "San Pablo Bay", "Suisun Bay", "Suisun Marsh", 
                       "Confluence", "Cache Slough", "Sacramento\nRiver", "San Joaquin River", "Napa River",
                       "Sacramento Ship Channel", "Cosumnes\nRiver", "Mokelumne\nRiver", "Carquinez Strait"), 
               Latitude=c(37.9, 38.07, 38.08, 38.2, 38.046, 38.24, 38.50000, 37.9, 38.23, 38.51892, 38.35944, 38.2, 38.045), 
               Longitude=c(-122.4, -122.4, -122.05, -122.05, -121.9, -121.69, -121.5600, -121.325, -122.3, -121.588, -121.3404, -121.335, -122.189866),
               label_lat=c(37.9, 38.11, 38.15, 38.25, 38, 38.2, 38.49785, 37.85, 38.25, 38.54994, 38.43199, 38.11588, 37.961428), 
               label_lon=c(-122.25, -122.38, -122.18, -122.18, -122, -121.8, -121.4, -121.42, -122.37, -121.8, -121.4, -121.4, -122.178691))%>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326)%>%
  st_transform(crs=26910)%>%
  mutate(X=st_coordinates(.)[,"X"], Y=st_coordinates(.)[,"Y"])%>%
  st_drop_geometry()%>%
  st_as_sf(coords=c("label_lon", "label_lat"), crs=4326)%>%
  st_transform(crs=26910)%>%
  mutate(label_X=st_coordinates(.)[,"X"], label_Y=st_coordinates(.)[,"Y"])%>%
  st_drop_geometry()

states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))%>%
  st_transform(crs=st_crs(base))
california<-filter(states, ID=="california")

base2<-base%>%
  st_crop(st_bbox(Stations))

station_lims<-st_bbox(Stations)

pout<-ggplot(states)+
  geom_sf(color="dodgerblue3")+
  geom_sf(data=base2, color="dodgerblue3", fill="dodgerblue3")+
  geom_rect(xmin = station_lims["xmin"]-22000, xmax = station_lims["xmax"]+22000, ymin = station_lims["ymin"]-22000, ymax = station_lims["ymax"]+22000, 
            fill = NA, colour = "black", size = 0.7)+
  coord_sf(xlim=c(st_bbox(california)["xmin"], st_bbox(california)["xmax"]), ylim=c(st_bbox(california)["ymin"], st_bbox(california)["ymax"]))+
  theme_bw()+
  theme(panel.background = element_rect(fill = "dodgerblue3"), axis.text.x=element_text(angle=45, hjust=1))
pout

p<-ggplot() +
  geom_sf(data=base, fill="slategray3", color="slategray4")+
  geom_segment(data=labels, aes(x=label_X, y=label_Y, xend=X, yend=Y), arrow=arrow(type="closed", length=unit(0.1, "inches")), size=1)+
  geom_label(data=labels, aes(label=label, x=label_X, y=label_Y))+
  geom_sf(data=Stations, aes(fill=Source, alpha=Type, shape=Source))+
  coord_sf(xlim=c(station_lims["xmin"], station_lims["xmax"]), ylim=c(station_lims["ymin"], station_lims["ymax"]))+
  scale_fill_brewer(type="qual", palette="Set1", name="Survey", 
                    guide=guide_legend(title.position = "top", title.hjust = 0.5))+
  scale_shape_manual(values=21:25, name="Survey")+
  scale_alpha_manual(values=c(0.2, 0.7), name="Station type", 
                     guide=guide_legend(override.aes=list(shape=19), title.position = "top", title.hjust = 0.5))+
  ylab("")+
  xlab("")+
  annotation_scale(location = "bl") +
  annotation_north_arrow(location = "bl", pad_y=unit(0.05, "npc"), which_north = "true")+
  theme_bw()+
  theme(legend.position = "bottom", legend.background=element_rect(color="black"))+
  annotation_custom(
    grob = ggplotGrob(pout),
    xmin = -Inf,
    xmax = 570000,
    ymin = 4240000,
    ymax = Inf
  )
p

ggsave("figures/Map.png", plot=p, device="png", width=8, height=8, units = "in")
