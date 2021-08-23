require(dplyr)
require(lubridate)
require(tidyr)
require(ggplot2)
require(zooper)
require(mgcv)
require(sf)
require(brms)
require(spacetools)
require(geoR)
require(corpcor)
require(gstat)
require(spacetime)
require(colorspace)
require(stringr)
source("Utility functions.R")

# Cluster nearby stations
Stations<-zooper::stations%>%
  mutate(Station=paste(Source, Station))%>%
  filter(Source!="YBFMP")%>%
  bind_rows(zooper::stationsEMPEZ%>%
              mutate(Station=paste("EMP", Station, year(Date), month(Date), day(Date))))%>%
  drop_na(Latitude, Longitude)

Stations_clust <- Pointcluster(Points=Stations, Distance=1000, Latitude_column=Latitude, 
                               Longitude_column=Longitude, PointID_column=Station, Expand = FALSE)

Stations_final <- unnest(Stations_clust, Station)

#Get in-water distance between stations
distance <- Waterdist(Water_map=spacetools::Delta, Points=Stations_clust, Latitude_column=Latitude, 
                      Longitude_column=Longitude, PointID_column=Clust)

AV<-zoopComb%>%
  filter(Taxname=="Acanthocyclops vernalis")%>%
  left_join(zoopEnvComb%>%
              select(-Source),
            by="SampleID")%>%
  filter(!(Year==2019 & Source=="twentymm"))%>% # Ensure same version of dataset is used as Bosmina analysis
  mutate(Date_num=as.numeric(Date))%>%
  mutate(SalSurf_l = log(SalSurf),
         Julian_day = yday(Date))%>%
  drop_na(Date, SalSurf, Temperature, Latitude, Longitude, CPUE)%>%
  mutate(Station=if_else(Station%in%unique(zooper::stationsEMPEZ$Station), paste(Source, Station, year(Date), month(Date), day(Date)), paste(Source, Station)))%>%
  left_join(Stations_final%>%
              select(Clust, Station),
            by="Station")%>%
  mutate(Date2=if_else(is.na(Datetime), parse_date_time(paste0(year(Date), "-", month(Date), "-", day(Date), " 12:00"), "%Y-%m-%d %H:%M", tz="America/Los_Angeles"), Datetime))%>%
  group_by_at(vars(-CPUE, -Volume, -SampleID))%>%
  summarize(CPUE=mean(CPUE), Volume=mean(Volume), .groups="drop")%>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=F)%>%
  st_transform(crs=st_crs(deltamapr::R_EDSM_Subregions_Mahardja))%>%
  st_join(deltamapr::R_EDSM_Subregions_Mahardja%>%
            select(SubRegion))%>%
  st_drop_geometry()%>%
  mutate(Count=round(CPUE*Volume),
         Year_fac=factor(Year),
         Month=month(Date, label = T),
         ID=1:nrow(.),
         Clust2=Clust)%>%
  mutate_at(vars(Date_num, SalSurf, SalSurf_l, Temperature, Latitude, Longitude, Julian_day, Year), list(s=~(.-mean(., na.rm=T))/sd(., na.rm=T)))

Sal_cat<-quantile(AV$SalSurf, probs=seq(0, 1, by=0.1))

Sal_cats<-tibble(percentile=names(Sal_cat), Min=lag(Sal_cat), Max=Sal_cat)%>%
  drop_na()%>%
  mutate(Min=if_else(Min==min(Min), 0, Min),
         Name=paste0(percentile, ": ", round(Min, 4), "-", round(Max, 4)))

AV<-AV%>%
  rowwise()%>%
  mutate(Sal_cat=filter(Sal_cats, Min<SalSurf & Max>=SalSurf)$Name)%>%
  ungroup()%>%
  mutate(Sal_cat=factor(Sal_cat, levels=Sal_cats$Name),
         Sal_cat_int=as.integer(Sal_cat))

distance2<-distance/max(distance)
distance3<-distance2[which(rownames(distance2)%in%unique(AV$Clust)), which(colnames(distance2)%in%unique(AV$Clust))]

distance_cov <-cov.spatial(distance3, cov.pars=c(1,2))
distance_cov <- make.positive.definite(distance_cov)
# Plot change over time
#AV_sum<-AV%>%

ggplot(AV, aes(x=Year, y=CPUE))+
  geom_point()+
  facet_wrap(~Month)

ggplot(AV, aes(x=SalSurf_l, y=CPUE))+
  geom_point()+
  facet_wrap(~Month)

ggplot(AV, aes(x=Temperature, y=CPUE))+
  geom_point()+
  facet_wrap(~Month)

AV_spatial_sum<-AV%>%
  group_by(SubRegion)%>%
  summarise(CPUE_mean=mean(CPUE), CPUE_sd=sd(CPUE), prop_zero=length(which(CPUE==0))/n(), N_stations=n_distinct(Station), N=n(), .groups="drop")

AV_spatial<-deltamapr::R_EDSM_Subregions_Mahardja%>%
  select(SubRegion)%>%
  left_join(AV_spatial_sum, by="SubRegion")

ggplot()+
  geom_sf(data=AV_spatial, aes(fill=log(CPUE_mean+1)))+
  scale_fill_viridis_c()

ggplot(AV_spatial_sum, aes(x=Month, y=CPUE_mean, ymax=CPUE_mean+CPUE_sd, ymin=CPUE_mean))+
  geom_pointrange()+
  facet_wrap(~SubRegion)

ggplot(AV, aes(x=Month, y=CPUE))+
  geom_point()+
  facet_wrap(~Sal_cat)


# models ------------------------------------------------------------------
iterations <- 5e3
warmup <- iterations/4

mb<-brm(bf(CPUE ~ poly(Julian_day_s, 2)*Sal_cat + (1|Year_fac) + (1|Clust), hu ~ Sal_cat),
          data=AV, family=hurdle_lognormal(),
          prior=prior(normal(0,10), class="Intercept")+
            prior(normal(0,5), class="b")+
            prior(cauchy(0,5), class="sigma"),
          chains=1, cores=1,
          iter = iterations, warmup = warmup,
          backend = "cmdstanr", threads = threading(4))
mb<-add_criterion(mb, c("loo", "waic"))

mb2<-brm(bf(CPUE ~ t2(Julian_day_s, SalSurf_l_s, Year_s, d=c(1,1,1), bs=c("cc", "cr", "cr"), k=c(13, 5, 5)) + (1|Clust), hu ~ s(SalSurf_l_s, bs="cr", k=5)),
          data=AV, family=hurdle_lognormal(),
          prior=prior(normal(0,10), class="Intercept")+
            prior(normal(0,5), class="b")+
            prior(cauchy(0,5), class="sigma"),
          chains=1, cores=1, control=list(adapt_delta=0.995, max_treedepth=15),
          iter = iterations, warmup = warmup,
          backend = "cmdstanr", threads = threading(4))
# No warnings
mb2<-add_criterion(mb2, c("loo", "waic"), moment_match = TRUE)

mb2_check<-pp(mb2)

mb2_full<-brm(bf(CPUE ~ t2(Julian_day_s, SalSurf_l_s, Year_s, d=c(1,1,1), bs=c("cc", "cr", "cr"), k=c(13, 5, 5)) + (1|Clust), hu ~ s(SalSurf_l_s, bs="cr", k=5)),
         data=AV, family=hurdle_lognormal(),
         prior=prior(normal(0,10), class="Intercept")+
           prior(normal(0,5), class="b")+
           prior(cauchy(0,5), class="sigma"),
         chains=3, cores=3, control=list(adapt_delta=0.995, max_treedepth=15),
         iter = iterations, warmup = warmup,
         backend = "cmdstanr", threads = threading(2))
# 1 divergent transition, which shouldn't be a problem https://stats.stackexchange.com/questions/432479/divergent-transitions-in-stan


mb3<-brm(bf(CPUE ~ t2(Julian_day_s, SalSurf_l_s, Year_s, d=c(1,1,1), bs=c("cc", "cr", "cr"), k=c(13, 5, 5)) + (1|gr(Clust, cov=distance_cov)), hu ~ s(SalSurf_l_s, bs="cr", k=5)),
         data=AV, data2=list(distance_cov=distance_cov),
         family=hurdle_lognormal(),
         prior=prior(normal(0,10), class="Intercept")+
           prior(normal(0,5), class="b")+
           prior(cauchy(0,5), class="sigma"),
         chains=1, cores=1, control=list(adapt_delta=0.95, max_treedepth=15),
         iter = iterations, warmup = warmup,
         backend = "cmdstanr", threads = threading(5))

 # Model checks ------------------------------------------------------------

pp(mb2_full)
mb2_full_vario<-zoop_vario(model=mb2_full, data=AV, yvar="CPUE", cores=5)
mb2_full_vario_plot<-zoop_vario_plot(mb2_full_vario)
ggsave(mb2_full_vario_plot, filename="figures/Acanthocyclops_variogram.png", device="png", width=8, height=5, units="in")

mb3_vario<-zoop_vario(model=mb3, data=AV, yvar="CPUE", cores=5)
mb3_vario_plot<-zoop_vario_plot(mb3_vario)
# predict -----------------------------------------------------------------

Data_effort<-AV%>%
  group_by(Month, Year)%>%
  summarise(N=n(), .groups="drop")%>%
  mutate(Month=as.integer(Month))

AV_preds<-zoop_predict(mb2_full, AV, confidence=99)%>%
  left_join(Data_effort, by=c("Month", "Year"))%>%
  mutate(across(c(Pred, lowerCI, upperCI), ~if_else(is.na(N), NA_real_, .x)))

AV_salinity<-zoop_plot(AV_preds, "salinity")
AV_year<-zoop_plot(AV_preds, "year")
AV_season<-zoop_plot(AV_preds, "season")

ggsave(AV_season, file="figures/Acanthocyclops_season.png", device="png", units = "in", width=8, height=6)

ggsave(AV_year, file="figures/Acanthocyclops_year.png", device="png", units = "in", width=8, height=6)

ggsave(AV_salinity, file="figures/Acanthocyclops_salinity.png", device="png", units = "in", width=8, height=6)


# Plot station intercepts -------------------------------------------------

p_intercepts<-zoop_stations(mb2_full, select(Stations_clust, Clust, Latitude, Longitude))

ggsave(p_intercepts, file="figures/Acanthocyclops_intercepts.png", device="png", units = "in", width=6, height=4)
