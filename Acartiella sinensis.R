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


# Prepare data ------------------------------------------------------------

# Cluster nearby stations
Stations<-zooper::stations%>%
  mutate(Source=recode(Source, twentymm="20mm"),
    Station=paste(Source, Station))%>%
  filter(Source!="YBFMP")%>%
  bind_rows(zooper::stationsEMPEZ%>%
              mutate(Station=paste("EMP", Station, year(Date), month(Date), day(Date))))%>%
  drop_na(Latitude, Longitude)

Stations_clust <- Pointcluster(Points=Stations, Distance=1000, Latitude_column=Latitude, 
                               Longitude_column=Longitude, PointID_column=Station, Expand = FALSE)

Stations_final <- unnest(Stations_clust, Station)


## Make data

## Still some NAS?

AS<-Zoopsynther("Taxa", Taxa="Acartiella")%>%
  filter(Taxlifestage%in%c("Acartiella_all_Meso Adult", "Acartiella_all_Meso Juvenile"))%>%
  mutate(Taxname=recode(Taxname, `Acartiella_all_Meso`="Acartiella"))%>%
  select(-Taxlifestage, -SizeClass, -Genus, -Family, -Order, -Class, -Phylum, -Taxatype, -Species, -Orphan, -Undersampled)%>%
  filter(!(Year==2019 & Source=="twentymm"))%>% # Ensure same version of dataset is used as Bosmina analysis
  mutate(Date_num=as.numeric(Date))%>%
  mutate(SalSurf_l = log(SalSurf),
         Julian_day = yday(Date))%>%
  drop_na(Date, SalSurf, Temperature, Latitude, Longitude, CPUE)%>%
  pivot_wider(names_from=c(Lifestage), values_from=CPUE)%>%
  filter(!is.na(Adult) & !is.na(Juvenile))%>%
  mutate(Station=if_else(Station%in%unique(zooper::stationsEMPEZ$Station), paste(Source, Station, year(Date), month(Date), day(Date)), paste(Source, Station)))%>%
  left_join(Stations_final%>%
              select(Clust, Station),
            by="Station")%>%
  mutate(Date2=if_else(is.na(Datetime), parse_date_time(paste0(year(Date), "-", month(Date), "-", day(Date), " 12:00"), "%Y-%m-%d %H:%M", tz="America/Los_Angeles"), Datetime))%>%
  group_by_at(vars(-Adult, -Juvenile, -Volume, -SampleID))%>%
  summarize(Adult=mean(Adult), Juvenile=mean(Juvenile), Volume=mean(Volume), .groups="drop")%>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326, remove=F)%>%
  st_transform(crs=st_crs(deltamapr::R_EDSM_Subregions_Mahardja))%>%
  st_join(deltamapr::R_EDSM_Subregions_Mahardja%>%
            select(SubRegion))%>%
  st_drop_geometry()%>%
  mutate(Adult_Count=round(Adult*Volume),
         Juvenile_Count=round(Juvenile*Volume),
         Year_fac=factor(Year),
         Month=month(Date, label = T),
         ID=1:nrow(.),
         Clust2=Clust)%>%
  mutate_at(vars(Date_num, SalSurf, SalSurf_l, Temperature, Latitude, Longitude, Julian_day, Year), list(s=~(.-mean(., na.rm=T))/sd(., na.rm=T)))

Sal_cat<-quantile(AS$SalSurf, probs=seq(0, 1, by=0.1))

Sal_cats<-tibble(percentile=names(Sal_cat), Min=lag(Sal_cat), Max=Sal_cat)%>%
  drop_na()%>%
  mutate(Min=if_else(Min==min(Min), 0, Min),
         Name=paste0(percentile, ": ", round(Min, 4), "-", round(Max, 4)))

AS<-AS%>%
  rowwise()%>%
  mutate(Sal_cat=filter(Sal_cats, Min<SalSurf & Max>=SalSurf)$Name)%>%
  ungroup()%>%
  mutate(Sal_cat=factor(Sal_cat, levels=Sal_cats$Name),
         Sal_cat_int=as.integer(Sal_cat))

AS_long<-pivot_longer(AS, c(Adult, Juvenile), names_to="Lifestage", values_to="CPUE")%>%
  mutate(Lifestage=factor(Lifestage))

# Plots -------------------------------------------------------------------



ggplot(AS, aes(x=Year, y=Adult))+
  geom_point()+
  facet_wrap(~Month)

ggplot(AS, aes(x=Year, y=Juvenile))+
  geom_point()+
  facet_wrap(~Month)

ggplot(AS, aes(x=SalSurf_l, y=Adult))+
  geom_point()+
  facet_wrap(~Month)

ggplot(AS, aes(x=SalSurf_l, y=Juvenile))+
  geom_point()+
  facet_wrap(~Month)


# Models ------------------------------------------------------------------

iterations <- 5e3
warmup <- iterations/4

# Try a bivariate model with both life stages. 

mb2<-brm(bf(mvbind(Adult, Juvenile) ~ t2(Julian_day_s, SalSurf_l_s, Year_s, d=c(1,1,1), bs=c("cc", "cr", "cr"), k=c(13, 5, 5)) + (1|p|Clust), hu ~ s(SalSurf_l_s, bs="cr", k=5))+
           set_rescor(FALSE),
         data=AS, family=hurdle_lognormal(),
         prior=prior(normal(0,10), class="Intercept", resp = "Adult")+
           prior(normal(0,10), class="Intercept", resp = "Juvenile")+
           prior(normal(0,5), class="b", resp = "Adult")+
           prior(normal(0,5), class="b", resp = "Juvenile")+
           prior(cauchy(0,5), class="sigma", resp = "Adult")+
           prior(cauchy(0,5), class="sigma", resp = "Juvenile"),
         chains=1, cores=1,
         iter = iterations, warmup = warmup,
         backend = "cmdstanr", threads = threading(5))
# Model fit looks bad, a few divergent transitions

mb3<-brm(bf(CPUE ~ t2(Julian_day_s, SalSurf_l_s, Year_s, Lifestage, d=c(1,1,1,1), bs=c("cc", "cr", "cr", "fs"), k=c(13, 5, 5)) + (1|Clust), hu ~ t2(SalSurf_l_s, Lifestage, bs=c("cr", "fs"), k=5)),
         data=AS_long, family=hurdle_lognormal(),
         prior=prior(normal(0,10), class="Intercept")+
           prior(cauchy(0,5), class="sigma"),
         chains=1, cores=1,
         iter = iterations, warmup = warmup,
         backend = "cmdstanr", threads = threading(5))
# 10 of 3750 (0.0%) transitions ended with a divergence.
# 2172 of 3750 (58.0%) transitions hit the maximum treedepth limit of 10 or 2^10-1 leapfrog steps