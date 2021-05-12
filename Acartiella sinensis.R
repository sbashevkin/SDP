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

# Introduced in 1993. Surveys started counting Adults in 1994 or the inception of the survey.
# Juveniles have only been counted since 2006

AS<-Zoopsynther("Taxa", Taxa="Acartiella")%>%
  filter((Taxlifestage=="Acartiella_all_Meso Adult" & Year>=1994)| (Taxlifestage=="Acartiella_all_Meso Juvenile" & Year>=2006))%>%
  mutate(Taxname=recode(Taxname, `Acartiella_all_Meso`="Acartiella"))%>%
  select(-Taxlifestage, -SizeClass, -Genus, -Family, -Order, -Class, -Phylum, -Taxatype, -Species, -Orphan, -Undersampled)%>%
  filter(!(Year==2019 & Source=="twentymm"))%>% # Ensure same version of dataset is used as Bosmina analysis
  mutate(Date_num=as.numeric(Date))%>%
  mutate(SalSurf_l = log(SalSurf),
         Julian_day = yday(Date))%>%
  drop_na(Date, SalSurf, Temperature, Latitude, Longitude, CPUE)%>%
  pivot_wider(names_from=c(Lifestage), values_from=CPUE)%>%
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


# Model Adults ------------------------------------------------------------------

iterations <- 5e3
warmup <- iterations/4

### Try a bivariate model with both life stages. 

mb2<-brm(bf(mvbind(Adult, Juvenile) ~ t2(Julian_day_s, SalSurf_l_s, Year_s, d=c(1,1,1), bs=c("cc", "cr", "cr"), k=c(13, 5, 5)) + (1|p|Clust), hu ~ s(SalSurf_l_s, bs="cr", k=5))+
           set_rescor(FALSE),
         data=filter(AS, !is.na(Adult) & !is.na(Juvenile)), family=hurdle_lognormal(),
         prior=prior(normal(0,10), class="Intercept", resp = "Adult")+
           prior(normal(0,10), class="Intercept", resp = "Juvenile")+
           prior(normal(0,5), class="b", resp = "Adult")+
           prior(normal(0,5), class="b", resp = "Juvenile")+
           prior(cauchy(0,5), class="sigma", resp = "Adult")+
           prior(cauchy(0,5), class="sigma", resp = "Juvenile"),
         chains=1, cores=1,
         iter = iterations, warmup = warmup,
         backend = "cmdstanr", threads = threading(5))
### Model fit looks bad, a few divergent transitions

mb3<-brm(bf(CPUE ~ t2(Julian_day_s, SalSurf_l_s, Year_s, Lifestage, d=c(1,1,1,1), bs=c("cc", "cr", "cr", "fs"), k=c(13, 5, 5)) + (1|Clust), hu ~ t2(SalSurf_l_s, Lifestage, bs=c("cr", "fs"), k=5)),
         data=filter(AS_long, !is.na(CPUE)), family=hurdle_lognormal(),
         prior=prior(normal(0,10), class="Intercept")+
           prior(cauchy(0,5), class="sigma"),
         chains=1, cores=1,
         iter = iterations, warmup = warmup,
         backend = "cmdstanr", threads = threading(5))
### 10 of 3750 (0.0%) transitions ended with a divergence.
### 2172 of 3750 (58.0%) transitions hit the maximum treedepth limit of 10 or 2^10-1 leapfrog steps
### Model fit looks a little better than mb2

### Juvenile data aren't fitting well.
### Probaly because people didn't start counting them until 2006

mb4<-brm(bf(Adult ~ t2(Julian_day_s, SalSurf_l_s, Year_s, d=c(1,1,1), bs=c("cc", "cr", "cr"), k=c(13, 5, 5)) + (1|Clust), hu ~ s(SalSurf_l_s, bs=c("cr"), k=5)),
         data=filter(AS, !is.na(Adult)), family=hurdle_lognormal(),
         prior=prior(normal(0,10), class="Intercept")+
           prior(normal(0,5), class="b")+
           prior(cauchy(0,5), class="sigma"),
         chains=1, cores=1, control=list(adapt_delta=0.95),
         iter = iterations, warmup = warmup,
         backend = "cmdstanr", threads = threading(5))
### No divergences or tree depth issues

### Weird pp results, might help to add another variable to hurdle, maybe year?

mb5<-brm(bf(Adult ~ t2(Julian_day_s, SalSurf_l_s, Year_s, d=c(1,1,1), bs=c("cc", "cr", "cr"), k=c(13, 5, 5)) + (1|Clust), 
            hu ~ t2(SalSurf_l_s, Julian_day_s, bs=c("cr", "cc"), k=c(5, 4))),
         data=filter(AS, !is.na(Adult)), family=hurdle_lognormal(),
         prior=prior(normal(0,10), class="Intercept")+
           prior(normal(0,5), class="b")+
           prior(cauchy(0,5), class="sigma"),
         chains=1, cores=1, control=list(adapt_delta=0.99),
         iter = iterations, warmup = warmup,
         backend = "cmdstanr", threads = threading(5))

mb5_full<-brm(bf(Adult ~ t2(Julian_day_s, SalSurf_l_s, Year_s, d=c(1,1,1), bs=c("cc", "cr", "cr"), k=c(13, 5, 5)) + (1|Clust), 
            hu ~ t2(SalSurf_l_s, Julian_day_s, bs=c("cr", "cc"), k=c(5, 4))),
         data=filter(AS, !is.na(Adult)), family=hurdle_lognormal(),
         prior=prior(normal(0,10), class="Intercept")+
           prior(normal(0,5), class="b")+
           prior(cauchy(0,5), class="sigma"),
         chains=3, cores=3, control=list(adapt_delta=0.995, max_treedepth=15),
         iter = iterations, warmup = warmup,
         backend = "cmdstanr", threads = threading(2))
# Finished with no treedepth or divergence issues 
pp(mb5_full)

## predict -----------------------------------------------------------------

AS_preds<-zoop_predict(mb5_full, filter(AS, !is.na(Adult)))

AS_salinity<-zoop_plot(AS_preds, "salinity")
AS_year<-zoop_plot(AS_preds, "year")
AS_season<-zoop_plot(AS_preds, "season")

ggsave(AS_season, file="C:/Users/sbashevkin/OneDrive - deltacouncil/Zooplankton synthesis/Species modeling/Figures/Acartiella_season.png", device="png", units = "in", width=8, height=6)

ggsave(AS_year, file="C:/Users/sbashevkin/OneDrive - deltacouncil/Zooplankton synthesis/Species modeling/Figures/Acartiella_year.png", device="png", units = "in", width=8, height=6)

ggsave(AS_salinity, file="C:/Users/sbashevkin/OneDrive - deltacouncil/Zooplankton synthesis/Species modeling/Figures/Acartiella_salinity.png", device="png", units = "in", width=8, height=6)


## Plot station intercepts -------------------------------------------------

p_intercepts<-zoop_stations(mb5_full, select(Stations_clust, Clust, Latitude, Longitude))

ggsave(p_intercepts, file="C:/Users/sbashevkin/OneDrive - deltacouncil/Zooplankton synthesis/Species modeling/Figures/Acartiella_intercepts.png", device="png", units = "in", width=9, height=7)


# Model juveniles----------------------------------------------------------

mb4_juv<-brm(bf(Juvenile ~ t2(Julian_day_s, SalSurf_l_s, Year_s, d=c(1,1,1), bs=c("cc", "cr", "cr"), k=c(13, 5, 3)) + (1|Clust), hu ~ s(SalSurf_l_s, bs=c("cr"), k=5)),
         data=filter(AS, !is.na(Juvenile)), family=hurdle_lognormal(),
         prior=prior(normal(0,10), class="Intercept")+
           prior(normal(0,5), class="b")+
           prior(cauchy(0,5), class="sigma"),
         chains=1, cores=1, control=list(adapt_delta=0.99),
         iter = iterations, warmup = warmup,
         backend = "cmdstanr", threads = threading(5))

mb5_juv<-brm(bf(Juvenile ~ t2(Julian_day_s, SalSurf_l_s, Year_s, d=c(1,1,1), bs=c("cc", "cr", "cr"), k=c(13, 5, 3)) + (1|Clust), 
                hu ~ t2(SalSurf_l_s, Year_s, d=c(1,1), bs=c("cr", "cr"), k=c(5, 3))),
             data=filter(AS, !is.na(Juvenile)), family=hurdle_lognormal(),
             prior=prior(normal(0,10), class="Intercept")+
               prior(normal(0,5), class="b")+
               prior(cauchy(0,5), class="sigma"),
             chains=1, cores=1, control=list(adapt_delta=0.99, max_treedepth=15),
             iter = iterations, warmup = warmup,
             backend = "cmdstanr", threads = threading(5))
# No divergence or treedepth issues
# 6.7 hours
# But doesn't fix weird pp_check pattern



pp_check(mb5_juv, type="error_scatter_avg_vs_x", x="Julian_day_s")
pp_check(mb5_juv, type="error_scatter_avg_vs_x", x="Year_s") # no year pattern
pp_check(mb4_juv, type="error_scatter_avg_vs_x", x="Year_s") # not even in the model without year in the hu formula
# Looks like lots of error in the middle of the year, so trying with Julian day in the hu formula instead of year

mb6_juv<-brm(bf(Juvenile ~ t2(Julian_day_s, SalSurf_l_s, Year_s, d=c(1,1,1), bs=c("cc", "cr", "cr"), k=c(13, 5, 3)) + (1|Clust), 
                hu ~ t2(SalSurf_l_s, Julian_day_s, d=c(1,1), bs=c("cr", "cc"), k=c(5, 4))),
             data=filter(AS, !is.na(Juvenile)), family=hurdle_lognormal(),
             prior=prior(normal(0,10), class="Intercept")+
               prior(normal(0,5), class="b")+
               prior(cauchy(0,5), class="sigma"),
             chains=1, cores=1, control=list(adapt_delta=0.99, max_treedepth=15),
             iter = iterations, warmup = warmup,
             backend = "cmdstanr", threads = threading(5))
# No divergence or treedepth issues
pp_check(mb6_juv)+scale_x_log10()
# pp_check plot looks way better

mb6_juv_full<-brm(bf(Juvenile ~ t2(Julian_day_s, SalSurf_l_s, Year_s, d=c(1,1,1), bs=c("cc", "cr", "cr"), k=c(13, 5, 3)) + (1|Clust), 
                hu ~ t2(SalSurf_l_s, Julian_day_s, d=c(1,1), bs=c("cr", "cc"), k=c(5, 4))),
             data=filter(AS, !is.na(Juvenile)), family=hurdle_lognormal(),
             prior=prior(normal(0,10), class="Intercept")+
               prior(normal(0,5), class="b")+
               prior(cauchy(0,5), class="sigma"),
             chains=3, cores=3, control=list(adapt_delta=0.99, max_treedepth=15),
             iter = iterations, warmup = warmup,
             backend = "cmdstanr", threads = threading(2))
# 1 divergent transition
pp(mb6_juv_full)
## predict -----------------------------------------------------------------

AS_juv_preds<-zoop_predict(mb6_juv_full, filter(AS, !is.na(Juvenile)))

AS_juv_salinity<-zoop_plot(AS_juv_preds, "salinity")
AS_juv_year<-zoop_plot(AS_juv_preds, "year")
AS_juv_season<-zoop_plot(AS_juv_preds, "season")

ggsave(AS_juv_season, file="C:/Users/sbashevkin/OneDrive - deltacouncil/Zooplankton synthesis/Species modeling/Figures/Acartiella_juv_season.png", device="png", units = "in", width=8, height=6)

ggsave(AS_juv_year, file="C:/Users/sbashevkin/OneDrive - deltacouncil/Zooplankton synthesis/Species modeling/Figures/Acartiella_juv_year.png", device="png", units = "in", width=8, height=6)

ggsave(AS_juv_salinity, file="C:/Users/sbashevkin/OneDrive - deltacouncil/Zooplankton synthesis/Species modeling/Figures/Acartiella_juv_salinity.png", device="png", units = "in", width=8, height=6)


## Plot station intercepts -------------------------------------------------

p_juv_intercepts<-zoop_stations(mb6_juv_full, select(Stations_clust, Clust, Latitude, Longitude))

ggsave(p_juv_intercepts, file="C:/Users/sbashevkin/OneDrive - deltacouncil/Zooplankton synthesis/Species modeling/Figures/Acartiella_juv_intercepts.png", device="png", units = "in", width=9, height=7)
