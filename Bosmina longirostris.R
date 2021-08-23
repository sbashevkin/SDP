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

# Try adding volume to the hu formula?

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

Delta_transitioned <- spacetools::Maptransitioner(spacetools::Delta)

#Get in-water distance between stations
distance <- Waterdist(Water_map=spacetools::Delta, Points=Stations_clust, Latitude_column=Latitude, 
                      Longitude_column=Longitude, PointID_column=Clust, Water_map_transitioned = Delta_transitioned)


# This was fit before 20mm 2019 data incorporated into zooper -------------


BL<-zoopComb%>%
  filter(Taxname=="Bosmina longirostris")%>%
  left_join(zoopEnvComb%>%
              select(-Source),
            by="SampleID")%>%
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
  mutate_at(vars(Date_num, SalSurf, SalSurf_l, Temperature, Latitude, Longitude, Julian_day, Year, Volume), list(s=~(.-mean(., na.rm=T))/sd(., na.rm=T)))

Sal_cat<-quantile(BL$SalSurf, probs=seq(0, 1, by=0.1))

Sal_cats<-tibble(percentile=names(Sal_cat), Min=lag(Sal_cat), Max=Sal_cat)%>%
  drop_na()%>%
  mutate(Min=if_else(Min==min(Min), 0, Min),
    Name=paste0(percentile, ": ", round(Min, 4), "-", round(Max, 4)))

BL<-BL%>%
  rowwise()%>%
  mutate(Sal_cat=filter(Sal_cats, Min<SalSurf & Max>=SalSurf)$Name)%>%
  ungroup()%>%
  mutate(Sal_cat=factor(Sal_cat, levels=Sal_cats$Name),
         Sal_cat_int=as.integer(Sal_cat))

distance2<-distance/max(distance)
distance3<-distance2[which(rownames(distance2)%in%unique(BL$Clust)), which(colnames(distance2)%in%unique(BL$Clust))]

distance_cov <-cov.spatial(distance3, cov.pars=c(1,2))
distance_cov <- make.positive.definite(distance_cov)
# Plot change over time
#BR_sum<-BL%>%
  
ggplot(BL, aes(x=Year, y=CPUE))+
  geom_point()+
  facet_wrap(~Month)

ggplot(BL, aes(x=SalSurf_l, y=CPUE))+
  geom_point()+
  facet_wrap(~Month)

ggplot(BL, aes(x=Temperature, y=CPUE))+
  geom_point()+
  facet_wrap(~Month)

BR_spatial_sum<-BL%>%
  group_by(SubRegion)%>%
  summarise(CPUE_mean=mean(CPUE), CPUE_sd=sd(CPUE), prop_zero=length(which(CPUE==0))/n(), N_stations=n_distinct(Station), N=n(), .groups="drop")

BR_spatial<-deltamapr::R_EDSM_Subregions_Mahardja%>%
  select(SubRegion)%>%
  left_join(BR_spatial_sum, by="SubRegion")

ggplot()+
  geom_sf(data=BR_spatial, aes(fill=log(CPUE_mean+1)))+
  scale_fill_viridis_c()

ggplot(BR_spatial_sum, aes(x=Month, y=CPUE_mean, ymax=CPUE_mean+CPUE_sd, ymin=CPUE_mean))+
  geom_pointrange()+
  facet_wrap(~SubRegion)

ggplot(BL, aes(x=Month, y=CPUE))+
  geom_point()+
  facet_wrap(~Sal_cat)
# Models ------------------------------------------------------------------
m1<-gam(Count ~ offset(Volume) + te(Longitude_s, Latitude_s, Julian_day_s, d=c(2,1), bs=c("cr", "cc"), k=c(20, 13)),
        family = nb, data=BL, control=list(trace=TRUE), method="REML")

m2<-gam(list(CPUE ~ te(Julian_day_s, Year, by=Sal_cat, d=c(1,1), bs=c("cc", "cr"), k=c(13, 15)), ~1), family=gammals, data=BL, 
        control=list(trace=TRUE), method="REML")

m3<-gam(list(CPUE ~ te(Julian_day_s, Year, SalSurf_l_s, bs=c("cc", "cr", "cr"), k=c(13, 15, 5)), ~s(SalSurf_l_s, k=5)), family=gammals, data=BL, 
        control=list(trace=TRUE), method="REML")

iterations <- 5e3
warmup <- iterations/4

mb<-brm(CPUE ~ t2(Julian_day_s, SalSurf_l_s, bs=c("cc", "cr"), k=c(13, 5)), 
        data=BL, family=hurdle_lognormal(),
        prior=prior(normal(0,10), class="Intercept")+
          prior(normal(0,5), class="b")+
          prior(cauchy(0,5), class="sigma"),
        chains=1, cores=1, control=list(adapt_delta=0.99),
        iter = iterations, warmup = warmup,
        backend = "cmdstanr", threads = threading(4))

mb2<-brm(CPUE ~ poly(Julian_day_s, 2)*Sal_cat + Year_fac, 
         data=BL, family=hurdle_lognormal(),
         prior=prior(normal(0,10), class="Intercept")+
           prior(normal(0,5), class="b")+
           prior(cauchy(0,5), class="sigma"),
         chains=1, cores=1,
         iter = iterations, warmup = warmup,
         backend = "cmdstanr", threads = threading(4))

mb2B<-brm(CPUE ~ poly(Julian_day_s, 2)*Sal_cat + Year_fac + (1|Clust), 
         data=BL, family=hurdle_lognormal(),
         prior=prior(normal(0,10), class="Intercept")+
           prior(normal(0,5), class="b")+
           prior(cauchy(0,5), class="sigma"),
         chains=1, cores=1,
         iter = iterations, warmup = warmup,
         backend = "cmdstanr", threads = threading(4))


mb2C<-brm(CPUE ~ poly(Julian_day_s, 2)*Sal_cat + (1|Year_fac) + (1|Clust), 
          data=BL, family=hurdle_lognormal(),
          prior=prior(normal(0,10), class="Intercept")+
            prior(normal(0,5), class="b")+
            prior(cauchy(0,5), class="sigma"),
          chains=1, cores=1,
          iter = iterations, warmup = warmup,
          backend = "cmdstanr", threads = threading(4))


mb2D<-brm(CPUE ~ -1 + poly(Julian_day_s, 2)*Sal_cat + (1|Year_fac) + (1|Clust), 
          data=BL, family=hurdle_lognormal(),
          prior=prior(normal(0,5), class="b")+
            prior(cauchy(0,5), class="sigma"),
          chains=1, cores=1,
          iter = iterations, warmup = warmup,
          backend = "cmdstanr", threads = threading(4))

mb2E<-brm(bf(CPUE ~ poly(Julian_day_s, 2)*Sal_cat + (1|Year_fac) + (1|Clust), hu ~ Sal_cat),
         data=BL, family=hurdle_lognormal(),
         prior=prior(normal(0,10), class="Intercept")+
           prior(normal(0,5), class="b")+
           prior(cauchy(0,5), class="sigma"),
         chains=1, cores=1,
         iter = iterations, warmup = warmup,
         backend = "cmdstanr", threads = threading(4))
mb2E<-add_criterion(mb2E, c("loo", "waic"))

mb2F<-brm(bf(CPUE ~ poly(Julian_day_s, 2)*Sal_cat + (1|Year_fac) + (1|Clust), hu ~ poly(SalSurf_l_s, 2)),
          data=BL, family=hurdle_lognormal(),
          prior=prior(normal(0,10), class="Intercept")+
            prior(normal(0,5), class="b")+
            prior(cauchy(0,5), class="sigma"),
          chains=1, cores=1,
          iter = iterations, warmup = warmup,
          backend = "cmdstanr", threads = threading(4))
# Weird parameter estimates and fitted vs measured values plot looks similar to mb2E so sticking with that model

mb2G<-brm(bf(CPUE ~ poly(Julian_day_s, 2)*Sal_cat + (1|Year_fac) + (1|Clust) + (1|ID), hu ~ Sal_cat),
          data=BL, family=hurdle_lognormal(),
          prior=prior(normal(0,10), class="Intercept")+
            prior(normal(0,5), class="b")+
            prior(cauchy(0,5), class="sigma"),
          chains=1, cores=1, control=list(adapt_delta=0.99),
          iter = iterations, warmup = warmup,
          backend = "cmdstanr", threads = threading(4))
# OLRE has divergent transitions even with adapt_delta=0.99

## Next try incorporating station distance into model

mb2H<-brm(bf(CPUE ~ poly(Julian_day_s, 2)*Sal_cat + (1|Year_fac) + (1|gr(Clust, cov=distance_cov)), hu ~ Sal_cat),
          data=BL, data2=list(distance_cov=distance_cov), 
          family=hurdle_lognormal(),
          prior=prior(normal(0,10), class="Intercept")+
            prior(normal(0,5), class="b")+
            prior(cauchy(0,5), class="sigma"),
          chains=1, cores=1,
          iter = iterations, warmup = warmup,
          backend = "cmdstanr", threads = threading(4))
mb2H<-add_criterion(mb2H, c("loo", "waic"))

mb2I<-brm(bf(CPUE ~ poly(Julian_day_s, 2)*Sal_cat + (1|Year_fac) + (1|gr(Clust, cov=distance_cov)) + (1|Clust2), hu ~ Sal_cat),
          data=BL, data2=list(distance_cov=distance_cov), 
          family=hurdle_lognormal(),
          prior=prior(normal(0,10), class="Intercept")+
            prior(normal(0,5), class="b")+
            prior(cauchy(0,5), class="sigma"),
          chains=1, cores=1,
          iter = iterations, warmup = warmup,
          backend = "cmdstanr", threads = threading(4))
mb2I<-add_criterion(mb2I, c("loo", "waic"))
## Very tiny loo_ic improvement from including the distance matrix, so moving ahead without it. 

# See if cubic polynomial improves fit
mb2J<-brm(bf(CPUE ~ poly(Julian_day_s, 3)*Sal_cat + (1|Year_fac) + (1|Clust), hu ~ Sal_cat),
          data=BL, family=hurdle_lognormal(),
          prior=prior(normal(0,10), class="Intercept")+
            prior(normal(0,5), class="b")+
            prior(cauchy(0,5), class="sigma"),
          chains=1, cores=1,
          iter = iterations, warmup = warmup,
          backend = "cmdstanr", threads = threading(4))
mb2J<-add_criterion(mb2J, c("loo", "waic"))

mb2K<-brm(bf(CPUE ~ s(Julian_day_s, bs="cc", k=13, by=Sal_cat) + Sal_cat + (1|Year_fac) + (1|Clust), hu ~ Sal_cat),
          data=BL, family=hurdle_lognormal(),
          prior=prior(normal(0,10), class="Intercept")+
            prior(normal(0,5), class="b")+
            prior(cauchy(0,5), class="sigma"),
          chains=1, cores=1,
          iter = iterations, warmup = warmup,
          backend = "cmdstanr", threads = threading(4))
mb2K<-add_criterion(mb2K, c("loo", "waic"))

mb2L<-brm(bf(CPUE ~ t2(Julian_day_s, SalSurf_l_s, bs=c("cc", "cr"), k=c(13, 5)) + (1|Year_fac) + (1|Clust), hu ~ s(SalSurf_l_s, bs="cr", k=5)),
          data=BL, family=hurdle_lognormal(),
          prior=prior(normal(0,10), class="Intercept")+
            prior(normal(0,5), class="b")+
            prior(cauchy(0,5), class="sigma"),
          chains=1, cores=1, control=list(adapt_delta=0.99),
          iter = iterations, warmup = warmup,
          backend = "cmdstanr", threads = threading(4))
mb2L<-add_criterion(mb2L, c("loo", "waic"))

mb2M<-brm(bf(CPUE ~ t2(Julian_day_s, SalSurf_l_s, Year_s, d=c(1,1,1), bs=c("cc", "cr", "cr"), k=c(13, 5, 5)) + (1|Clust), hu ~ s(SalSurf_l_s, bs="cr", k=5)),
          data=BL, family=hurdle_lognormal(),
          prior=prior(normal(0,10), class="Intercept")+
            prior(normal(0,5), class="b")+
            prior(cauchy(0,5), class="sigma"),
          chains=1, cores=1, control=list(adapt_delta=0.99),
          iter = iterations, warmup = warmup,
          backend = "cmdstanr", threads = threading(4))
mb2M<-add_criterion(mb2M, c("loo", "waic"))


mb2M_full<-brm(bf(CPUE ~ t2(Julian_day_s, SalSurf_l_s, Year_s, d=c(1,1,1), bs=c("cc", "cr", "cr"), k=c(13, 5, 5)) + (1|Clust), hu ~ s(SalSurf_l_s, bs="cr", k=5)),
          data=BL, family=hurdle_lognormal(),
          prior=prior(normal(0,10), class="Intercept")+
            prior(normal(0,5), class="b")+
            prior(cauchy(0,5), class="sigma"),
          chains=3, cores=3, control=list(adapt_delta=0.99, max_treedepth = 15),
          iter = iterations, warmup = warmup,
          backend = "cmdstanr", threads = threading(2))

mb2N<-brm(bf(CPUE ~ t2(Julian_day_s, SalSurf_l_s, Year_s, d=c(1,1,1), bs=c("cc", "cr", "cr"), k=c(13, 5, 5)) + (1|Clust), hu ~ s(SalSurf_l_s, bs="cr", k=5) + Volume_s),
          data=BL, family=hurdle_lognormal(),
          prior=prior(normal(0,10), class="Intercept")+
            prior(normal(0,5), class="b")+
            prior(cauchy(0,5), class="sigma"),
          chains=1, cores=1, control=list(adapt_delta=0.99),
          iter = iterations, warmup = warmup,
          backend = "cmdstanr", threads = threading(4))
mb2N<-add_criterion(mb2N, c("loo", "waic"))

# Spatio-temporal variogram -----------------------------------------------

mb2M_full_vario<-zoop_vario(model=mb2M_full, data=BL, yvar="CPUE", cores=5)
mb2M_full_vario_plot<-zoop_vario_plot(mb2M_full_vario)
ggsave(mb2M_full_vario_plot, filename="figures/Bosmina_variogram.png", device="png", width=8, height=5, units="in")


# Prediction plots --------------------------------------------------------

Data_effort<-BL%>%
  group_by(Month, Year)%>%
  summarise(N=n(), .groups="drop")%>%
  mutate(Month=as.integer(Month))

BL_preds<-zoop_predict(mb2M_full, BL, confidence=99)%>%
  left_join(Data_effort, by=c("Month", "Year"))%>%
  mutate(across(c(Pred, lowerCI, upperCI), ~if_else(is.na(N), NA_real_, .x)))

BL_salinity<-zoop_plot(BL_preds, "salinity")
BL_year<-zoop_plot(BL_preds, "year")
BL_season<-zoop_plot(BL_preds, "season")

ggsave(BL_season, file="figures/Bosmina_season.png", device="png", units = "in", width=8, height=6)

ggsave(BL_year, file="figures/Bosmina_year.png", device="png", units = "in", width=8, height=6)

ggsave(BL_salinity, file="figures/Bosmina_salinity.png", device="png", units = "in", width=8, height=6)

# Plot station intercepts -------------------------------------------------

p_intercepts<-zoop_stations(mb2M_full, select(Stations_clust, Clust, Latitude, Longitude))

ggsave(p_intercepts, file="figures/Bosmina_intercepts.png", device="png", units = "in", width=6, height=4)
