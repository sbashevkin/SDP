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

pp <- function(model){
  prop_zero <- function(x) mean(x == 0)
  
  p<-list()
  
  p$zero<-pp_check(model, type="stat", stat=prop_zero)
  
  p$dist<-pp_check(model)+scale_x_log10()
  
  return(p)
}

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
  mutate_at(vars(Date_num, SalSurf, SalSurf_l, Temperature, Latitude, Longitude, Julian_day, Year), list(s=~(.-mean(., na.rm=T))/sd(., na.rm=T)))

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

# Spatio-temporal variogram -----------------------------------------------

resid_mb2M<-residuals(mb2M_full, type="pearson")
resid_mb2MB<-residuals(mb2M_full, type="ordinary")

Data_vario<-BL%>%
  mutate(Resid=resid_mb2M[,"Estimate"],
         ResidB=resid_mb2MB[,"Estimate"])

Data_coords<-Data_vario%>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326)%>%
  st_transform(crs=26910)%>%
  st_coordinates()%>%
  as_tibble()%>%
  mutate(across(c(X,Y), ~(.x-mean(.x))/1000))

Data_vario<-bind_cols(Data_vario%>%
                        select(Date, Resid, ResidB), Data_coords)
sp<-SpatialPoints(coords=data.frame(X=Data_vario$X, Y=Data_vario$Y))
sp2<-STIDF(sp, time=Data_vario$Date, 
           data=data.frame(Residuals=Data_vario$Resid, ResidualsB=Data_vario$ResidB))
mb2M_vario<-variogramST(Residuals~1, data=sp2, tunit="weeks", cores=5, tlags=seq(0,30, by=2))
mb2M_varioB<-variogramST(ResidualsB~1, data=sp2, tunit="weeks", cores=5, tlags=seq(0,30, by=2))

p_time<-ggplot(mb2M_vario, aes(x=timelag, y=gamma, color=spacelag, group=spacelag))+
  geom_line()+
  geom_point()+
  scale_color_viridis_c(name="Distance (km)")+
  scale_x_continuous()+
  xlab("Time difference (weeks)")+
  theme_bw()+
  theme(legend.justification = "left")

p_space<-ggplot(mb2M_vario, aes(x=spacelag, y=gamma, color=timelag, group=timelag))+
  geom_line()+
  geom_point()+
  scale_color_viridis_c(name="Time difference\n(weeks)")+
  xlab("Distance (km)")+
  theme_bw()+
  theme(legend.justification = "left")

p_variogram<-p_time/p_space+plot_annotation(tag_levels="A")

#ggsave(p_variogram, filename="C:/Users/sbashevkin/OneDrive - deltacouncil/Discrete water quality analysis/Manuscripts/Climate change/Figures/variogram.png",
#       device="png", width=8, height=5, units="in")


# Prediction plots --------------------------------------------------------

jdays<-expand_grid(Year=2001, Month=1:12, Day=seq(1, 26, by=5))%>%
  mutate(Julian_day=yday(ymd(paste(Year, Month, Day, sep="-"))))%>%
  filter(Julian_day>=min(BL$Julian_day) & Julian_day<=max(BL$Julian_day))

newdata<-expand_grid(Salinity=quantile(BL$SalSurf, probs=seq(0.05, 0.95, by=0.05)), 
                     Julian_day=jdays$Julian_day,
                     Year=unique(BL$Year))%>%
  mutate(Year_s=(Year-mean(BL$Year))/sd(BL$Year),
         SalSurf_l_s=(log(Salinity)-mean(BL$SalSurf_l))/sd(BL$SalSurf_l),
         Julian_day_s=(Julian_day-mean(BL$Julian_day))/sd(BL$Julian_day))%>%
  left_join(jdays%>%
            select(Julian_day, Month, Day),
            by="Julian_day")

pred<-fitted(mb2M, newdata=newdata, re_formula=NA, scale="response")

newdata_pred<-newdata%>%
  mutate(Pred=pred[,"Estimate"],
         l95=pred[,"Q2.5"],
         u95=pred[,"Q97.5"])

p_season<-ggplot(filter(newdata_pred, Salinity%in%unique(newdata$Salinity)[seq(1,19, by=6)] & Year%in%seq(1975, 2020, by=5)),
       aes(x=Julian_day, y=Pred, ymin=l95, ymax=u95, fill=Salinity, group=Salinity))+
  geom_ribbon(alpha=0.4)+
  geom_line(aes(color=Salinity))+
  facet_wrap(~Year, scales = "free_y")+
  scale_color_viridis_c(aesthetics = c("color", "fill"), trans="log", 
                        breaks=round(unique(newdata$Salinity)[seq(1,19, by=6)], 3), 
                        limits=round(unique(newdata$Salinity)[c(1,19)], 3))+
  ylab("CPUE")+
  xlab("Julian day")+
  theme_bw()

p_year<-ggplot(filter(newdata_pred, Salinity%in%unique(newdata$Salinity)[seq(1,19, by=6)] & Day==16),
       aes(x=Year, y=Pred, ymin=l95, ymax=u95, fill=Salinity, group=Salinity))+
  geom_ribbon(alpha=0.4)+
  geom_line(aes(color=Salinity))+
  facet_wrap(~month(Month, label=T), scales = "free_y")+
  scale_color_viridis_c(aesthetics = c("color", "fill"), trans="log", 
                        breaks=round(unique(newdata$Salinity)[seq(1,19, by=6)], 3), 
                        limits=round(unique(newdata$Salinity)[c(1,19)], 3))+
  ylab("CPUE")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45, hjust=1))

p_salinity<-ggplot(filter(newdata_pred, Year%in%seq(1975, 2020, by=10) & Day==16),
       aes(x=Salinity, y=Pred, ymin=l95, ymax=u95, fill=Year, group=Year))+
  geom_ribbon(alpha=0.4)+
  geom_line(aes(color=Year))+
  facet_wrap(~month(Month, label=T), scales = "free_y")+
  scale_x_continuous(trans="log", breaks=round(exp(seq(log(min(newdata$Salinity)), log(max(newdata$Salinity)), length.out=5)), 3),  minor_breaks = NULL)+
  scale_color_viridis_c(aesthetics = c("color", "fill"))+
  ylab("CPUE")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45, hjust=1))

ggsave(p_season, file="Figures/Bosmina_season.png", device="png", units = "in", width=8, height=6)

ggsave(p_year, file="Figures/Bosmina_year.png", device="png", units = "in", width=8, height=6)

ggsave(p_salinity, file="Figures/Bosmina_salinity.png", device="png", units = "in", width=8, height=6)
