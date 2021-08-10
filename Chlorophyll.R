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

load("Bosmina Chlorophyll.Rdata")

iterations <- 5e3
warmup <- iterations/4

mchl<-brm(bf(CPUE ~ t2(Chl_s, Julian_day_s, k=c(5,13), bs=c("cr", "cc")) + (1|Clust), hu ~ s(SalSurf_l_s, bs="cr", k=5)),
          data=BL_chl, family=hurdle_lognormal(),
          prior=prior(normal(0,10), class="Intercept")+
            prior(normal(0,5), class="b")+
            prior(cauchy(0,5), class="sigma"),
          chains=1, cores=1, control=list(adapt_delta=0.99),
          iter = iterations, warmup = warmup,
          backend = "cmdstanr", threads = threading(4))

mchl2<-brm(bf(CPUE ~ t2(Chl_s, SalSurf_l_s, Julian_day_s, k=c(5,5,6), bs=c("cr", "cr", "cc")) + (1|Clust), hu ~ s(SalSurf_l_s, bs="cr", k=5)),
          data=BL_chl, family=hurdle_lognormal(),
          prior=prior(normal(0,10), class="Intercept")+
            prior(normal(0,5), class="b")+
            prior(cauchy(0,5), class="sigma"),
          chains=1, cores=1, control=list(adapt_delta=0.99, max_treedepth=15),
          iter = iterations, warmup = warmup,
          backend = "cmdstanr", threads = threading(4))
