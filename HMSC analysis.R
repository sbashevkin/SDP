require(dplyr)
require(zooper)
require(lubridate)
require(hms)
require(tidyr)
require(Hmsc)

PF<-zoopComb%>%
  filter(Taxname=="Pseudodiaptomus forbesi")%>%
  left_join(zoopEnvComb%>%
              select(-Source),
            by="SampleID")%>%
  filter(Year>=1995)%>% # Introduced in 1987, give them 8 years to get established
  mutate(Date_num=as.numeric(Date))%>%
  mutate(SalSurf_l = log(SalSurf),
         Julian_day = yday(Date))%>%
  drop_na(Date, SalSurf, Temperature, Latitude, Longitude, CPUE)%>%
  mutate(Station=if_else(Station%in%unique(zooper::stationsEMPEZ$Station), paste(Source, Station, year(Date), month(Date), day(Date)), paste(Source, Station)))%>%
  mutate(Date2=if_else(is.na(Datetime), parse_date_time(paste0(year(Date), "-", month(Date), "-", day(Date), " 12:00"), "%Y-%m-%d %H:%M", tz="America/Los_Angeles"), Datetime))%>%
  group_by_at(vars(-CPUE, -Volume, -SampleID))%>%
  summarize(CPUE=mean(CPUE), Volume=mean(Volume), .groups="drop")%>%
  mutate(Count=round(CPUE*Volume),
         Year_fac=factor(Year),
         Month=month(Date, label = T),
         SampleID=1:n())%>%
  mutate_at(vars(Date_num, SalSurf, SalSurf_l, Temperature, Latitude, Longitude, Julian_day), list(s=~(.-mean(., na.rm=T))/sd(., na.rm=T)))

Y <- as.matrix(PF$CPUE)
XData <- data.frame(x = PF$SalSurf_l)
m <- Hmsc(Y = Y, XData = XData, XFormula = ~x, distr="lognormal poisson")

nChains<-2
thin<-5
samples<-1000
transient<-500*thin
verbose<-500*thin

m<-sampleMcmc(m, thin = thin, samples = samples, transient = transient,
               nChains = nChains, verbose = verbose)

mpost <- convertToCodaObject(m)
summary(mpost$Beta)

preds <- computePredictedValues(m, expected = FALSE)
evaluateModelFit(hM=m, predY=preds)

plot(mpost$Beta)
effectiveSize(mpost$Beta)
gelman.diag(mpost$Beta,multivariate=FALSE)$psrf


preds.mean <- apply(preds, FUN=mean, MARGIN=1)
nres<-scale(Y-preds.mean)
par(mfrow=c(1,2))
hist(nres, las = 1)
plot(preds.mean,nres, las = 1)
abline(a=0,b=0)

## Spatial model
Stations<-distinct(PF, Longitude, Latitude, Station)
xycoords <- as.matrix(select(Stations, Longitude, Latitude))
rownames(xycoords) <- Stations$Station
colnames(xycoords) <- c("x-coordinate","y-coordinate")

studyDesign <- data.frame(Station = PF$Station)
rL <- HmscRandomLevel(sData = xycoords, sMethod = 'NNGP', nNeighbours = 2)
rL = setPriors(rL,nfMin=1,nfMax=1)
m <- Hmsc(Y=Y, XData=XData, XFormula=~x,
         studyDesign=studyDesign, ranLevels=list("Station"=rL), distr="lognormal poisson")

m <- sampleMcmc(m, thin = thin, samples = samples, transient = transient,
               nChains = nChains, verbose = verbose)
