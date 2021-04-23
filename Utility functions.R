pp <- function(model, resp=NULL){
  prop_zero <- function(x) mean(x == 0)
  
  p<-list()
  
  p$zero<-pp_check(model, type="stat", stat=prop_zero, resp=resp)
  
  p$dist<-pp_check(model, resp=resp)+scale_x_log10()
  
  p$scatter<-pp_check(model, type="scatter_avg", resp=resp)+scale_y_log10()+scale_x_log10()
  
  return(p)
}

zoop_predict<-function(model, data){
  jdays<-expand_grid(Year=2001, Month=1:12, Day=seq(1, 26, by=5))%>%
    mutate(Julian_day=yday(ymd(paste(Year, Month, Day, sep="-"))))%>%
    filter(Julian_day>=min(data$Julian_day) & Julian_day<=max(data$Julian_day))
  
  newdata<-expand_grid(Salinity=quantile(data$SalSurf, probs=seq(0.05, 0.95, by=0.05)), 
                       Julian_day=jdays$Julian_day,
                       Year=unique(data$Year))%>%
    mutate(Year_s=(Year-mean(data$Year))/sd(data$Year),
           SalSurf_l_s=(log(Salinity)-mean(data$SalSurf_l))/sd(data$SalSurf_l),
           Julian_day_s=(Julian_day-mean(data$Julian_day))/sd(data$Julian_day))%>%
    left_join(jdays%>%
                dplyr::select(Julian_day, Month, Day),
              by="Julian_day")
  
  pred<-fitted(model, newdata=newdata, re_formula=NA, scale="response")
  
  newdata_pred<-newdata%>%
    mutate(Pred=pred[,"Estimate"],
           l95=pred[,"Q2.5"],
           u95=pred[,"Q97.5"])%>%
    mutate(Month2=month(Month, label=T))
  
  return(newdata_pred)
}

zoop_plot<-function(data, type){
  
  require(ggplot2)
  require(dplyr)
  require(stringr)
  
  if(!type%in%c("season", "year", "salinity")){
    stop('Valid types are "season", "year", or "salinity."')
  }
  
  
  if(type=="season"){
    data<-filter(data, Salinity%in%unique(data$Salinity)[seq(1,19, by=6)] & Year%in%seq(1975, 2020, by=5))
  } else{
    if(type=="year"){
      data<-filter(data, Salinity%in%unique(data$Salinity)[seq(1,19, by=6)] & Day==16)
    }else{
      data<-filter(data, Year%in%seq(1975, 2020, by=10) & Day==16)
    }
  }
  
  xvar<-case_when(type=="season" ~ "Julian_day", 
                  type=="year" ~ "Year", 
                  type=="salinity" ~ "Salinity")
  
  fillvar<-case_when(type=="season" ~ "Salinity", 
                     type=="year" ~ "Salinity", 
                     type=="salinity" ~ "Year")
  
  facetvar<-case_when(type=="season" ~ "Year", 
                      type=="year" ~ "Month2", 
                      type=="salinity" ~ "Month2")
  
  xlabel<-str_replace(xvar, "_", " ")
  
  if(type=="season"){
    scales<-list(scale_color_viridis_c(aesthetics = c("color", "fill"), trans="log", 
                                       breaks=unique(data$Salinity), 
                                       labels=round(unique(data$Salinity), 3),
                                       limits=range(data$Salinity)),
                 scale_x_continuous(breaks=c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335,
                                             15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349), 
                                    labels=c(rep("", 12), as.character(month(1:12, label = T))), limits=c(0,366),
                                    expand=expansion(0,0)),
                 theme(axis.ticks.x = element_line(color = c(rep("black", 12), rep(NA, 12))), axis.text.x=element_text(angle=45, hjust=1),
                       panel.grid.minor=element_blank(), panel.grid.major.x=element_line(color = c(rep("grey92", 12), rep(NA, 12)))))
  }else{
    if(type=="year"){
      scales<-list(scale_color_viridis_c(aesthetics = c("color", "fill"), trans="log", 
                                         breaks=unique(data$Salinity), 
                                         labels=round(unique(data$Salinity), 3),
                                         limits=range(data$Salinity)),
                   theme(axis.text.x=element_text(angle=45, hjust=1)))
      
    }else{
      scales<-list(scale_color_viridis_c(aesthetics = c("color", "fill")),
                   scale_x_continuous(trans="log", breaks=round(exp(seq(log(min(data$Salinity)), log(max(data$Salinity)), length.out=5)), 3),  minor_breaks = NULL),
                   theme(axis.text.x=element_text(angle=45, hjust=1)))
    }
    
  }
  
  p<-ggplot(data, aes(x=.data[[xvar]], y=Pred, ymin=l95, ymax=u95, fill=.data[[fillvar]], group=.data[[fillvar]]))+
    geom_ribbon(alpha=0.4)+
    geom_line(aes(color=.data[[fillvar]]))+
    facet_wrap(~.data[[facetvar]], scales = "free_y")+
    scale_y_continuous(expand=c(0,0))+
    ylab("CPUE")+
    xlab(xlabel)+
    theme_bw()+
    scales
  
  
  return(p)
}

zoop_vario<-function(model, data, cores=4){
  resids<-residuals(model, type="pearson")
  
  Data_vario<-data%>%
    mutate(Resid=resids[,"Estimate"])
  
  Data_coords<-Data_vario%>%
    st_as_sf(coords=c("Longitude", "Latitude"), crs=4326)%>%
    st_transform(crs=26910)%>%
    st_coordinates()%>%
    as_tibble()%>%
    mutate(across(c(X,Y), ~(.x-mean(.x))/1000))
  
  Data_vario<-bind_cols(Data_vario%>%
                          select(Date, Resid), Data_coords)
  sp<-SpatialPoints(coords=data.frame(X=Data_vario$X, Y=Data_vario$Y))
  sp2<-STIDF(sp, time=Data_vario$Date, 
             data=data.frame(Residuals=Data_vario$Resid))
  mb2M_vario<-variogramST(Residuals~1, data=sp2, tunit="weeks", cores=cores, tlags=seq(0,30, by=2))
  
  return(mb2M_vario)
}

zoop_vario_plot<-function(vario){
  require(patchwork)
  p_time<-ggplot(vario, aes(x=timelag, y=gamma, color=spacelag, group=spacelag))+
    geom_line()+
    geom_point()+
    scale_color_viridis_c(name="Distance (km)")+
    scale_x_continuous()+
    xlab("Time difference (weeks)")+
    theme_bw()+
    theme(legend.justification = "left")
  
  p_space<-ggplot(vario, aes(x=spacelag, y=gamma, color=timelag, group=timelag))+
    geom_line()+
    geom_point()+
    scale_color_viridis_c(name="Time difference\n(weeks)")+
    xlab("Distance (km)")+
    theme_bw()+
    theme(legend.justification = "left")
  
  p_variogram<-p_time/p_space+plot_annotation(tag_levels="A")
  return(p_variogram)
}

zoop_stations<-function(model, stations){
  
  station_intercepts<-posterior_summary(model, pars="r_Clust")%>%
    as_tibble(rownames="Clust")%>%
    mutate(Clust=as.numeric(str_extract(Clust, "(\\d)+")))%>%
    left_join(stations,
              by="Clust")
  
  p_intercepts<-ggplot()+
    geom_sf(data=deltamapr::WW_Delta%>%st_transform(crs=4326), fill="gray90", color="gray90")+
    geom_point(data=station_intercepts, aes(x=Longitude, y=Latitude, fill=Estimate), color="black", shape=21, alpha=0.7)+
    scale_fill_continuous_diverging(guide=guide_colorbar(barheight=10), palette="Blue-Red 3", name="Varying\nintercept\nestimate")+
    theme_bw()
  
  return(p_intercepts)
}
