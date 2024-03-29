pp <- function(model, resp=NULL){
  prop_zero <- function(x) mean(x == 0)
  
  p<-list()
  
  p$zero<-pp_check(model, type="stat", stat=prop_zero, resp=resp)
  
  p$dist<-pp_check(model, resp=resp)+scale_x_log10()
  
  p$scatter<-pp_check(model, type="scatter_avg", resp=resp)+scale_y_log10()+scale_x_log10()
  
  return(p)
}

zoop_predict<-function(model, data, confidence=95){
  prob_lower<-(100-confidence)/200
  probs<-c(prob_lower, 1-prob_lower)
  
  quantiles<-paste0("Q", probs*100)
  
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
  
  pred<-fitted(model, newdata=newdata, re_formula=NA, scale="response", probs=probs)
  
  newdata_pred<-newdata%>%
    mutate(Pred=pred[,"Estimate"],
           lowerCI=pred[,quantiles[1]],
           upperCI=pred[,quantiles[2]])%>%
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
  
  if(min(data$Year)>1975){
    if(min(data$Year)>2000){
      plot_years1<-seq(2006, 2019, by=3)
      plot_years2<-seq(2006, 2019, by=3)
    }else{
      plot_years1<-seq(1975, 2020, by=5)
      plot_years2<-seq(1995, 2015, by=5)
    }
  }else{
    plot_years1<-seq(1975, 2020, by=5)
    plot_years2<-seq(1975, 2020, by=10)
  }
  
  
  if(type=="season"){
    data<-filter(data, Salinity%in%unique(data$Salinity)[seq(1,19, by=6)] & Year%in%plot_years1)
  } else{
    if(type=="year"){
      data<-filter(data, Salinity%in%unique(data$Salinity)[seq(1,19, by=6)] & Day==16)
    }else{
      data<-filter(data, Year%in%plot_years2 & Day==16)
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
  
  xlabel<-case_when(xvar=="Julian_day" ~ "Day of year", 
                    xvar=="Salinity" ~ "Salinity",
                    TRUE ~ xvar)
  
  data_orphans<-data%>%
    group_by(.data[[fillvar]], .data[[facetvar]])%>%
    mutate(Pred_lag=lag(Pred, order_by = .data[[xvar]]), Pred_lead=lead(Pred, order_by = .data[[xvar]]))%>%
    ungroup()%>%
    filter(is.na(Pred_lag) & is.na(Pred_lead) & !is.na(Pred))%>%
    select(-Pred_lag, -Pred_lead)
  
  if(type=="season"){
    scales<-list(scale_color_viridis_d(breaks=unique(factor(data$Salinity)),
                                       labels=round(unique(data$Salinity), 3),
                                       aesthetics = c("color", "fill"), 
                                       end=0.8, 
                                       name="Salinity"),
                 scale_linetype(breaks=unique(factor(data$Salinity)),
                                labels=round(unique(data$Salinity), 3),
                                name="Salinity"),
                 scale_x_continuous(breaks=c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335,
                                             15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349), 
                                    labels=c(rep("", 12), as.character(month(1:12, label = T))), limits=c(0,366),
                                    expand=expansion(0,0)),
                 theme(axis.ticks.x = element_line(color = c(rep("black", 12), rep(NA, 12))), axis.text.x=element_text(angle=45, hjust=1),
                       panel.grid.minor=element_blank(), panel.grid.major.x=element_line(color = c(rep("grey92", 12), rep(NA, 12)))))
  }else{
    if(type=="year"){
      scales<-list(scale_color_viridis_d(breaks=unique(factor(data$Salinity)),
                                         labels=round(unique(data$Salinity), 3),
                                         aesthetics = c("color", "fill"), 
                                         end=0.8, 
                                         name="Salinity"),
                   scale_linetype(breaks=unique(factor(data$Salinity)),
                                  labels=round(unique(data$Salinity), 3),
                                  name="Salinity"),
                   theme(axis.text.x=element_text(angle=45, hjust=1)))
      
    }else{
      scales<-list(scale_color_viridis_d(aesthetics = c("color", "fill"), end=0.8, name="Year"),
                   scale_linetype(name="Year"),
                   scale_x_continuous(trans="log", breaks=round(exp(seq(log(min(data$Salinity)), log(max(data$Salinity)), length.out=5)), 3),  minor_breaks = NULL),
                   theme(axis.text.x=element_text(angle=45, hjust=1)))
    }
    
  }
  
  p<-ggplot(data, aes(x=.data[[xvar]], y=Pred, ymin=lowerCI, ymax=upperCI, 
                      fill=factor(.data[[fillvar]]), group=.data[[fillvar]]))+
    geom_ribbon(alpha=0.4)+
    geom_line(aes(color=factor(.data[[fillvar]]), linetype=factor(.data[[fillvar]])))+
    geom_pointrange(data=data_orphans, aes(color=factor(.data[[fillvar]])), 
                    shape=21, position=position_dodge(width=2), size=0.2,
                    show.legend=F)+
    facet_wrap(~.data[[facetvar]], scales = "free_y")+
    scale_y_continuous(expand=c(0,0), limits = c(0, NA)) +
    ylab(bquote("Abundance (# "*m^-3*")"))+
    xlab(xlabel)+
    theme_bw()+
    scales
  
  
  return(p)
}

zoop_vario<-function(model, data, yvar, resid_type="standardized", cores=4){
  require(sp)
  require(gstat)
  require(spacetime)
  require(brms)
  require(dplyr)
  require(sf)
  
  if(resid_type=="standardized"){
    resids<-residuals(model, method="posterior_predict")/sd(data[[yvar]])
  }else{
    if(resid_type=="deviance"){
      resids<-residuals(model, method="posterior_predict")^2 # This seems to be the method used to calculate deviance residuals by mgcv
    }else
      stop("Only 'standardized' or 'deviance' resid_type values are accepted")
  }
  
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
              by="Clust")%>%
    st_as_sf(coords=c("Longitude", "Latitude"), crs=4326)%>%
    st_transform(crs=st_crs(deltamapr::WW_Delta))
  
  map<-deltamapr::WW_Delta%>%
    st_crop(station_intercepts%>%
              st_union()%>%
              st_buffer(dist = units::set_units(1, "km")))
  
  p_intercepts<-ggplot()+
    geom_sf(data=map, fill="gray90", color="gray90")+
    geom_sf(data=station_intercepts, aes(fill=Estimate), color="black", shape=21, alpha=0.7)+
    scale_fill_continuous_diverging(guide=guide_colorbar(barwidth=10, title.hjust=0.5, direction="horizontal", title.position = "top"), 
                                    palette="Blue-Red 3", name="Clustered station estimate")+
    coord_sf(expand=FALSE)+
    theme_bw()+
    theme(legend.position=c(0.3, 0.85), legend.background=element_rect(color="black"), legend.margin = margin(10, 10, 10, 10))
  
  return(p_intercepts)
}
