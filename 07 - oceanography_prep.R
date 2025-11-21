rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #only works in rstudio / posit
source(file.path(getwd(),'_SRC','_startup.R'))

#### ----- Oceanography prep ---- #####
require(signal)
require(oce)
require(grid)
require(ggplot2)

GFX_OCE<-file.path(GFX_DIR,'OCE')
GFX_CTD<-file.path(GFX_OCE,'CTD')
dir.create(GFX_CTD,recursive=T,showWarnings = F)

out<-read.csv(MSM115_CTD_DATA,sep='\t',dec='.')
LL<-sf::st_as_sf(out,coords=c('lon','lat'),crs=WGS84)
LL<-sf::st_transform(LL,IBCSO)
xy<-as.data.frame(sf::st_coordinates(LL))
names(out)[which(names(out) %in% c('lon','lat'))]<-c('x','y')
out$x<-xy$X
out$y<-xy$Y

#### create unified and equispaced data
resolution<-.5
depth_pp <- seq(-200, 0, resolution)

# critical frequency for the filter (in ratios of Nyquist)
W <- resolution / 2
f1 <- signal::butter(1, W)
f2 <- signal::butter(2, W)

depth_col<-which(names(out) == 'depth_m')
measurement_cols<-which(names(out) %in% c("temp1_C", "temp2_C", "salinity", "oxygen_ml_l", "oxygen_saturation_perc", "fluorescence_mg_m3", "turbidity_NTU")) #excluding columns 15: sound_velocity_m_s and 16: bottle_sequence

unified_detrended<-data.frame()
for (id in unique(out$station_id)){
  sta<-out[out$station_id==id,]
  for (cast_type in unique(sta$cast)){
    print(id)
    print(cast_type)
    st<-subset(sta, cast==cast_type)
    new_data<-cbind(st[rep(1,length(depth_pp)),c(1:4)], data.frame(depth_m_unified=depth_pp))
    
    for (idx in measurement_cols){
      varname<-names(out)[idx]
      y<-st[,idx]
      
      equi_name<-paste0(varname,'_equispaced')
      detrend_name<-paste0(varname,'_detrended')
      f1_name<-paste0(varname,'_filter1')
      f2_name<-paste0(varname,'_filter2')
      
      new_data[[equi_name]]<-new_data[[detrend_name]]<-new_data[[f1_name]]<-new_data[[f2_name]]<-NA
      
      new_data[[equi_name]] <- Y0<- approx(st[,depth_col], y, new_data$depth_m_unified)$y
      Yd <- oce::detrend(new_data$depth_m_unified, Y0)
      non_na_idx<-which(!is.na(Yd$Y))
      Yd$Y<-Yd$Y[non_na_idx]
      new_data[[detrend_name]][non_na_idx] <- Yd$Y+ Yd$a + Yd$b * new_data$depth_m_unified[non_na_idx]
      new_data[[f1_name]][non_na_idx] <- signal::filtfilt(f1, Yd$Y) + Yd$a + Yd$b * new_data$depth_m_unified[non_na_idx]
      new_data[[f2_name]][non_na_idx] <- signal::filtfilt(f2, Yd$Y)+ Yd$a + Yd$b * new_data$depth_m_unified[non_na_idx]
    }
    unified_detrended<-rbind(unified_detrended,new_data)
  }
}
row.names(unified_detrended)<-NULL

filtered_data<-list(filter_resNULLfiltered_data<-list(filter_res=resolution, f1=f1, f2=f2, W=W, data=unified_detrended))

stations<-sf::st_read(MSM115_SPATIAL_DATA, layer='MSM115_stations') # Just the coordinates and id of each station
stations<-subset(stations,!is.na(station_id))

output_data<-list(stations = stations, station_data = out, filtered_data=filtered_data)

MSM115_OCE_DATA<-list('combined'=output_data,'unified_detrended'=unified_detrended)
save(MSM115_OCE_DATA, file = OCE_DATA_FILE, compress='gzip')

#### plotting ####
## added multiple x - axis
start<-which(names(unified_detrended)=='depth_m_unified')+1
vars<-names(unified_detrended[start:ncol(unified_detrended)])

df<-reshape2::melt(unified_detrended,
                   id.vars=c('station_id','depth_m_unified', 'cast'),
                   measure.vars=vars)
df$group[grepl('temp[1|2]_.*',df$variable)]<-'Temp [°C]'
df$group[grepl('salinity_.*',df$variable)]<-'Salinity'
df$group[grepl('oxygen_ml_.*',df$variable)]<-'Oxygen [ml/L]'
df$group[grepl('oxygen_saturation_.*',df$variable)]<-'Oxygen saturation [%]'
df$group[grepl('fluorescence_*',df$variable)]<-'Fluorescence [mg/m³]'
df$group[grepl('turbidity_*',df$variable)]<-'Turbidity [NTU]'
df$group<-factor(df$group, levels=c('Temp [°C]','Salinity','Oxygen [ml/L]','Oxygen saturation [%]','Fluorescence [mg/m³]',
                                    'Turbidity [NTU]'))

df$var_type[grepl('.*_filter2',df$variable)]<-'Filter 2'
df$var_type[grepl('.*_filter1',df$variable)]<-'Filter 1'
df$var_type[grepl('.*_detrended',df$variable)]<-'Detrended'
df$var_type[grepl('.*_equispaced',df$variable)]<-'Equispaced raw'
df$var_type<-factor(df$var_type, levels=c('Equispaced raw','Detrended','Filter 1','Filter 2'))

for (st in unique(df$station_id)){
  data<-df[df$station_id==st,]
  data$station_id<-paste0('MSM115 Station ',st)
  sub_data<-subset(data, var_type=='Filter 1')
  for (daynight in unique(sub_data$cast)){
    scaled_data<-subset(sub_data, var_type=='Filter 1' & cast == daynight)
    for (v in unique(scaled_data$variable)){
      for (cast_type in unique(scaled_data$cast)){
        idx<-which(scaled_data$variable == v & scaled_data$cast==cast_type)
        values<-scaled_data$value[idx]
        scaled_data$value_scaled[idx] <- range01(values)
        scaled_data$colour[idx]<-CTD_COLORMAP[[v]]
      }
    }
    
    scaled_data<-scaled_data[-which(scaled_data$variable=='temp2_C_filter1'),]
    
    ### one plot for all augmentation
    p<-ggplot(data=data, aes(y=value, x=depth_m_unified, col=var_type, group=group))
    p<-p+geom_line()+coord_flip()
    p<-p+facet_grid(var_type~station_id+cast+group, scales='free_x')+MAINTHEME
    
    png(file.path(GFX_CTD,paste0('MSM115_station_', st, '_detrending.png')),1000,1000,res=150)
    print(p)
    graphics.off()
    
    ### one plot per station and all vars
    # using grid to stack plots and take axis
    # need to use fixed colors so we can use them everywhere
    table(scaled_data$group, scaled_data$colour)
    
    p<-ggplot(data=scaled_data, aes(y=value_scaled, x=depth_m_unified, colour=group))
    p<-p+geom_line(lwd=1)+coord_flip()
    p<-p+scale_color_manual(values=unique(scaled_data$colour))
    p<-p+facet_wrap(station_id+var_type~cast, scales='free_x')
    p<-p+labs(col=NULL)+ guides(colour = guide_legend(nrow = 1))
    p<-p+theme_linedraw()+ theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
    p<-p+theme(
      legend.position = c(0.5, 0),
      legend.justification = c("center", "bottom"),
      legend.margin = ggplot2::margin(2, 2, 2, 2),
      legend.key = element_rect(fill =  "transparent")
    )
    
    plot_obj_dark<-ggplotGrob(p+theme_ctd())
    
    for (v in unique(scaled_data$variable)){
      pretty_name<-unique(scaled_data$group[scaled_data$variable==v])
      current_color=CTD_COLORMAP[[v]]
      values<-round(as.numeric(quantile(scaled_data$value[scaled_data$variable==v], na.rm=T)),4)
      f<-ggplot(data=data.frame(x=c(0,1), y=c(0,1)), aes(x=x,y=y))
      f<-f+scale_x_continuous(NULL, breaks = seq(0,1,.25), labels = values, limits=c(0,1))
      f<-f+ coord_fixed(0/5)
      f<-f+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),axis.title.y=element_blank(),
                 legend.position="none",
                 panel.background=element_blank(),
                 plot.background=element_blank())
      f<-f+theme(axis.text.x = element_text( color=current_color, face=3),
                 axis.ticks.x = element_line(size = 1.5, color=current_color),
                 axis.line = element_line(size = 1.5, colour = current_color, linetype=1)
      )
      plot_obj_dark<-rbind(plot_obj_dark, ggplotGrob(f+theme_axis(current_color)))
    }
    
    png(file.path(GFX_CTD,paste0('MSM115_station_', st, '_',daynight,'_overview.png')),2500,3000,res=300)
    grid::grid.newpage()
    grid::grid.draw(plot_obj_dark)
    graphics.off()
  }
}

