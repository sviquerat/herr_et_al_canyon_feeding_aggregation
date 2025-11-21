rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #only works in rstudio / posit
source(file.path(getwd(),'_SRC','_startup.R'))

#### ----- Impact of currents ---- #####
require(dplyr)
require(ncdf4)
require(terra)

GFX_CURRENTS<-file.path(GFX_DIR,'CURRENTS')
dir.create(GFX_CURRENTS,showWarnings = F)

strata<-sf::st_read(MSM115_STRATA_FILE,layer='simplified_strata')
strata<-sf::st_minimum_rotated_rectangle(strata)
strata<-sf::st_transform(strata,WGS84)

shiptrack<-sf::st_read(MSM115_EXPEDITION_DATA, layer='MSM115_Ship_Track_Points')
shiptrack<-sf::st_transform(shiptrack, WGS84)

load(MSM115_CLUSTER_DATA)
posteriori_cluster<-MSM115_CLUSTER$posteriori

### need to add station dates back
abu <- openxlsx::read.xlsx(MSM115_KRILL_ABU_FILE,sheet=1)
abu$date2<-as.Date(ceiling(abu$START.DATE), origin = "1899-12-30")
idx<-grep('night',abu$FINWAP.STATION)
abu<-abu[-idx,]
stations<-sqldf::sqldf('select station,date2 from abu group by station,date2')

idx<-which(posteriori_cluster$recordType=='station')
st<-posteriori_cluster[idx,]
st<-sqldf::sqldf('select * from st as a left join stations as b on a.id = b.STATION')
posteriori_cluster$date[idx]<-st$date2

aggr<-posteriori_cluster[posteriori_cluster$recordType=='aggr',]
aggr<-aggr[,which(names(aggr) %in% c('date','x','y'))]
aggr$source<-'a priori classification'
#names(aggr)[which(names(aggr) == 'recordType')]<-'classification'

pc<-posteriori_cluster[posteriori_cluster$cluster_name=='aggregations',]
pc<-pc[,which(names(pc) %in% c('date','x','y'))]
pc$source<-'a posteriori classification'

data<-rbind(aggr,pc)
data$jday<-as.numeric(format(data$date,'%j'))
data$year<-as.numeric(format(data$date,'%Y'))

##### spatial manipulation & block creation
PP<-sf::st_as_sf(data,coords=c('x','y'),crs=IBCSO)
PP<-sf::st_transform(PP,WGS84)
PP_buffer<-sf::st_buffer(PP,25000)
PP_buffer<-sf::st_join(PP_buffer,sf::st_buffer(strata,10000))

derived_strata <- PP_buffer %>% 
  group_by(ShortName) %>% 
  summarize(geometry = sf::st_as_sfc(sf::st_bbox(sf::st_union(geometry))))

derived_strata$ShortName[is.na(derived_strata$ShortName)]<-'beyond strata'

p<-ggplot()+geom_sf(data=strata)
p<-p+geom_sf(data=derived_strata,mapping=aes(fill=ShortName))
p<-p+geom_sf(data=PP,mapping=aes(pch=source,col=source),cex=2)
p<-p+scale_colour_discrete(name = "Data source") + scale_shape_discrete(name = "Data source")
p+labs(fill='Reference block')

# EI is a bit too narrow
poly<-derived_strata[derived_strata$ShortName=='EI',]
coords <- sf::st_coordinates(poly$geometry)
coords[3:4,2]<-coords[3:4,2]+.3
updated_geometry<-sf::st_polygon(list(coords[,1:2]))
derived_strata$geometry[derived_strata$ShortName=='EI'] <- updated_geometry

shiptrack<-sf::st_join(shiptrack,derived_strata)
shiptrack<-subset(shiptrack, !is.na(ShortName))

p<-ggplot()+geom_sf(data=strata,mapping=aes(fill=ShortName),alpha=.4)
p<-p+geom_sf(data=derived_strata,mapping=aes(fill=ShortName),col='black')
p<-p+geom_sf(data=PP_buffer,mapping=aes(pch=source,col=source),cex=2)
p<-p+geom_sf(data=PP,mapping=aes(pch=source,col=source),cex=2)
p<-p+scale_colour_discrete(name = "Data source") + scale_shape_discrete(name = "Data source")
p<-p+labs(fill='Reference block')

filename<-file.path(PUB_PACKAGE_DIR,'supp_Figure_8.png')
png(filename,6400,3200,res=300)
print(p+MAINTHEME+theme_large_font(14)+theme(legend.position='bottom'))
graphics.off()

sf::st_write(derived_strata, dsn=CURRENTS_STRATA_GPKG,layer='derived_strata',append=F, delete_dsn=T)
sf::st_write(strata, dsn=CURRENTS_STRATA_GPKG,layer='original_strata',append=T)

PP<-sf::st_join(PP,derived_strata)

# currents for 2023
# NOTE: the first time through will create all neccessary data at CURRENTS_DATA_DIR. On subsequent runs, the old data will be reused
dayrange<-10 # consider +- these dates
depth_layer=1 # depth_layer 1 is at 0m, depth_layer 2 at 15m

out<-data.frame()
for (block_id in seq_along(unique(derived_strata$ShortName))){
  block<-unique(derived_strata$ShortName)[block_id]
  stratum<-subset(derived_strata,ShortName==block)
  obs<-subset(PP,ShortName==block)
  
  days<-unique(obs$date)
  start_day<-min(days)-dayrange
  end_day<-max(days)+dayrange
  
  for (target_day in seq(start_day,end_day,1)){
    day<-as.Date(target_day)
    is_aggregation_day=F
    source<-NA
    observer_present<-F
    if (day %in% obs$date){
      is_aggregation_day=T
      source<-paste0(unique(obs$source[obs$date==day]),collapse = ' and ')
    }
    
    if (day %in% shiptrack$date){
      observer_present<-T
    }
    stack<-create_currents_raster(COPERNICUS_CURRENTS_NC,date=day,depth_layer=depth_layer,output_folder=CURRENTS_DATA_DIR, force_return = T)
    names(stack)<-c('current_direction','current_velocity')
    
    points <- terra::as.points(stack)
    
    # Convert to sf object for easier handling
    pp <- sf::st_as_sf(points)
    
    p_strata<-sf::st_crop(pp,stratum)
    avg_dir<-yamartino_avg(p_strata$current_direction)
    sd_dir<-yamartino_std(p_strata$current_direction)
    se_dir<-sd_dir/sqrt(nrow(p_strata))
    
    avg_velo<-mean(p_strata$current_velocity)
    sd_velo<-sd(p_strata$current_velocity)
    se_velo<-sd_velo/sqrt(nrow(p_strata))
    new_line<-data.frame(block=block,day=day, aggregation=is_aggregation_day,source=source, observer_present=observer_present,
                         avg_surface_dir=avg_dir, sd_surface_dir=sd_dir, se_surface_dir=se_dir,
                         avg_velocity = avg_velo, sd_velocity=sd_velo, se_velocity=se_velo)
    out<-rbind(out, new_line)
  }
}

out$source[is.na(out$source)]<-'no observation'

vars<-c('avg_surface_dir','sd_surface_dir','avg_velocity', 'sd_velocity')
pretty_covars<-c('surface current direction [°]',
                 'standard deviation of surface current direction',
                 'surface current velocity [m/s]',
                 'standard deviation of surface velocity')
pretty_covars<-c(expression(paste(theta[average], ' [°]')),
                 expression(paste(sigma,'(',theta,')')),
                 expression(paste(nu[average], ' [m/s]')),
                 expression(paste(sigma,'(',nu,')'))
)

for (block in unique(out$block)){
  data<-out[out$block==block,]
  aggregation_days<-subset(data,aggregation==T)
  observer_days<-subset(data,observer_present==T)
  visit<-1
  observer_days$visit<-visit
  for (r in 2:nrow(observer_days)){
    if (observer_days$day[r]-1==observer_days$day[r-1]){
      visit<-observer_days$visit[r-1]
    }else{
      visit<-visit+1
    }
    observer_days$visit[r]<-visit
  }
  observer_days<-sqldf::sqldf('select visit,min(day) as start, max(day) as end from observer_days group by visit')
  observer_days$start<-as.Date(observer_days$start)
  observer_days$end<-as.Date(observer_days$end)
  observer_days$pretty_visit<-paste0(observer_days$visit, 'st: ', observer_days$start, ' - ', observer_days$end)
  observer_days$pretty_visit<-as.factor(observer_days$pretty_visit)
  
  p<-ggplot(data=data,aes(x=day))
  p<-p+geom_rect(data = observer_days, aes(xmin = start-.5, xmax = end+.5,fill=pretty_visit), alpha = 0.3, col="grey", inherit.aes = FALSE,ymin=0,ymax=.2)
  p<-p+scale_fill_manual(values=obs_presence_colors)
  p<-p+geom_vline(data=aggregation_days, aes(xintercept=day),col='black',alpha=.2,lty=2,lwd=1)
  p<-p+geom_point(data=aggregation_days,mapping=aes(col=source),y=.3,pch='\U02605',cex=8)
  p<-p+scale_colour_manual(values=cluster_source_colors)
  p<-p+annotate(geom = "text", x = data$day, y = .5, label = '\U02191', color = "red",angle = 360-data$avg_surface_dir,cex=8,hjust=.5,vjust=0)
  p<-p+labs(colour='Source',fill='Observer presence')+ylab('')+xlab('')
  p<-p+theme_minimal()+theme_large_font()
  p<-p+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.text.y=element_blank())
  p<-p+theme(panel.grid.minor.x = element_blank(),panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank())
  p<-p+scale_x_date(breaks=unique(data$day))
  p<-p+ylim(c(0,1))
  p
  
  cols<-which(names(data) %in% c('day', 'is_aggregation','source'))
  plotlist=list()
  
  for (i in seq_along(vars)){ #this outside loop - it's always the same
    var<-vars[i]
    
    df<-data[,cols]
    df$value<-data[[var]]
    
    se<-F
    if (grepl('avg_.*',var)){
      se<-T
      df$se_value<-data[[paste0('se_',strsplit(var,'avg_')[[1]][2])]]  
      poly<-data.frame(x=c(df$day,rev(df$day)),y=c(df$value+qnorm(0.975)*df$se_value,rev(df$value-qnorm(0.975)*df$se_value)))
    }
    
    moving_average<-ma(df$value,3)
    ma_data<-data.frame(x=df$day,y=moving_average)
    
    g<-ggplot(data=df, aes(x=day,y=value))
    g<-g+geom_rect(data = observer_days, aes(xmin = start-.5, xmax = end+.5,fill=pretty_visit), alpha = 0.3, col="grey", inherit.aes = FALSE,ymin=min(df$value),ymax=max(df$value))
    g<-g+scale_fill_manual(values=obs_presence_colors)
    
    g<-g+labs(y=pretty_covars[i]) + theme_minimal()
    
    if (se){
      g<-g+geom_polygon(data=poly,inherit.aes=F,aes(x=x,y=y),fill='grey',col='darkgrey',alpha=.8)
    }
    
    g<-g+geom_vline(data=aggregation_days, aes(xintercept=day),col='black',alpha=.2,lty=2,lwd=1)
    g<-g+scale_colour_manual(values=obs_presence_colors)
    
    g<-g+geom_line(col='black',lwd=1.5)
    g<-g+geom_line(data=ma_data[complete.cases(ma_data),], mapping=aes(x=x,y=y),col='red',lwd=1.25,alpha=.7)
    g<-g+scale_x_date(breaks=unique(data$day),labels=format(unique(data$day),'%m-%d'))
    g<-g+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    if (i != length(vars)){
      g<-g+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
    }else{
    }

    if (i %%2 ==0){
      print('')
    }
    plotlist[[i]]<-g+theme_large_font()+theme(plot.margin=unit(c(0,0,0,0),'cm')) 
  }
  plotlist[[length(plotlist)]]<-plotlist[[length(plotlist)]]+labs(x='')
  
  fpl<-list(p)
  fpl<-append(fpl,plotlist)
  leg<-ggpubr::get_legend(p + theme(legend.key.size = unit(1, 'cm'),legend.text = element_text(size=10)))
  
  f<-ggpubr::ggarrange(plotlist=fpl, ncol=1,align='v',heights=c(4,rep(3,length(vars))),legend.grob=leg,legend='right',labels='auto')

  filename<-file.path(GFX_CURRENTS,paste0('MSM115_currents_',block,'.png'))
  png(filename,6400,4200,res=300)
  print(f)
  graphics.off()
  
  if (block=='EI'){
    filename<-file.path(PUB_PACKAGE_DIR,'Figure_16.png')
    png(filename,6400,3600,res=300)
    print(f)
    graphics.off()
    
    filename<-file.path(PUB_PACKAGE_DIR,'supp_Figure_9.png')
    png(filename,6400,3600,res=300)
    print(f)
    graphics.off()
  }
  if (block=='SOI'){
    filename<-file.path(PUB_PACKAGE_DIR,'supp_Figure_10.png')
    png(filename,6400,3600,res=300)
    print(f)
    graphics.off()
  }
  if (block=='SSI'){
    filename<-file.path(PUB_PACKAGE_DIR,'supp_Figure_11.png')
    png(filename,6400,3600,res=300)
    print(f)
    graphics.off()
  }
  if (block=='beyond strata'){
    filename<-file.path(PUB_PACKAGE_DIR,'supp_Figure_12.png')
    png(filename,6400,3600,res=300)
    print(f)
    graphics.off()
  }
}
