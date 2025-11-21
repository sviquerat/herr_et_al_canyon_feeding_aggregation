#### Assortment of various functions

draw_components<-function(spp, multimodality_data, length_data,abbreviate_spp=T, opdf_colour='#e3adea',dip_test_size=10){
  spp_data<-subset(multimodality_data,species==spp)
  
  if (nrow(spp_data)==0){return(NULL)}
  
  spp_length<-subset(length_data,species==spp & station_id %in% spp_data$station_id)
  spp_length$stratum<-sqldf::sqldf('select * from spp_length as a left join station_strata_table as b on a.station_id=b.station_id')$stratum
  spp_length$station<-paste0('Station ',formatC(spp_length$station_id,width=2,flag='0'),spp_length$dip_sig)
  
  spp_lbl=spp
  if (abbreviate_spp){
    lbl<-strsplit(spp_lbl,' ')[[1]]
    lbl[1]<-paste0(substr(lbl[1],1,1),'.')
    spp_lbl<-paste0(lbl, collapse=' ')
  }
  pdfs<-data.frame()
  for (st in unique(spp_data$station_id)){
    st_data<-spp_data[spp_data$station_id==st,]
    for (cmp in st_data$component){
      mu<-st_data$mu[st_data$component==cmp]
      std<-st_data$stdev[st_data$component==cmp]
      max_length=mu+3*std
      x<-seq(0,max_length,.1)
      y<-dnorm(x,mu, sd = std)
      new_line<-data.frame(station=unique(st_data$station),stratum=unique(st_data$stratum),station_id=st,component=cmp,x=x,y=y)
      pdfs<-rbind(pdfs,new_line)
    }
  }
  pdfs$component<-as.factor(pdfs$component)
  spp_data$component<-as.factor(spp_data$component)
  spp_length$fill_col<-paste0('Pooled distribution for \n',spp_lbl)
  
  p<-ggplot(data=spp_length)
  p<-p+geom_histogram(mapping=aes(x=length,y=after_stat(density)),alpha=.8, fill='#6b6b6b', show.legend=F)
  p<-p+geom_density(mapping=aes(x=length, fill=fill_col),col=opdf_colour,lwd=1)
  p<-p+scale_fill_manual(name='',values=c(adjustcolor(opdf_colour,.5)))
  p<-p+geom_line(data=pdfs,mapping=aes(x=x,y=y, group=component, color=component),lwd=1.5)
  p<-p+scale_color_manual(name='Mixture\ncomponent',values=component_colors[1:max(spp_data$n_components)])
  p<-p+geom_vline(data=spp_data,mapping=aes(xintercept=mu,color=component),lty=2,lwd=1, show.legend=F)
  p<-p+geom_text(data=spp_data,mapping=aes(label=paste0('Dip test p: ', round(dip_test_p,2), dip_sig)),x=0,y=Inf, hjust=0,vjust=1, size=dip_test_size)
  
  p<-p+xlab('Body length [mm]') + labs(spp)
  return(p)
}

translate<-function(vector,named_translation_vector,na_replace=NULL){
  for (i in seq_along(named_translation_vector)){
    item <- named_translation_vector[i]
    original <- names(named_translation_vector)[i]
    vector[which(vector==original)]<-item
  }
  if (!is.null(na_replace)){
    vector[is.na(vector)]<-na_replace
  }
  return(vector)
}

create_currents_raster<-function(currents_netcdf_file, date, 
                              depth_layer=1,
                              output_folder=NULL, 
                              create_if_not_exist=F, 
                              force_return=F,
                              overwrite=F,
                              var_lon = "longitude",
                              var_lat = "latitude",
                              var_u = "uo",
                              var_v = "vo"
){
  require(terra)
  require(ncdf4)
  
  date<-as.Date(date)
  year<-as.numeric(format(date,"%Y"))
  jday<-as.numeric(format(date,"%j"))
  
  print(paste0('Extracting currents data for ',date))
  
  if (!is.null(output_folder)){
    filename<-file.path(output_folder, paste0(date, '_at_layer_',depth_layer,'.tif'))
    if (file.exists(filename) & !overwrite){
      print('Warnig! Using existing version of file.')
      return(terra::rast(filename)[[c('direction','speed')]])
    }
  }
  
  nc_data <- ncdf4::nc_open(currents_netcdf_file)
  lon <- ncdf4::ncvar_get(nc_data, var_lon)
  lat <- ncdf4::ncvar_get(nc_data, var_lat)
  u   <- ncdf4::ncvar_get(nc_data, var_u)
  v   <- ncdf4::ncvar_get(nc_data, var_v)
  ncdf4::nc_close(nc_data)
  
  # for values, it's in slices:
  # first 2 slices are the values, 3rd is depth band 4th is julian day
  uo<-u[,,depth_layer,jday]
  vo<-v[,,depth_layer,jday]
  
  grid <- expand.grid(lon = lon, lat = lat)
  grid$u<-as.vector(uo)
  grid$v<-as.vector(vo)
  
  grid$speed <- sqrt(grid$u^2 + grid$v^2)  # Calculate current speed
  grid$direction<-(atan2(grid$u,grid$v)*180/pi+360) %% 360 # Calculate current direction
  r <- terra::rast(grid, type='xyz',crs='epsg:4326')
  
  if (!is.null(output_folder)){
    if (!file.exists(output_folder) & create_if_not_exist){
      print(paste0('will create output folder at: ',output_folder))
      dir.create(output_folder,recursive=T, showWarnings = F)
    }
    terra::writeRaster(r[[c('direction','speed')]], filename = filename, overwrite=overwrite)
    if (force_return){
      return(r[[c('direction','speed')]])
    }
  }else{
    return(r)
  }
}

stars.pval <- function(x,lvls=c(0, 0.01, 0.05, 0.10, 1),stars = c("***", "**", "*", "")){
  i <- findInterval(x, lvls, left.open = T, rightmost.closed = T)
  return(stars[i])
}