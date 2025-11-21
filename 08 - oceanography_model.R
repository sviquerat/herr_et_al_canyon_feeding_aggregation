rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source(file.path(getwd(),'_SRC','_startup.R'))

#### ----- Oceanography modelling ---- #####
GFX_OCE<-file.path(GFX_DIR,'OCE')
dir.create(GFX_OCE,recursive=T, showWarnings = F)

ENV_STACK<-terra::rast(ENV_STACK_FILE)

### modelling ###
library(mgcv)

load(OCE_DATA_FILE)
load(SOAP_PARAMETERS)

ud<-MSM115_OCE_DATA$unified_detrended[MSM115_OCE_DATA$unified_detrended$cast=='day',]
ud<-ud[complete.cases(ud),]

depth_name<-'depth_m_unified'
covars<-c('temp1_C', 'salinity','oxygen_ml_l','oxygen_saturation_perc','fluorescence_mg_m3')
processing_type<-'detrended' #filter1 includes detrending, but introduces NA due to moving windows
covars<-paste0(covars,'_',processing_type)

surface_limit=c(-5,-50,-200) #everything from 0 to ...
gradients<-sqldf::sqldf('select station_id,x,y from ud group by station_id;')

for (st in unique(ud$station_id)){
  for (var_name in covars){
    for (max_depth in surface_limit){
      row_idx<-which(ud$station_id==st & ud[[depth_name]] >= max_depth)
      y<-ud[row_idx, var_name]
      x<-ud[row_idx, depth_name]
      
      v_name<-gsub(paste0('_',processing_type),paste0('_',abs(max_depth),'m'),var_name)
      
      gradients[gradients$station_id==st,paste0(v_name,'_mean')]<-mean(y,na.rm = T)
      gradients[gradients$station_id==st,paste0(v_name,'_se')]<-sd(y,na.rm = T)/sqrt(length(y[!is.na(y)]))
    }
  }
}
gradients<-gradients[order(formatC(gradients$station_id,2,2,flag="0")),]
v<-terra::extract(ENV_STACK[[-c(1,2)]],sf::st_as_sf(gradients,coords=c('x','y'),crs=IBCSO))
gradients<-cbind(gradients, v[,-1])

averages<-names(gradients)[grep('*_mean',names(gradients))]

bound<-MSM115_SOAP_PARAMETERS$boundary
knots<-MSM115_SOAP_PARAMETERS$knots

p_stack<-ENV_STACK

model_list<-lapply(averages, function(var_name){
  print(paste0('Processing: ',var_name))
  m<-mgcv::gam(gradients[[var_name]] ~ s(x, y, bs = "so", xt = list(bnd = bound)), knots=knots, data=gradients, family=gaussian(link='identity'))
  
  p<-terra::predict(p_stack, m, type='response', se.fit=T)
  return(p$fit)
})

names(model_list)<-averages

OCE_STACK<-list()
for (var_name in averages){
  p<-model_list[[var_name]]
  OCE_STACK[[var_name]]<-p
}

OCE_STACK<-terra::rast(OCE_STACK)

names(OCE_STACK) <- averages

terra::writeRaster(OCE_STACK,OCE_STACK_FILE,overwrite=T)

#### plotting ####
OCE_STACK<-terra::rast(OCE_STACK_FILE)

g<-list()
for (var in c('temp1','salinity','oxygen_ml','oxygen_saturation','fluorescence')){
  temp<-grep(paste0(var,'_.*'),names(OCE_STACK))
  extra<-''
  if (grepl('saturation',var)){
    extra = ' saturation '
  }
  temp_stack<-OCE_STACK[[temp]]
  var<-names(temp_stack)[1]
  var<-gsub('salinity_.*','salinity [psu]', var)
  var<-gsub('_C_.*',' [°C]', var)
  var<-gsub('_ml_l.*',' [ml / l]', var)
  var<-gsub('_saturation_.*',' [%]', var)
  var<-gsub('_mg_m3_.*',' [mg / m³]', var)
  
  pure_var<-strsplit(var,' ')[[1]][1]
  
  names(temp_stack)<-c('5m','50m','200m')
  p<-ggplot()+tidyterra::geom_spatraster(data=temp_stack) +  scale_fill_distiller(palette = "RdBu")
  p<-p+geom_sf(data=MSM115_OCE_DATA$combined$stations,cex=.7,col=adjustcolor('black',.6))
  p<-p+facet_wrap(~lyr)+labs(fill=var) #+ theme(plot.margin = margin(1,2,1,1, "cm"))
  g[[var]]<-p+MAINTHEME + theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold")) + theme(strip.text = element_text(size = 12))
  
}

png(file.path(PUB_PACKAGE_DIR,paste0('Figure 4.png')),4000,3200,res=300)
gridExtra::grid.arrange(grobs = g,ncol=2,widths=c(1,1))
graphics.off()
