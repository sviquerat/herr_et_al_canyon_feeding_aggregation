rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #only works in rstudio / posit
source(file.path(getwd(),'_SRC','_startup.R'))

#### ----- env stack creation ---- #####
keys<-c('slope','aspect',  'roughness')
depth<-terra::rast(DEPTH_FILE)

rasterfiles<-list(depth=depth)
for (key in keys){
  rasterfiles[[key]]<-terra::terrain(depth,v=key, unit='radians')
}
rasterfiles[['dist2shelf_km']]<-terra::rast(DIST2SHELF_FILE)

x<-depth
coords<-terra::xyFromCell(x,1:terra::ncell(x))
x[1:terra::ncell(x)]<-coords[,1]
y<-x
y[1:terra::ncell(x)]<-coords[,2]

env_stack<-c(x,y)
for (r in rasterfiles){
  env_stack<-c(env_stack,r)
}
names(env_stack)<-c('x','y',names(rasterfiles))

mask_keys<-cenv_stackmask_keys<-c('shelf_mask', 'island_mask')
rasterfiles<-list()
mask_stack<-c()
for (key in mask_keys){
  rasterfiles[[key]]<-list.files(IBCSO_DATA_DIR,paste0('MSM115_IBCSO_',key,'.tif$'), full.names = T)
  mask_stack<-c(mask_stack,terra::rast(rasterfiles[[key]]))
}
names(mask_stack)<-names(rasterfiles)
mask_stack<-terra::rast(mask_stack)
mask_stack[mask_stack==0]<-NA

env_stack<-env_stack*mask_stack$island_mask

ENV_STACK<-env_stack
MASK_STACK<-mask_stack

ENV_COORDS_COVARS<-names(ENV_STACK)[1:2]
ENV_COVARS<-names(ENV_STACK)[3:length(names(ENV_STACK))]

##### save stuff
terra::writeRaster(ENV_STACK,ENV_STACK_FILE, overwrite=T)
terra::writeRaster(MASK_STACK,ENV_STACK_MASK_FILE, overwrite=T)
MSM115_ENV_DATA<-list('coords_covars'=ENV_COORDS_COVARS, 'covars'= ENV_COVARS)
save(MSM115_ENV_DATA, file=ENV_DATA_FILE, compress='gzip')
