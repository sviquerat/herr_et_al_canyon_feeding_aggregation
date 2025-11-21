rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #only works in rstudio / posit
source(file.path(getwd(),'_SRC','_startup.R'))

#### ----- Combine net trawl and cetacean data with oceanography and geomorphology ---- #####
ENV_STACK<-terra::rast(ENV_STACK_FILE)
OCE_STACK<-terra::rast(OCE_STACK_FILE)

env_col_names<-names(ENV_STACK)[-c(1,2)] # 1,2 are x,y
oce_col_names<-names(OCE_STACK)

load(MSM115_CET_ANALYSIS_DATA)
load(MSM115_KRILL_COMPLETE_DATA)
load(OCE_DATA_FILE)

sigs<-MSMS115_CET$sig_data
krill<-KRILL_FINAL$abu_age

# add covars from oce and topo stack and done, maybe plotting
KRILL<-sf::st_as_sf(krill,coords=c('x','y'),crs=IBCSO,remove=F)
v<-terra::extract(ENV_STACK[[-c(1:2)]],KRILL)
v<-cbind(v,terra::extract(OCE_STACK,KRILL))
KRILL<-sf::st_drop_geometry(cbind(KRILL,v))

CET<-sf::st_as_sf(sigs,coords=c('x','y'),crs=IBCSO,remove=F)
v<-terra::extract(ENV_STACK[[-c(1:2)]],CET)
v<-cbind(v,terra::extract(OCE_STACK,CET))
CET<-sf::st_drop_geometry(cbind(CET,v))

keep_vars<-which(names(KRILL) %in% c('x','y',env_col_names,oce_col_names))
pca_st<-KRILL[,c(1,keep_vars)]
keep_vars<-which(names(CET) %in% c('id','species','x','y','date',env_col_names,oce_col_names))
pca_sigs<-CET[,keep_vars]

pca_st<-cbind(pca_st$station_id,rep(NA,nrow(pca_st)),rep('station',nrow(pca_st)),pca_st[,2:length(names(pca_st))])
names(pca_st)[1:3]<-names(pca_sigs)[1:3]
pca_dat<-rbind(pca_sigs,pca_st)
names(pca_dat)[3]<-'recordType'

#check distances to stations
ctd_stations<-MSM115_OCE_DATA$unified_detrended[,1:3] # this should go somewhere?
ctd<-sqldf::sqldf('select * from ctd_stations group by station_id, "cast"')
PP<-sf::st_as_sf(pca_dat,coords=c('x','y'),crs=IBCSO)
CTD<-sf::st_as_sf(ctd,coords=c('x','y'),crs=IBCSO)
distances<-sf::st_distance(PP,CTD)
ndist<-c()
for (j in 1:nrow(distances)){
  x<-min(distances[j,])
  ndist<-c(ndist,x)
}
dist_break<-.9
max_dist<-quantile(ndist,dist_break)
INCLUSION_ZONE<-sf::st_buffer(CTD,dist = max_dist)
INCLUSION_ZONE<-sf::st_union(INCLUSION_ZONE)
j<-sf::st_intersects(PP,INCLUSION_ZONE,sparse = F)
PP$include<-as.numeric(j[,1])
INCLUSION_ZONE<-sf::st_as_sf(INCLUSION_ZONE)

pca_dat$include<-as.numeric(j[,1])
pca_dat$fancy_recordType<-'Solitary'
pca_dat$fancy_recordType[pca_dat$recordType=='aggr']<-'Aggregation'
pca_dat$fancy_recordType[pca_dat$recordType=='station']<-'Station'
table(pca_dat$recordType,pca_dat$include)

xy_cols<-which(names(pca_dat) %in% c('x','y'))
id_columns<-names(pca_dat)[c(1:3,xy_cols,which(names(pca_dat)%in%c('include','fancy_recordType')))]

MSM115_FULL<-list(pca_dat=pca_dat, id_columns=id_columns,xy_cols=xy_cols,krill=KRILL, sigs=CET)
save(MSM115_FULL, file=MSM115_MERGED_CET_NET_OCE_TOPO,compress='gzip')

#### plotting ####
require(ggplot2)
exclude_color='#a40f0a'
include_color='#2a689e'
buffer_color<-'#ffb602'

PP$pretty_include<-'excluded'
PP$pretty_include[PP$include==1]<-'included'
PP$fancy_recordType<-'Solitary'
PP$fancy_recordType[PP$recordType=='aggr']<-'Aggregation'
PP$fancy_recordType[PP$recordType=='station']<-'Station'

M<-table(PP$fancy_recordType,PP$pretty_include)
M<-as.data.frame(M)
md<-format(round(as.numeric(max_dist), 1), nsmall=0, big.mark=",")
coords<-data.frame(x=c(max_dist,max_dist,0,max_dist),y=c(0,rep(dist_break,3)), grp=c(rep(md,2),rep(dist_break,2)))
coords$grp<-as.factor(coords$grp)

p<-ggplot(data.frame(x=ndist),aes(x=x))+stat_ecdf()
p<-p+geom_line(data=coords,aes(x=x,y=y,group=grp,col=grp),lwd=1.5,alpha=1)
p<-p+geom_point(data=coords[2,],aes(x=x,y=y),col='orange',alpha=1,cex=3)
p<-p+scale_color_manual(values=c(exclude_color,buffer_color))
p<-p + scale_y_continuous(breaks=c(0,.25,.5,.75,.9,1.),expand = c(0, 0), limits = c(0, NA))
p<-p + scale_x_continuous(expand = c(0, 0), limits = c(0, NA))
p<-p+labs(col='',x='maximum distance to CTD station [m]',y='ecdf')+MAINTHEME
p<-p+theme(legend.position='inside',legend.justification.inside = c(1, 0),legend.background=element_blank())

g<- ggplot() + tidyterra::geom_spatraster(data=ENV_STACK$depth)
g<- g + scale_fill_gradientn( colours = adjustcolor(pal_depth(10),.7), na.value = "transparent", guide='none')
g<-g+geom_sf(data=INCLUSION_ZONE, fill=buffer_color)
g<-g+geom_sf(data=PP,mapping=aes(col=pretty_include), show.legend=T)+scale_colour_manual(values=c(exclude_color,include_color))
g<-g+MAINTHEME
g<-g+labs(col='Data')
g<-g+theme(legend.position='inside',legend.justification.inside = c(1, 0),legend.background=element_blank())

f<-ggplot(M[M$Freq>0,], aes(x=Var1, y = Freq, fill=Var2)) + geom_bar(stat="identity")
f<-f+scale_fill_manual(values=c(exclude_color,include_color))+labs(x='',y='',fill='')
f<-f+geom_text(aes(x = Var1, y = Freq, label = Freq, group = Var2),position = position_stack(vjust = .5),col='white')
f<-f+scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
f<-f+MAINTHEME+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank()) 
f<-f+theme(legend.position='inside',legend.justification.inside = c(0, 1),legend.background=element_blank())

m<-rbind(c(1, 2),
         c(1, 3))

png(file.path(PUB_PACKAGE_DIR,'supp_Figure_2.png'),3200,3200,res=300)
gridExtra::grid.arrange(grobs=list(g+labs(tag='A')+theme_large_font(),f+labs(tag='B')+theme_large_font(),p+labs(tag='C')+theme_large_font()),layout_matrix = m)
graphics.off()
