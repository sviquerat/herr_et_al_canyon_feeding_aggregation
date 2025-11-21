rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #only works in rstudio / posit
source(file.path(getwd(),'_SRC','_startup.R'))

#### ----- Clustering ---- #####
# require(gridExtra)
GFX_CLUSTER<-file.path(GFX_DIR,'CLUSTER')
GFX_ASS<-file.path(GFX_CLUSTER,'ASSOCIATION')
dir.create(GFX_ASS,showWarnings = F, recursive=T)

#### net trawl data association ####
load(file = MSM115_CLUSTER_DATA)

ENV_STACK<-terra::rast(ENV_STACK_FILE)
OCE_STACK<-terra::rast(OCE_STACK_FILE)

station_cluster<-MSM115_CLUSTER$station_cluster
k_clusters<-MSM115_CLUSTER$k_clusters

mm_data<-sf::st_read(MSM115_KRILL_MULTIMODALITY_GPKG,layer='detailed')
mm_data<-sf::st_drop_geometry(mm_data)

#filter out everything but salps and euphausiids
# salps don't have a shannon index, there's only one species....
species<-grep('H_.*krill|H_.*salp|.*krill.*|.*salp.*|.*Euphausia.*|.*Thysanoessa.*|.*Salpa.*',names(station_cluster))
oce<-which(names(station_cluster) %in% names(OCE_STACK))
env<-which(names(station_cluster) %in% names(ENV_STACK)[-c(1,2)]) #1 & 2 are x,y
static<-which(names(station_cluster) %in% c('station_id','x','y','cluster_name','cluster'))
keep<-c(static,species,oce,env)
stc<-station_cluster[,keep]
stc$cluster[stc$cluster==5]<-'5*'
stc$cluster<-as.factor(stc$cluster)

x<-sqldf::sqldf('select * from stc as a left join mm_data as b on a.station_id = b.station_id')
mm_cols<-grep('is_multimodal.*',names(x))
x[is.na(x)]<-0
cols<-names(x)[mm_cols]
stc[,cols]<-x[,cols]
stc[,cols][stc[,cols]=='1']<-'yes'
stc[,cols][stc[,cols]=='0']<-'no'

H<-grep('H_.*',names(stc))
multimodal<-grep('is_multimodal.*',names(stc))
ratios<-grep('.*_ratio',names(stc))
krill<-grep('krill_.*',names(stc))
oce<-which(names(stc) %in% names(OCE_STACK))
env<-which(names(stc) %in% names(ENV_STACK)[-c(1,2)]) #1 & 2 are x,y

stc[,multimodal]<-apply(stc[,multimodal],2,FUN='as.factor')

df<-stc[,c('cluster','cluster_name',names(stc)[multimodal])]
df<-reshape2::melt(df, id.vars=c('cluster','cluster_name'))
df$variable<-gsub('is_multimodal_','',df$variable)
df$variable<-gsub('\\.',' ',df$variable)
df$variable[df$variable == 'is_multimodal']<-'Any species'
df$variable<-factor(df$variable, levels=c('Euphausia superba','Euphausia frigida','Euphausia triacantha','Thysanoessa macrura','Any species'))

#### all clusters
fac_data<-sqldf::sqldf('select cluster, variable as species, COUNT(*) AS total, SUM(CASE WHEN value = "yes" THEN 1 ELSE 0 END) AS multimodal_yes, 0 as multimodal_no from df group by cluster, species order by species, cluster')
fac_data$multimodal_no<-fac_data$total-fac_data$multimodal_yes

out<-data.frame(cluster=unique(fac_data$cluster))
for (spp in unique(fac_data$species)){
  dat<-fac_data[fac_data$species==spp,]
  new_names<-paste0(spp,'_',names(dat)[3:5])
  new_names_relative<-paste0(spp,'_r_',names(dat)[3:5])
  dat[,new_names]<-dat[,names(dat[3:5])]
  dat[,new_names_relative]<-dat[,names(dat[3:5])]/dat[,3]
  dat<-dat[,-c(2:5)]
  out<-sqldf::sqldf('select * from out as a left join dat as b on a.cluster=b.cluster')
  kill<-which(names(out) == 'cluster')
  if (length(kill)>1){
    out<-out[,-kill[2]]
  }
}
openxlsx::write.xlsx(out, MSM115_KRILL_MULTIMODALITY_CLUSTER_FILE)

fac_data<-sqldf::sqldf('select cluster_name, variable as species, COUNT(*) AS total, SUM(CASE WHEN value = "yes" THEN 1 ELSE 0 END) AS multimodal_yes, 0 as multimodal_no from df group by cluster_name, species order by species, cluster_name')
fac_data$multimodal_no<-fac_data$total-fac_data$multimodal_yes

out<-data.frame(cluster_name=unique(fac_data$cluster_name))
for (spp in unique(fac_data$species)){
  dat<-fac_data[fac_data$species==spp,]
  new_names<-paste0(spp,'_',names(dat)[3:5])
  new_names_relative<-paste0(spp,'_r_',names(dat)[3:5])
  dat[,new_names]<-dat[,names(dat[3:5])]
  dat[,new_names_relative]<-dat[,names(dat[3:5])]/dat[,3]
  dat<-dat[,-c(2:5)]
  out<-sqldf::sqldf('select * from out as a left join dat as b on a.cluster_name=b.cluster_name')
  kill<-which(names(out) == 'cluster_name')
  if (length(kill)>1){
    out<-out[,-kill[2]]
  }
}
openxlsx::write.xlsx(out, MSM115_KRILL_MULTIMODALITY_CLUSTER_NAME_FILE)

### association before reclassification of clusters
a<-gg_pairs(num_dat=stc,data_scale=F,groups = stc$cluster,columns=H,group_name='Assigned cluster',force_single_plot = F,
            colours=cluster_colors,show_outliers=F, theme = MAINTHEME)
b<-gg_pairs(num_dat=stc,data_scale=F,groups = stc$cluster,columns=ratios,group_name='Assigned cluster',force_single_plot = F,
            colours=cluster_colors,show_outliers=F, theme = MAINTHEME)
c<-gg_pairs(num_dat=stc,data_scale=F,groups = stc$cluster,columns=krill,group_name='Assigned cluster',force_single_plot = F,
            colours=cluster_colors,show_outliers=F, theme = MAINTHEME)
d<-gg_pairs(num_dat=stc,data_scale=F,groups = stc$cluster,columns=oce,group_name='Assigned cluster',force_single_plot = F,
            colours=cluster_colors,show_outliers=F, theme = MAINTHEME)
e<-gg_pairs(num_dat=stc,data_scale=F,groups = stc$cluster,columns=env,group_name='Assigned cluster',force_single_plot = F,
            colours=cluster_colors,show_outliers=F, theme = MAINTHEME)

f<-ggally_count(data=df,aes(x=value,y=cluster,fill=cluster,group=cluster),show.legend=F)
f<-f+scale_fill_manual(values=cluster_colors)+labs(x='Multimodality',y='Assigned cluster')
f<-f+facet_grid(.~variable)+MAINTHEME

plots<-c('diversity'=a,'ratios'=b,'krill'=c,'oce'=d,'topo'=e,'multimodal'=1)

for (i in seq_along(plots)){
  name<-names(plots)[i]
  print(name)
  
  if (name=='multimodal'){
    p<-ggally_count(data=df,aes(x=value,y=cluster_name,fill=cluster_name,group=cluster_name),show.legend=F)
    p<-p+scale_fill_manual(values=cluster_name_colors)+labs(x='Multimodality',y='Assigned cluster')
    p<-p+facet_grid(.~variable)+MAINTHEME
    
    filename<-file.path(GFX_ASS,paste0('MSM115_all_cluster_vs_',name,'.png'))
    png(filename,2500,2500,res=300)
    print(p)
    graphics.off()    
    next()
  }
  
  p<-plots[i]
  n<-sqrt(length(p[[1]]$plots))
  filename<-file.path(GFX_ASS,paste0('MSM115_all_cluster_vs_',names(plots)[i],'.png'))
  
  png(filename,n*700,n*700,res=300)
  print(p)
  graphics.off() 
  
}

#### reassigned clusters
a<-gg_pairs(num_dat=stc,data_scale=F,groups = stc$cluster_name,columns=H,group_name='Assigned cluster',force_single_plot = F,
            colours=cluster_name_colors,show_outliers=F, theme = MAINTHEME)
b<-gg_pairs(num_dat=stc,data_scale=F,groups = stc$cluster_name,columns=ratios,group_name='Assigned cluster',force_single_plot = F,
            colours=cluster_name_colors,show_outliers=F, theme = MAINTHEME)
c<-gg_pairs(num_dat=stc,data_scale=F,groups = stc$cluster_name,columns=krill,group_name='Assigned cluster',force_single_plot = F,
            colours=cluster_name_colors,show_outliers=F, theme = MAINTHEME)
d<-gg_pairs(num_dat=stc,data_scale=F,groups = stc$cluster_name,columns=oce,group_name='Assigned cluster',force_single_plot = F,
            colours=cluster_name_colors,show_outliers=F, theme = MAINTHEME)
e<-gg_pairs(num_dat=stc,data_scale=F,groups = stc$cluster_name,columns=env,group_name='Assigned cluster',force_single_plot = F,
            colours=cluster_name_colors,show_outliers=F, theme = MAINTHEME)

plots<-c('diversity'=a,'ratios'=b,'krill'=c,'oce'=d,'topo'=e,'multimodal'= 1)

for (i in seq_along(plots)){
  name<-names(plots)[i]
  print(name)
  
  if (name=='multimodal'){
    p<-ggally_count(data=df,aes(x=value,y=cluster_name,fill=cluster_name,group=cluster_name),show.legend=F)
    p<-p+scale_fill_manual(values=cluster_name_colors)+labs(x='Multimodality',y='Assigned cluster')
    p<-p+facet_grid(.~variable)+MAINTHEME
    
    filename<-file.path(GFX_ASS,paste0('MSM115_cluster_vs_',name,'.png'))
    png(filename,2500,2500,res=300)
    print(p)
    graphics.off()    
    next()
  }
  
  p<-plots[i]
  n<-sqrt(length(p[[1]]$plots))
  filename<-file.path(GFX_ASS,paste0('MSM115_cluster_vs_',names(plots)[i],'.png'))
  
  png(filename,n*700,n*700,res=300)
  print(p)
  graphics.off() 

}

# focussed on post clustering selection
columns<-c("oxygen_ml_l_200m_mean","oxygen_saturation_perc_200m_mean","fluorescence_mg_m3_200m_mean","salinity_50m_mean", names(stc)[grep('temp1_*',names(stc))])
keys<-LETTERS[1:length(columns)]
names(keys)<-columns
df<-stc
names(df)[names(df) %in% columns]<-keys
df<-df[,which(names(df) %in% c(keys, 'cluster_name'))]
df<-cbind(df,stc[,columns])
df$cluster<-stc$cluster

legend<-data.frame(key=keys,variable=names(keys))
row.names(legend)<-NULL
write.table(legend, file=file.path(PUB_PACKAGE_DIR,'Figure 14 legend.txt'), sep='\t',row.names = F)

a<-gg_pairs(num_dat=df,data_scale=F,groups = df$cluster_name,columns=keys,group_name='Assigned cluster',force_single_plot = T,
            colours=cluster_name_colors,show_outliers=F, theme = MAINTHEME+theme_large_font(16),point.size=4,diag.size=5)

b<-gg_pairs(num_dat=df,data_scale=F,groups = df$cluster,columns=keys,group_name='Assigned cluster',force_single_plot = T,
            colours=cluster_colors,show_outliers=F, theme = MAINTHEME+theme_large_font(16),point.size=4,diag.size=5)

png(file.path(PUB_PACKAGE_DIR,'Figure 14.png'),7000,5000,res=300)
print(a)
graphics.off()

png(file.path(PUB_PACKAGE_DIR,'supp_Figure_7.png'),5000,5000,res=300)
print(b)
graphics.off()
