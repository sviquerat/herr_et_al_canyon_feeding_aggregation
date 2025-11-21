rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #only works in rstudio / posit
source(file.path(getwd(),'_SRC','_startup.R'))

#### ----- PCA dimensionality reduction ---- #####
ENV_STACK<-terra::rast(ENV_STACK_FILE)
OCE_STACK<-terra::rast(OCE_STACK_FILE)

env_col_names<-names(ENV_STACK)[-c(1,2)] # 1,2 are x,y
oce_col_names<-names(OCE_STACK)

load(MSM115_MERGED_CET_NET_OCE_TOPO)
pca_dat<-MSM115_FULL$pca_dat

#subset based on distance to ctds
pca_dat<-pca_dat[pca_dat$include==T,]
table(pca_dat$recordType)

ids_original<-pca_dat[,MSM115_FULL$id_columns]
pca_dat_o<-pca_dat[,!(names(pca_dat) %in% MSM115_FULL$id_columns)]

# new thing: all oxygen surface data are constant (model result)
kill<-c('oxygen_ml_l_5m_mean','oxygen_saturation_perc_5m_mean')
pca_dat<-pca_dat_o[,-which(names(pca_dat_o) %in% kill)]

oce_cols<-which(names(pca_dat) %in% oce_col_names)
env_cols<-which(names(pca_dat) %in% env_col_names)

# dim reduction
reduction<-brute_force_dropper_pca(pca_dat)
reduction$data

png(file.path(PUB_PACKAGE_DIR,'supp_Figure_3.png'),3200,2200,res=300)
reduction$p + theme_large_font()
graphics.off()

drops<-c('depth','aspect','roughness')
pca_dat<-pca_dat[,-which(names(pca_dat) %in% drops)]
reduction<-brute_force_dropper_pca(pca_dat)
reduction$p
ggTotalLoading(pca_dat,axis=1:3,sorted=T,colours=pal_cluster(ncol(pca_dat)))

#### final
axis<-1:3

pca_res<-gg_pca_summary(pca_dat,add_cor=T)

barcol<-colorRampPalette(c('#5571a9','#8badf3'))
n_comp<-min(10,length(pca_res$sdev))
var<-pca_res$sdev^2
total_var<-sum(var)
percentage_explained_variance <-( (var / total_var) * 100)[1:n_comp]
pexpv<-paste0(round(percentage_explained_variance,2),'%')

p<-factoextra::fviz_screeplot(pca_res, geom='bar',addlabels = F, ncp=n_comp,barcol=barcol(n_comp),barfill=barcol(n_comp))
p<-p+geom_line(col='orange',lwd=2)#+geom_point(col='black',size=4)
p<-p+geom_text(vjust=-1,label=pexpv)
p<-p+labs(title=NULL)+MAINTHEME + theme_large_font()

png(file.path(PUB_PACKAGE_DIR,'supp_Figure_4.png'),3200,3200,res=300)
print(p)
graphics.off()

png(file.path(PUB_PACKAGE_DIR,'Figure 12.png'),3500,3500,res=300)
gg_pca_contrib_all(pca_dat,axes=1:3, addlabels = TRUE)+MAINTHEME+theme_large_font()
graphics.off()

p<-ggAxisLoads(pca_dat,axis=axis, which='axis',colours = pal_cluster(15))
p<-p+MAINTHEME+theme(legend.position='right')+theme_large_font()

g<-ggTotalLoading(pca_dat,axis=axis,sorted=T,colours=pal_cluster(ncol(pca_dat)))
g<-g + MAINTHEME + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))
g<-g + theme_large_font()+coord_flip()

m<-matrix(c(1,1,2,1,1,2),ncol=2)
grobs<-gridExtra::arrangeGrob(grobs=list(p+labs(tag='A'),g+labs(tag='B')), layout_matrix=m)

png(file.path(PUB_PACKAGE_DIR,'supp_Figure_5.png'),4400,4400,res=300)
gridExtra::grid.arrange(grobs)
graphics.off()

s_pca_dat<-scale(pca_dat, center=T, scale=T)
centers<-attr(s_pca_dat, 'scaled:center')
scales<-attr(s_pca_dat, 'scaled:scale')

pca_res<-prcomp(s_pca_dat,scale=F,center=F)
imp<-summary(pca_res)
prediction<-predict(pca_res)
pca_dat_transformed<-cbind(ids_original,prediction)

MSM115_PCA<-list(pca_dat=pca_dat, s_pca_dat=s_pca_dat, centers=centers, scales=scales, pca_dat_transformed=pca_dat_transformed,axis=axis,importance=imp,id_columns=ids_original,pca_res=pca_res)
save(MSM115_PCA, file=MSM115_PCA_DATA,compress='gzip')
