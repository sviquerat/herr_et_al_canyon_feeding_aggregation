rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #only works in rstudio / posit
source(file.path(getwd(),'_SRC','_startup.R'))

#### ----- Clustering ---- #####
max_clusters=10

GFX_CLUSTER<-file.path(GFX_DIR,'CLUSTER')
dir.create(GFX_CLUSTER,showWarnings = F)

load(MSM115_PCA_DATA)
load(MSM115_KRILL_COMPLETE_DATA)

ENV_STACK<-terra::rast(ENV_STACK_FILE)
OCE_STACK<-terra::rast(OCE_STACK_FILE)

stations<-KRILL_FINAL$abu_age

# add covars from oce and topo stack and done, maybe plotting
stations<-sf::st_as_sf(stations,coords=c('x','y'),crs=IBCSO,remove=F)
v<-terra::extract(ENV_STACK[[-c(1:2)]],stations)
v<-cbind(v,terra::extract(OCE_STACK,stations))
stations<-sf::st_drop_geometry(cbind(stations,v))

ids_original<-MSM115_PCA$id_columns
axis<-MSM115_PCA$axis

s_pca_dat<-MSM115_PCA$pca_dat_transformed
table(s_pca_dat$recordType)
s_pca_dat<-s_pca_dat[,-which(names(s_pca_dat) %in% names(ids_original))]
s_pca_dat<-s_pca_dat[,axis]

n_clusters<-NbClust::NbClust(data = s_pca_dat, distance = "euclidean", min.nc = 2, max.nc = max_clusters, method = "kmeans")

bars<-as.data.frame(t(n_clusters$Best.nc))
names(bars)<-c('n_clusters','value_idx')
bars$index<-rownames(bars)
bars$n_clusters_f<-as.factor(bars$n_clusters)

f<-ggplot(data=subset(bars,!is.na(value_idx) & n_clusters>0),aes(x=n_clusters_f,fill=n_clusters))
f<-f+geom_bar(stat='count')+xlab('Proposed number of clusters')+ylab('Frequency among all indices')+guides(fill='none')
f<-f+labs(title='Optimal number of clusters',subtitle='Pooled frequency of clusters')+ labs(tag='A')+MAINTHEME

png(file.path(PUB_PACKAGE_DIR,'supp_Figure_6.png'),3500,3500,res=300)
print(f+labs(tag=''))+theme_large_font()
graphics.off()

max_clusters<-6 #after checking the pooled vote, it doesn't make much sense to go higher

#### fuzzy c - means ####
model_idx<-1:(max_clusters-1)

fanny_clusters<-lapply(2:max_clusters,function(i)cluster::fanny(s_pca_dat,k=i,metric='SqEuclidean',maxit=1000))
cmeans_clusters<-lapply(2:max_clusters,function(i)e1071::cmeans(s_pca_dat, centers=i, method="cmeans", dist='euclidean', iter.max=1000))

dunn_fanny = sapply(model_idx,function(i)Dunn(fanny_clusters[[i]]))
dunn_cmeans = sapply(model_idx,function(i)Dunn(cmeans_clusters[[i]]))

dunn_data<-data.frame(n_clusters=rep(2:max_clusters,2),
                      dunn=c(dunn_fanny[2,],dunn_cmeans[2,]),
                      algorithm=c(rep('fanny',max_clusters-1),rep('c-means',max_clusters-1)),
                      version='regular')

dunn_data<-rbind(dunn_data,data.frame(n_clusters=rep(2:max_clusters,2),
                                      dunn=c(dunn_fanny[1,],dunn_cmeans[1,]),
                                      algorithm=c(rep('fanny',max_clusters-1),rep('c-means',max_clusters-1)),
                                      version='normalised')
)

p<-ggplot(data=dunn_data,aes(x=n_clusters,y=dunn,col=algorithm))
p<-p+geom_line(lwd=1.5) + labs(y="Dunn's coefficient (normalised)",x='Number of clusters')
p<-p+MAINTHEME+facet_wrap(.~version)

png(file.path(GFX_CLUSTER,'MSM115_fuzzy_cluster_dunn_coeff.png'),2500,2500,res=300)
print(p+theme_large_font())
graphics.off()

#### final clustering
k_clusters<-6

cm<-cmeans_clusters[[k_clusters-1]]
fa<-fanny_clusters[[k_clusters-1]]

cluster_data<-data.frame(
  cm=cm$cluster,
  fa=fa$cluster
)

ids<-cbind(ids_original,cluster_data)

tables<-list(
  cm = table(ids$cm,ids$fancy_recordType),
  fa = table(ids$fa,ids$fancy_recordType)
)

### forced reclassification:
#group names now based on proportion of total number of aggregations

for (method in names(tables)){
  tbl<-tables[[method]]
  i<-which(colnames(tbl) == 'Aggregation')
  j<-which(colnames(tbl) == 'Solitary')
  k<-which(colnames(tbl) == 'Station')
  
  prop_agg<-tbl[,i]/sum(tbl[,i])
  prop_sin<-tbl[,j]/sum(tbl[,j])
  prop_sta<-tbl[,k]/sum(tbl[,k])
  vs<-data.frame(n_agg=tbl[,i],p_agg=prop_agg,n_single=tbl[,j],p_single=prop_sin,n_station=tbl[,k],p_stat=prop_sta)
  vs[[method]]<-as.numeric(rownames(tbl))
  vs$pretty_cluster_name_long<-paste0("A: ",round(vs$p_agg,4)*100,'% | S: ',round(vs$p_single,4)*100,'%')
  vs$pretty_cluster_name<-'nondescript'
  vs$pretty_cluster_name[vs$n_agg==max(vs$n_agg)]<-'aggregations'
  method_pretty<-paste0(method,'_pretty')
  cluster_data[[method_pretty]]<-dplyr::left_join(cluster_data, vs,by=c(method))$pretty_cluster_name
}

ids<-cbind(ids_original,cluster_data)

members<-cm$membership

idx<-which(ids$recordType == 'aggr')
m<-as.matrix(members)[idx,]
rownames(m)<-1:nrow(m)

png(file.path(GFX_CLUSTER,'MSM115_fuzzy_membership.png'),3200,3200,res=300)
p<-ggHeatmap(m, legend_name='fuzzy membership',sig_digits=4,removeTri=F,lower_limit=0,upper_limit=1,force_symmetry=F,low_col=low_colour,mid_col = mid_colour, high_col = high_colour)
p+labs(title = 'Fuzzy membership of in situ observations',subtitle = 'Aggregations')
graphics.off()

tables<-list(
  cm = table(ids$cm,ids$fancy_recordType),
  fa = table(ids$fa,ids$fancy_recordType)
)

cmeans_table_raw<-table(ids$cm,ids$fancy_recordType)
cmeans_table<-table(ids$cm_pretty,ids$fancy_recordType)
fanny_table<-table(ids$cm_pretty,ids$fancy_recordType)

l_raw<-ggHeatmap(cmeans_table_raw, removeTri=F,lower_limit=0,force_symmetry=F, low_col=low_colour,mid_col = mid_colour,high_col = high_colour)+labs(title='Fuzzy clustering',subtitle='C-Means')

n<-ggHeatmap(apply(cmeans_table,FUN = col_ratio,MARGIN=2), sig_digits=2,removeTri=F,lower_limit=0,upper_limit=100,force_symmetry=F,low_col=low_colour,mid_col = mid_colour,high_col = high_colour)+labs(title='Fuzzy clustering',subtitle='C-Means [%]')
o<-ggHeatmap(apply(fanny_table,FUN = col_ratio,MARGIN=2), sig_digits=2, removeTri=F,lower_limit=0,upper_limit=100,force_symmetry=F, low_col=low_colour,mid_col = mid_colour,high_col = high_colour)+labs(title='Fuzzy clustering',subtitle='Fanny [%]')

png(file.path(PUB_PACKAGE_DIR,'Figure 13.png'),3200,3200,res=300)
gridExtra::grid.arrange(grobs=list(l_raw+labs(tag='A')+scale_y_continuous(labels=1:k_clusters, breaks=1:k_clusters)+theme_large_font(12),
                                   n+labs(tag='B')+theme_large_font(16)),ncol=2)
graphics.off()

ids$cluster_name<-ids$cm_pretty
ids$cluster<-ids$cm

station_cluster<-subset(ids,recordType=='station')
station_cluster<-sqldf::sqldf('select * from stations as a left join station_cluster as b on a.station_id = b.id')
station_cluster$cluster_name<-factor(station_cluster$cluster_name,ordered=T)
d<-cmeans_table

krill<-grep('krill_.*',names(station_cluster))

#### should be fixed in krill_abundance script!
station_cluster[,krill][is.na(station_cluster[,krill])]<-0
station_cluster<-station_cluster[,colSums(station_cluster != 0,na.rm=T) > 0] # there are all 0 columns - we don't need them

table(station_cluster$cluster_name,station_cluster$cluster)
agg_cluster<-station_cluster$cluster[which(station_cluster$cluster_name == 'aggregations')][1]

classification_table<-data.frame(value=1:k_clusters, literal='nondescript', colors=cluster_colors[1:k_clusters])
classification_table$literal[agg_cluster]<-'aggregations'
openxlsx::write.xlsx(classification_table,file=MSM115_CLUSTER_CLASSIFICATION_FILE)

MSM115_CLUSTER<-list(ids_original=ids_original,posteriori=ids,k_clusters=k_clusters, models=list(cm=cm,fa=fa),chosen_model='cm',tables=tables,station_cluster=station_cluster,classification_table=classification_table,membership=members)
save(MSM115_CLUSTER,file = MSM115_CLUSTER_DATA,compress='gzip')

#### converting to geopackage ####
LL<-sf::st_as_sf(ids,coords=c('x','y'),crs=IBCSO)
LL$fancy_recordType<-as.factor(LL$fancy_recordType)
LL$cluster_name<-as.factor(LL$cluster_name)

##### data export #####
sf::st_write(LL,dsn=MSM115_CLUSTER_RESULTS,layer='clustered data',overwrite=T, delete_dsn=T)
