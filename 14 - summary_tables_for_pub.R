rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #only works in rstudio / posit
source(file.path(getwd(),'_SRC','_startup.R'))

# In case we'll need to add a stratum to a station...
# join via station_id
stations<-sf::st_read(MSM115_SPATIAL_DATA, layer='MSM115_stations') # Just the coordinates and id of each station

STRATA<-sf::st_read(MSM115_STRATA_FILE,layer='simplified_strata')
stations$stratum<-sf::st_join(stations,STRATA)$name
stations$stratum[is.na(stations$stratum)]<-'Between strata'
stations$stratum<-factor(stations$stratum,levels=c('South Shetland Islands (SSI)','Elephant Island (EI)','South Orkney Island (SOI)','Between strata'))
stations<-sf::st_drop_geometry(stations)

#### ----- Table creation for publication ---- ####

#### Table 1 detection function bulk settings ####
table<-openxlsx::read.xlsx(MSM115_DS_MODEL_COV_XLSX)
w<-unique(table$truncation)
w<-formatC(w,big.mark=',',format='f',drop0trailing = T,digits=0)
w<-paste0(w,collapse='\n')

adj<-unique(table$adjustment)
translation<-c(cos='cosine', herm='hermite',poly='polynomial')
adj<-translate(adj,translation,na_replace='NULL')
adj<-paste0(adj,collapse = '\n')

cov<-unique(table$formula)
translation<-c(glare_combined='glare', seastate_combined='sea state',subj_combined='subjective conditions')
cov<-translate(cov,translation,na_replace='NULL')
cov<-paste0(cov,collapse='\n')

keys<-unique(table$key)
translation<-c(hn='Half normal', hr='Hazard-rate',unif='Uniform')
keys<-translate(keys,translation)
table<-data.frame(Truncation_width_m=rep(w,length(keys)), Adjustments=rep(adj,length(keys)), Model_covariates=c(rep(cov,length(keys)-1),'None'))
table$Keys<-keys
openxlsx::write.xlsx(table,file.path(PUB_PACKAGE_DIR,'Table 1.xlsx'))

#### Table 2 is manual ####

#### Table 3 cetacean sightings ####
# TAKE CARE - There was one aggregation recorded during the visual survey.
# see script 05 - cet_prep.R ll. 21
# We decided to remove them from the visual survey and put them into the opp sightings data
# Therefore we'll need to add it back into the table, becuse it will not get separated correctly
new_line<-data.frame(survey='dedicated',species='aggregation',G=15,I=20)

table<-openxlsx::read.xlsx(MSM115_CET_SIGHTINGS)
table$data_source[table$data_source %in% c('opp_sightings','aggregations')]<-'opportunistic'
table$data_source[table$data_source!='opportunistic']<-'dedicated'
table<-table[table$species %in% c('aggr', 'bphy'),]
table$species[table$species=='aggr']<-'aggregation'
table$species[table$species=='bphy']<-'solitary group'
names(table)[1]<-'survey'
table<-rbind(table,new_line)
openxlsx::write.xlsx(table,file.path(PUB_PACKAGE_DIR,'Table 3.xlsx'))
# In the publication, we used more approximate numbers for aggregations since the numbers in situ could only be estimated 

#### Table 4 distance candidate models ####
table<-openxlsx::read.xlsx(MSM115_CET_DS_CANDIDATE_MODELS)
table<-sqldf::sqldf("select modelname as name, key, formula, truncation as w, CvM_p, CvM_w, pa_hat, pa_hat_se, N_sigs as G from 'table';")
table$w<-formatC(table$w, format='f',big.mark=',', drop0trailing = T,digits=0)
table[,5:8]<-round(table[,5:8],4)
names(table)[3]<-'covariates'
for (r in 1:nrow(table)){
  covar<-strsplit(table$covariates[r],'+',fixed=T)[[1]][2]
  if (is.na(covar)){
    table$covariates[r]<-'None'
    next()
  }
  table$covariates[r]<-covar
}
table$covariates<-gsub('_combined','',table$covariates)
table$covariates<-gsub('seastate','sea state',table$covariates)
openxlsx::write.xlsx(table,file.path(PUB_PACKAGE_DIR,'Table 4.xlsx'))

#### Table 5 survey results summary ####
load(MSM115_DS_RESULT)
m<-MSM115_DS_RESULTS$model

nL<-m$dht$individuals$summary #enc rates etc.
table<-sqldf::sqldf('select Region as stratum, Area as area_km2, Effort as effort_km, -1 as G,-1 as I,-1 as gs,ER as nL_groups_km, -1 as Ni, -1 as Di, -1 as CV from nL')
table$area_km2<-formatC(table$area_km2,format='f',big.mark = ',', drop0trailing = T,digits=0)

table$G<-m$dht$individuals$summary$n

abu<-m$dht$individuals$N #abundance
abu<-sqldf::sqldf('select Label as stratum, 0 as Ni, Estimate, se, lcl, ucl,cv from abu')
lcl<-formatC(abu$lcl,format='f', big.mark=',', drop0trailing = T,digits=0)
ucl<-formatC(abu$ucl,format='f', big.mark=',', drop0trailing = T,digits=0)
abu$Ni<-formatC(abu$Estimate,format='f', big.mark=',', drop0trailing = T,digits=0)
abu$Ni<-paste0(abu$Ni, ' ± ', formatC(round(abu$se,0),format='f', big.mark=',', drop0trailing = T))
abu$Ni<-paste0(abu$Ni, '\n(',lcl , ' - ',ucl,')' )
abu$cv<-round(abu$cv,4)*100
table$Ni<-abu$Ni
table$CV<-abu$cv

D<-m$dht$individuals$D #density
lcl<-formatC(D$lcl,format='f', big.mark=',', drop0trailing = T,digits=4)
ucl<-formatC(D$ucl,format='f', big.mark=',', drop0trailing = T,digits=4)
D$Di<-formatC(D$Estimate,format='f', big.mark=',', drop0trailing = F,digits=4)
D$Di<-paste0(D$Di, ' ± ', formatC(D$se,format='f', big.mark=',', drop0trailing = F,digits=4))
D$Di<-paste0(D$Di, '\n(',lcl , ' - ',ucl,')' )
table$Di<-D$Di

gs<-m$dht$Expected.S
table$gs<-round(gs$Expected.S,2)

I<-m$dht$individuals$bysample
I<-sqldf::sqldf('select Region as stratum, sum(n) as I from I group by Region')
I<-rbind(I, data.frame(stratum='Total',I=sum(I$I)))
table$I<-I$I
openxlsx::write.xlsx(table,file.path(PUB_PACKAGE_DIR,'Table 5.xlsx'))

#### Table 6 cluster vs multimodality ####
table<-openxlsx::read.xlsx(MSM115_KRILL_MULTIMODALITY_CLUSTER_NAME_FILE)
names(table)[1]<-'cluster'
names(table)[2]<-'N'

kill<-grep('.*_r_multimodal_[total|yes|no]',names(table))
table<-table[,-kill]
keep<-grep('cluster|N|_multimodal_yes|_multimodal_no',names(table))
table<-table[,keep]

out<-matrix(nrow=4,ncol=12)
out[3:4,]<-as.matrix(table)
out[2,1:12]<-c('','N',rep(c('multi','uni'),5))
out[1,1:4]<-c('','','Any species','Any species')
out[1,grep('frigida',names(table))]<-c('E. frigida','')
out[1,grep('superba',names(table))]<-c('E. superba','')
out[1,grep('triacantha',names(table))]<-c('E. triacantha','')
out[1,grep('macrura',names(table))]<-c('E. macrura','')
openxlsx::write.xlsx(out,file.path(PUB_PACKAGE_DIR,'Table 6.xlsx'), colNames=F)

#### ----- Table creation for supplement ---- ####

#### Table S1 bulk distance sampling results ####
model_data<-openxlsx::read.xlsx(MSM115_DS_MODEL_XLSX)
model_data[is.na(model_data)]<-''

for (col in c('CvM_p','CvM_W','pa_hat','pa_hat_se')){
  model_data[[col]]<-formatC(model_data[[col]],format='f',big.mark=',',drop0trailing = F, digits=4)
}

model_data$truncation<-formatC(model_data$truncation,format='f',big.mark=',',digits=0)
model_data$AIC<-formatC(model_data$AIC,format='f',big.mark=',',drop0trailing = F, digits=2)
model_data$N_sigs<-formatC(model_data$N_sigs,format='f',big.mark=',',drop0trailing = F, digits=0)

translation<-c(hn='Half normal', hr='Hazard-rate',unif='Uniform')
model_data$key<-translate(model_data$key,translation)
translation<-c(cos='cosine', herm='hermite',poly='polynomial')
model_data$used_adj<-translate(model_data$used_adj,translation)
model_data$ord<-model_data$order
model_data<-sqldf::sqldf('select modelname, key, formula, truncation, used_adj as adj, ord as "order", AIC, CvM_p, CvM_W, pa_hat,pa_hat_se,N_sigs from model_data')
openxlsx::write.xlsx(model_data,file.path(PUB_PACKAGE_DIR,'Table S1.xlsx'))

#### Table S2 station activities ####
load(MSM115_KRILL_COMPLETE_DATA)
load(OCE_DATA_FILE)

krill_stations<-KRILL_FINAL$abustationskrill_stations<-KRILL_FINAL$abu_age$station_id
ctd_stations<-unique(MSM115_OCE_DATA$unified_detrended$station_id)

station_summary<-sqldf::sqldf('select stratum, NULL as running_number, station_id as station, original_id from stations group by stratum, station_id')
station_summary$ctd<-station_summary$trawl<-''
station_summary$ctd[station_summary$station %in% ctd_stations]<-'x'
station_summary$trawl[station_summary$station %in% krill_stations]<-'x'

station_summary<-station_summary[order(station_summary$stratum,station_summary$station,station_summary$ctd,station_summary$trawl),]
kill<-which(names(station_summary) == 'original_id')
station_summary<-station_summary[,-kill]
station_summary<-station_summary[!is.na(station_summary$station),]
station_summary$running_number<-1:nrow(station_summary)
station_summary$stratum<-as.character(station_summary$stratum)
idx<-grep('.*(SSI)',station_summary$stratum)
station_summary$stratum[idx]<-'SSI'
idx<-grep('.*(EI)',station_summary$stratum)
station_summary$stratum[idx]<-'EI'
idx<-grep('.*(SOI)',station_summary$stratum)
station_summary$stratum[idx]<-'SOI'
idx<-grep('Between.*',station_summary$stratum)
station_summary$stratum[idx]<-'BT'

openxlsx::write.xlsx(station_summary,file.path(PUB_PACKAGE_DIR,'Table S2.xlsx'))

#### Table S3 GMM results ####
gmm<-openxlsx::read.xlsx(MSM115_KRILL_MULTIMODALITY_FILE)
gmm$species<-gsub('Euphausia','E.', fixed=T, gmm$species)
gmm$species<-gsub('Thysanoessa','T.', fixed=T, gmm$species)
table<-sqldf::sqldf('select stratum, station_id, species, n, BIC, df, loglik, dip_test_p as Dp, dip_sig as D_alpha, n_components, model_type, component, mixing_prop as p_mix, cumulative_prop as p_cum, mu, stdev from gmm')

idx<-grep('.*(SSI)',table$stratum)
table$stratum[idx]<-'SSI'
idx<-grep('.*(EI)',table$stratum)
table$stratum[idx]<-'EI'
idx<-grep('.*(SOI)',table$stratum)
table$stratum[idx]<-'SOI'
idx<-grep('Between.*',table$stratum)
table$stratum[idx]<-'BT'

for (col in c('BIC','loglik','Dp','p_mix','p_cum','mu','stdev')){
  table[[col]]<-formatC(table[[col]],format='f',big.mark=',',drop0trailing = F, digits=4)
}

openxlsx::write.xlsx(table,file.path(PUB_PACKAGE_DIR,'Table S3.xlsx'))
