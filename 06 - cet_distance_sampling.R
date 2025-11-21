rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #only works in rstudio / posit
source(file.path(getwd(),'_SRC','_startup.R'))

#### ----- Bulk distance sampling modelling ---- #####
require(Distance)

#### In case this has been run previously, you can simply go directly to line 66

GFX_DS<-file.path(GFX_DIR,'DS')
dir.create(GFX_DS, recursive = T, showWarnings = F)

load(file=MSM115_CET_ANALYSIS_DATA)

ds_data<-MSMS115_CET$distance_data$ds_data
region_table<-MSMS115_CET$distance_data$region_table
sample_table<-MSMS115_CET$distance_data$sample_table
obs_table<-MSMS115_CET$distance_data$obs_table

ds_summary<-sqldf::sqldf('select "Region.Label" as stratum, sum(Effort) as effort_km from sample_table group by "Region.Label"')
sig_summary<-sqldf::sqldf('select * from obs_table as a left join ds_data as b on a.object = b.object')
sig_summary<-sig_summary[,-duplicated(names(sig_summary))]
sig_summary<-sqldf::sqldf('select "Region.Label" as stratum, count(*) as G, sum(size) as I from sig_summary group by "Region.Label"')

ds_summary<-sqldf::sqldf('select * from region_table as a left join ds_summary as b on a."Region.Label" = b.stratum')
ds_summary<-sqldf::sqldf('select stratum, Area as area_km2, effort_km, G, I from (select * from ds_summary as a left join sig_summary as b on a.stratum = b.stratum)')
ds_summary$nL<-ds_summary$G/ds_summary$effort_km
ds_summary$gs<-ds_summary$I/ds_summary$G
openxlsx::write.xlsx(ds_summary, file=MSM115_CET_DS_SUMMARY_FILE)

# histogram inspection
hist(ds_data$distance, breaks=seq(0, 6000,250))
TRUNC<-c(1500,1750,seq(3000,4000,250))
COV<-c(NA,'seastate','swell','subj','glare')

check_covar_levels(ds_data,TRUNC,COV)
ds_data$subj_combined<-as.character(ds_data$subj)
ds_data$subj_combined[ds_data$subj%in%c('m','p')]<-'mp'
ds_data$glare_combined<-as.character(ds_data$glare)
ds_data$glare_combined[ds_data$glare%in%c(1,2)]<-'12'
ds_data$seastate_combined<-as.character(ds_data$seastate)
ds_data$seastate_combined[ds_data$seastate%in%c(4,5)]<-'45'

COV<-c(NA,'seastate_combined','glare_combined','subj_combined')
check_covar_levels(ds_data,TRUNC,COV)

cov.mat<-create_ds_settings(truncations=TRUNC,adjustments=c(NA,'cos','herm','poly'), keys=c('hn','hr','unif'),covars=COV)

#### bulk distance sampling
#### when adjustments are chosen, the model chosen within the adjustment series is based on AIC.
#### Some of these adjusted models may be not strictly monotonic
ds_model_summary<-bulk_ds(ds_data=ds_data,cov.mat=cov.mat)

model_data<-ds_model_summary[ds_model_summary$N_sigs>=20,]
model_data<-model_data[model_data$pa_hat_se<=1,]
model_data<-model_data[model_data$CvM_p>0.05,]
model_data<-model_data[model_data$pa_hat>=0.4,]

identical<-paste(model_data$key,model_data$formula,model_data$truncation,model_data$used_adj, model_data$order,'|')
model_data<-model_data[!duplicated(identical),]
model_data
openxlsx::write.xlsx(model_data, file=MSM115_DS_MODEL_XLSX)
openxlsx::write.xlsx(cov.mat, file=MSM115_DS_MODEL_COV_XLSX)

#### Start here in case you've run the long bulkd ds processing step
model_data<-openxlsx::read.xlsx(MSM115_DS_MODEL_XLSX)

candidate_models<-data.frame()
for (tr in unique(model_data$truncation)){
  x<-model_data[model_data$truncation==tr,]
  x$dAIC<-x$AIC-min(x$AIC)
  x<-x[order(x$dAIC),]
  x<-x[x$dAIC<=10,]
  x<-x[order(-x$CvM_p),]
  row.names(x)<-1:nrow(x)
  candidate_models<-rbind(candidate_models,x)

  for (i in 1:nrow(x)){
    cov<-x[i,]
    modelname<-cov$modelname
    adj<-cov$used_adj
    
    if (is.na(adj)){
      adj<-'cos'
      nadj<-NULL
    }else{
      nadj<-strsplit(gsub('\\(|\\)','',cov$order),',')[[1]]
      nadj<-as.numeric(nadj)
    }
    print(modelname)
    
    m<-Distance::ds(data=ds_data, region_table = region_table, sample_table=sample_table, obs_table=obs_table, 
                    key=cov$key, formula=as.formula(cov$formula), truncation=cov$truncation, adjustment=adj,nadj=nadj,
                    convert_units=0.001, quiet=T)
    
    filename=file.path(GFX_DS,paste0('MSM115_detfct_',modelname,'.png'))
    N_hat<-m$dht$individuals$N
    N_hat[,2:3]<-format(round(N_hat[,2:3],0),big.mark=',')
    lgd_txt<-c('Abundances')
    pad=max(nchar(N_hat[,1]))+2
    
    for (row in 1:nrow(N_hat)){
      l<-nchar(N_hat[row,1])
      lgd_txt<-c(lgd_txt,paste0(N_hat[row,1],':', strrep(' ',pad-l),N_hat[row,2],' +- ',N_hat[row,3]))
    }
    
    p<-gg_DS_DetFct(m,size.esw=4,thousands_sep=',')+MAINTHEME+theme_large_font()
    g<-gg_DS_GOF(m)+MAINTHEME+theme_large_font()
    g<-g+geom_label(aes(x = 0, y = 1,label = paste0(lgd_txt,collapse='\n')), vjust=1,hjust=0,fill='#d1e2f1')

    png(filename, 2400,2000,res=150)
    gridExtra::grid.arrange(grobs=list(p,g))
    graphics.off()
    
    if (!cov$is_mono){
      filename=file.path(GFX_DS,paste0('MSM115_detfct_monotonicity_',modelname,'.png'))
      png(filename, 1200,2000,res=150)
      check.mono(m$ddf, plot=TRUE, n.pts=100)
      graphics.off()
    }
  }
}
openxlsx::write.xlsx(candidate_models,MSM115_CET_DS_CANDIDATE_MODELS)

final_model<-Distance::ds(data=ds_data, region_table = region_table, sample_table=sample_table, obs_table=obs_table, 
                          key='hn', formula=~1, truncation=1750, adjustment=NULL,nadj=NULL,
                          convert_units=0.001, quiet=T)

MSM115_DS_RESULTS<-list(model=final_model,ds_data=ds_data, region_table=region_table, sample_table=sample_table, obs_table=obs_table,abundance=final_model$dht)
save(MSM115_DS_RESULTS,file=MSM115_DS_RESULT,compress='gzip')

openxlsx::write.xlsx(final_model$dht$individuals$N,MSM115_CET_DS_RESULTS)

#### plotting ####
load(file=MSM115_DS_RESULT)

m<-MSM115_DS_RESULTS$model
N_hat<-m$dht$individuals$N
N_hat[,2:3]<-format(round(N_hat[,2:3],0),big.mark=',')
lgd_txt<-c('Abundances')
pad=max(nchar(N_hat[,1]))+2

for (row in 1:nrow(N_hat)){
  l<-nchar(N_hat[row,1])
  lgd_txt<-c(lgd_txt,paste0(N_hat[row,1],':', strrep(' ',pad-l),N_hat[row,2],' +- ',N_hat[row,3]))
}

p<-gg_DS_DetFct(m, thousands_sep = ',',n_hist_bins=6,auto_adjust_breaks = F,size.esw=5)+MAINTHEME+theme_large_font()
g<-gg_DS_GOF(m,add_lines=T)+MAINTHEME+theme_large_font()
g<-g+geom_label(aes(x = 0, y = 1,label = paste0(lgd_txt,collapse='\n')), vjust=1,hjust=0,fill='#d1e2f1')

png(file.path(GFX_DS,'MSM115_final_det_fct.png'), 2400,2000,res=150)
gridExtra::grid.arrange(grobs=list(p,g))
graphics.off()

png(file.path(PUB_PACKAGE_DIR,'Figure 3.png'),3200,3200,res=300)
print(p)
graphics.off()
