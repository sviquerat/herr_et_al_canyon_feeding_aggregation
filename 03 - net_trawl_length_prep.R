rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #only works in rstudio / posit
source(file.path(getwd(),'_SRC','_startup.R'))

#### ----- Net trawl length data preparation ---- #####
stations<-sf::st_read(MSM115_SPATIAL_DATA, layer='MSM115_stations') # Just the coordinates and id of each station
stations<-sf::st_drop_geometry(stations)

lf<-openxlsx::read.xlsx(MSM115_KRILL_LENGTH_FILE,sheet=1)
LF<-sqldf::sqldf('select * from lf as a left join stations as b on a.STATION = b.station_id')

### quick reordering
station_cols<-which(names(LF) %in% names(stations))
LF<-cbind(LF[,station_cols],LF[,-station_cols])

LF<-subset(LF, !is.na(LF$scientificName))
LF<-LF[LF$Valid=='TRUE',]
LF<-LF[LF$COLLECTOR==8,] # is this relevant?
names(LF)
kill_cols<-c("AREA","TRANSECT.NO","FINWAP.STATION","STATION","Valid","COLLECTOR","COLLECTOR.NO","DEVICE.OPERATION","DEVICE","START.TIME","START.LAT.DEC",
             "START.LONG.DEC" ,"DEPTH.START","DEPTH.END","TRAWL.START","TRAWL.END","TRAWL.DISTANCE","TRAWL.SPEED",
             "OccurenceID","individualCount_Qualitycheck","EVENT.ID","COLLECTOR.ID","internal_occurence_id","occurrenceStatus",
             "organismQuantityType","organismQuantity","n_1000m3",'MAX.FISHING.DEPTH', 'DESIRED_DEPTH', 'FILTERED.VOLUME','original_id')
LF<-LF[,-which(names(LF) %in% kill_cols)]

names(LF)<-tolower(names(LF))
LF$start.date<-as.Date(LF$start.date, origin = "1899-12-30")

# remove random species not suitable for this net
species<-tolower(LF$scientificname)
LF$species<-LF$scientificname
LF$age.from.length<-LF$old_

LF<-LF[,-which(names(LF) %in% c('scientificname','old_'))]
unique(LF$species)

LF$stage<-factor(LF$stage,levels=rev(sort(unique(LF$stage))),ordered=T)
LF$stage2<-gsub('2 A.*','2A',LF$stage)
LF$stage2<-gsub('2 B.*','2B',LF$stage2)
LF$stage2<-gsub('3 A.*','3A',LF$stage2)
LF$stage2<-gsub('3 B.*','3B',LF$stage2)
LF$stage2<-gsub('3 C.*','3C',LF$stage2)
LF$stage2<-gsub('3 D.*','3D',LF$stage2)
LF$stage2<-gsub('3 E.*','3E',LF$stage2)
LF$stage2<-factor(LF$stage2,levels=rev(sort(unique(LF$stage2))),ordered=T)

LF$stage3<-as.character(LF$stage)
LF$stage3[substr(LF$stage,1,1)=='2']<-'2'
LF$stage3[substr(LF$stage,1,1)=='3']<-'3'
LF$stage3<-factor(LF$stage3,levels=rev(sort(unique(LF$stage3))),ordered=T)
LF$stage3

LF$age.group[is.na(LF$age.group)]<-'UNDEF'
LF$age.group<-factor(LF$age.group,levels=c('JUVENILE','SUBADULT','ADULT','UNDEF'),ordered=T)

LF$sex[is.na(LF$sex)]<-'UNDEF'
LF$sex<-factor(LF$sex,levels=c('JUV','M','F','UNDEF'),ordered=T)

### extraction
load(MSM115_KRILL_ABU_DATA)
abundance<-KRILL$MSM115_krill_final

SOI<-unique(LF$species)[grep('.*superba|.*frigida|.*triacantha|thysanoessa.*',tolower(unique(LF$species)))]
combined_abundance<-data.frame()
for (id in unique(LF$station_id)){
  df<-subset(LF,station_id==id&length>=0)
  abu<-subset(abundance,station_id==id)
  ratios_table<-data.frame()
  for (sp in SOI){
    x<-subset(df,species==sp)
    div<-0
    if (nrow(x) >0){
      div<-1/nrow(x)
    }
    ratios<-table(x$age.group)*div
    
    abu_col<-which(names(abu) %in% sp)
    
    new_data<-data.frame(t(as.matrix(ratios*abu[,abu_col])))
    names(new_data)<-paste0(sp,'_',names(new_data))
    new_data[[paste0(sp,'_','juv_ratio')]]<-ratios[['JUVENILE']]
    q<-quantile(x$length)
    names(q)<-paste0('length_',sp,'_',c('q0','q25','q50','q75','q100'))
    new_data<-cbind(new_data,as.data.frame(t(q)))
    mu<-mean(x$length,na.rm=T)
    sigma_sq<-var(x$length,na.rm=T)
    std<-sd(x$length,na.rm=T)
    new_data[[paste0('length_',sp,'_mu')]]<-mu
    new_data[[paste0('length_',sp,'_sigma')]]<-sigma_sq
    new_data[[paste0('length_',sp,'_std')]]<-std
    new_data[[paste0('length_',sp,'_N')]]<-nrow(x[!is.na(x$length),])
    if (nrow(ratios_table)==0){
      ratios_table<-data.frame(station_id=id,new_data)
    }else{
      ratios_table<-cbind(ratios_table,new_data) 
    }
  }
  combined_abundance<-rbind(combined_abundance,ratios_table)
}

juv_cols<-grep('.*JUVENILE',names(combined_abundance))
abu_cols<-grep('.*_JUVENILE|ADULT|SUBADULT|UNDEF',names(combined_abundance))
l50_cols<-grep('.*q50',names(combined_abundance))
l25_cols<-grep('.*q25',names(combined_abundance))
l75_cols<-grep('.*q75',names(combined_abundance))

combined_abundance$krill_N<-apply(combined_abundance[,abu_cols],FUN=sum,MARGIN=1)
combined_abundance$krill_N_juv<-apply(combined_abundance[,juv_cols],FUN=sum,MARGIN=1)
combined_abundance$krill_ratio_juv<-combined_abundance$krill_N_juv/combined_abundance$krill_N
combined_abundance$krill_ratio_juv[is.na(combined_abundance$krill_r_juv)]<-0
combined_abundance$krill_length_25<-apply(combined_abundance[,l25_cols],FUN=median,MARGIN=1,na.rm=T)
combined_abundance$krill_length_50<-apply(combined_abundance[,l50_cols],FUN=median,MARGIN=1,na.rm=T)
combined_abundance$krill_length_75<-apply(combined_abundance[,l75_cols],FUN=median,MARGIN=1,na.rm=T)
combined_abundance$krill_length_IQR<-combined_abundance$krill_length_75-combined_abundance$krill_length_25

juv_abu<-sqldf::sqldf('select * from abundance as a left join combined_abundance as b on a.station_id = b.station_id')
kill<-max(which(names(juv_abu) == 'station_id'))
juv_abu<-juv_abu[,-kill]
juv_abu$krill_D_juv<-juv_abu$krill_N_juv/juv_abu$filtered_volume*1000
juv_abu$krill_D<-juv_abu$krill_N/juv_abu$filtered_volume*1000

juv<-juv_abu

### station - strata table
#this should go into a function
STRATA<-sf::st_read(MSM115_STRATA_FILE,layer='simplified_strata')
LL<-sf::st_as_sf(juv,coords=c('x','y'),crs=IBCSO, remove=F)
juv$stratum<-sf::st_join(LL,STRATA)$name
juv$stratum[is.na(juv$stratum)]<-'Between strata'
juv$stratum<-factor(juv$stratum,levels=c('South Shetland Islands (SSI)','Elephant Island (EI)','South Orkney Island (SOI)','Between strata'))
station_strata_table<-juv[,c('station_id','stratum','x','y')]

KRILL_FINAL<-list(abu_age=juv, lf=LF,SOI=SOI, station_strata_table=station_strata_table)
save(KRILL_FINAL,file=MSM115_KRILL_COMPLETE_DATA,compress='gzip')
