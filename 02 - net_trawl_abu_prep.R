rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #only works in rstudio / posit
source(file.path(getwd(),'_SRC','_startup.R'))

#### ----- Net trawl abundance data preparation ---- #####
stations<-sf::st_read(MSM115_SPATIAL_DATA, layer='MSM115_stations') # Just the coordinates and id of each station
stations<-sf::st_drop_geometry(stations)

abu <- openxlsx::read.xlsx(MSM115_KRILL_ABU_FILE,sheet=1)
unique(abu$STATION)
length(unique(abu$STATION))
ABU<-sqldf::sqldf('select * from abu as a left join stations as b on a.STATION = b.station_id')

### quick reordering
station_cols<-which(names(ABU) %in% names(stations))
ABU<-cbind(ABU[,station_cols],ABU[,-station_cols])

ABU<-subset(ABU, !is.na(ABU$scientificName))
ABU<-ABU[ABU$Valid=='TRUE',]
length(unique(ABU$STATION))
ABU<-ABU[ABU$COLLECTOR==8,]
length(unique(ABU$STATION))

kill_cols<-c("AREA","TRANSECT.NO","FINWAP.STATION","STATION","Valid","COLLECTOR","COLLECTOR.NO","DEVICE.OPERATION","DEVICE","START.TIME","START.LAT.DEC",
             "START.LONG.DEC" ,"DEPTH.START","DEPTH.END","TRAWL.START","TRAWL.END","TRAWL.DISTANCE","TRAWL.SPEED",
             "OccurenceID","individualCount_Qualitycheck","EVENT.ID","COLLECTOR.ID","internal_occurence_id","occurrenceStatus",
             "organismQuantityType","organismQuantity","n_1000m3",'original_id')
ABU<-ABU[,-which(names(ABU) %in% kill_cols)]

names(ABU)<-tolower(names(ABU))
ABU$individualcount[is.na(ABU$individualcount)]<-0
unique(ABU$scientificname)
ABU$start.date<-as.Date(ABU$start.date, origin = "1899-12-30")

ABU$n_1000m3<-ABU$individualcount/ABU$filtered.volume*1000

# remove random species not suitable for this net
species<-tolower(ABU$scientificname)
ABU<-ABU[!(ABU$scientificname=="Euphausia superba furcilia"),]
kill_species=species[grep('chaetognatha.*|spongiobranchaea.*',species)]
ABU<-ABU[!(ABU$scientificname %in% kill_species),]
unique(ABU$scientificname)

# renaming
spp_idx<-which(ABU$scientificname %in% c('Electrona','Beroe','Clio','Cyllopus','Notolepis','Tomopteris','Vibilia','Gymnoscopelus'))
ABU$scientificname[spp_idx]<-paste0(ABU$scientificname[spp_idx],' spp.')
ABU$scientificname[ABU$scientificname=='Fish larvae Unid_']<-'Unidentified fish larvae'
ABU$scientificname[ABU$scientificname=='Fish Unid_']<-'Unidentified fish'
unique(ABU$scientificname)

head(ABU)

# move species into own columns
out<-data.frame(station_id=unique(ABU$station_id))
out$x<-NA
out$y<-NA
out$filtered_volume<-NA
out[,unique(ABU$scientificname)] <- NA

for (st in unique(ABU$station_id)){
  df<-ABU[ABU$station_id==st,]
  species<-sqldf::sqldf('select scientificname, x,y, sum(individualcount) as individualcount, n_1000m3 from df group by scientificname')
  counts<-species$individualcount
  names(counts)<-species$scientificname
  
  out[out$station_id==st,names(counts)]<-counts
  out[out$station_id==st,c('x','y')]<-c(species$x[1],species$y[1])
  out[out$station_id==st,'filtered_volume']<-sum(unique(df$filtered.volume),na.rm=T)
  
}
out[,4:ncol(out)][is.na(out[,4:ncol(out)])]<-0

#### add diversity indices
for (i in 1:nrow(out)){
  counts<-out[i,unique(ABU$scientificname)]
  out$H_shannon[i]<-vegan::diversity(counts,'shannon')
  out$H_simpson[i]<-vegan::diversity(counts,'simpson')
}

species<-tolower(colnames(out))
taxa<-list(
  amphipods=species[grep('amphipoda.*|cyllopus.*|hyperia .*|hyperiidea.*|themisto.*|vibilia.*',species)],
  pteropods=species[grep('clio.*|clione.*',species)],
  hydrozoa=species[grep('calycopsis.*|beroe.*|diphyes.*|hydromedusae.*|periphylla.*|scyphomedusae.*|spongiobranchaea.*',species)],
  tunicates=species[grep('salpa.*|ihlea.*',species)],
  other=species[grep('chaetognatha.*|phytoplankton.*|tomopteris.*|vanadis.*|zooplankton.*|acanthephyra.*|mysidacea.*|galiteuthis.*',species)],
  fish=species[grep('electrona.*|chaenodraco.*|channichthyid.*|.*fish.*|gymnoscopelus.*|notolepis.*',species)],
  krill=species[grep('euphausia.*|thysanoessa.*',species)]
)

for (i in 1:length(taxa)){
  idx<-which(species %in% taxa[[i]])
  for (r in 1:nrow(out)){
    out[[paste0('N_grp_',names(taxa)[i])]][r]<-sum(out[r,idx])
    out[[paste0('H_shannon_',names(taxa)[i])]][r]<-vegan::diversity(out[r,idx],'shannon')
    out[[paste0('H_simpson_',names(taxa)[i])]][r]<-vegan::diversity(out[r,idx],'simpson')
  }
}

NET_TRAWLS<-out

NET_TRAWLS_sf<-sf::st_as_sf(NET_TRAWLS, coords = c('x','y'),crs=IBCSO,remove=F)
ABU_sf<-sf::st_as_sf(ABU, coords = c('x','y'),crs=IBCSO,remove=F)

sf::st_write(NET_TRAWLS_sf,dsn=MSM115_KRILL_ABU_GPKG,layer='MSM115_net_trawls_horizontal_format', append=F, delete_dsn=T)
sf::st_write(ABU_sf,dsn=MSM115_KRILL_ABU_GPKG,layer='MSM115_net_trawls_vertical_format', append=T)

names(NET_TRAWLS)
station_columns<-list()
cols <- c('station_id','x','y','filtered_volume','H_shannon','H_simpson')
cols<-which(names(NET_TRAWLS) %in% cols)
for (i in cols){
  station_columns[names(NET_TRAWLS)[i]]<-i  
}

species_columns<-list()
for (i in c(5:39)){
  species_columns[names(NET_TRAWLS)[i]]<-i  
}
group_columns<-list()
for (i in c(42:62)){
  group_columns[names(NET_TRAWLS)[i]]<-i  
}
column_list<-list(station_columns=station_columns,species_columns=species_columns,group_columns=group_columns)

ABU$tax_group<-NA
species<-tolower(ABU$scientificname)

for (i in 1:length(names(taxa))){
  ABU$tax_group[species %in% taxa[[i]]]<-names(taxa)[[i]]
}
ABU_summary<-sqldf::sqldf('select sum(individualCount) as caught_individuals,"" as length_measurements,scientificname, tax_group from ABU group by scientificname,tax_group')

grouping<-table(ABU$tax_group,ABU$scientificname)

ggHeatmap(grouping)

tabulated<-sqldf::sqldf('select count(*) as N_stations, sum(individualcount) as N_ind, scientificname as species, tax_group from ABU group by species, tax_group')
openxlsx::write.xlsx(tabulated,MSM115_KRILL_ABU_SUMMARY_FILE)

KRILL<-list(MSM115_krill_final=NET_TRAWLS,
            MSM115_krill_abundance_summary=ABU_summary,
            MSM115_krill_abundance=ABU, 
            MSM115_krill_final_columns=column_list)
save(KRILL, file=MSM115_KRILL_ABU_DATA, compress='gzip')

### plotting ####
STRATA<-sf::st_read(MSM115_STRATA_FILE,layer='simplified_strata')
ABU$stratum<-sf::st_join(ABU_sf,STRATA)$ShortName
ABU$stratum<-factor(ABU$stratum,levels=c('SSI','EI','SOI','Between strata'))
ABU$stratum<-sf::st_join(ABU_sf,STRATA)$name
ABU$stratum[is.na(ABU$stratum)]<-'Between strata'
ABU$stratum<-factor(ABU$stratum,levels=c('South Shetland Islands (SSI)','Elephant Island (EI)','South Orkney Island (SOI)','Between strata'))

d<-subset(ABU,tax_group %in% c('krill','tunicates'))
d$tax_group<-toupper(d$tax_group)
TOTALS<-sqldf::sqldf('select stratum, station_id as stid,sum(individualcount) as N_total from d group by stratum,station_id')
N<-sqldf::sqldf('select * from d as D left join TOTALS as T on D.station_id=T.stid')
N<-N[,-which(names(N) == 'stid')]
N$ratio<-N$individualcount/N$N_total
N$species<-factor(N$scientificname,levels=c('Euphausia superba','Euphausia triacantha','Euphausia frigida','Thysanoessa macrura','Salpa thompsoni'))

N2<-subset(N,tax_group!='TUNICATES')
N2$type<-'Without salps'
N$type<-'All species'
N<-rbind(N,N2)

p<-ggplot(N, aes(fill=species, y=ratio, x=factor(station_id))) + geom_bar(position="fill", stat="identity")
p<-p+labs(x='')+scale_fill_manual(values=c(krill_colors,tunicate_colors))+labs(fill='')
p<-p+MAINTHEME+theme(legend.position = "right")

png(file.path(PUB_PACKAGE_DIR,'Figure 6.png'),7000,3000,res=300)
p+facet_grid(type~stratum, scales='free_x') + theme_large_font()
graphics.off()
