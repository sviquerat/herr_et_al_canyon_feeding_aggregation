rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #only works in rstudio / posit
source(file.path(getwd(),'_SRC','_startup.R'))

#### ----- Distance sampling data prep ---- #####
ds_survey<-read.csv(MSM115_CET_SURVEY_DATA,dec='.',sep='\t')
LL<-sf::st_as_sf(ds_survey,coords=c('lon','lat'),crs=WGS84)
LL<-sf::st_transform(LL,IBCSO)
cc<-sf::st_coordinates(LL)
ds_survey$x<-sf::st_coordinates(LL)[,1]
ds_survey$y<-sf::st_coordinates(LL)[,2]

# all sightings recorded on effort - except aggregation
survey<-ds_survey[!is.na(ds_survey$species),]
survey$orig_id<-survey$MSM115_id
survey$data_source<-'vis_survey'
survey<-subset(survey,species %in% c('bphy','aggr'))
survey<-data.frame(id=0, orig_id = survey$orig_id, date=survey$date, datetime=survey$datetime, tz = survey$timezone, species=survey$species, best=survey$best_number, calves=survey$calves, com=survey$com_sight,lon=survey$lon,lat=survey$lat,x=survey$x, y= survey$y,data_source=survey$data_source)

# There was one aggregation on effort - take care - they'll be added back to table 3 for documentation purposes. It does not matter for the analysis where they appear from
survey_agg<-sqldf::sqldf('select orig_id as id, date, datetime, tz, species, best, calves, com, lon,lat,x,y from survey where species=="aggr"')
survey<-subset(survey, species!='aggr')

opp_sightings<-read.csv(MSM115_CET_OPP_DATA,sep='\t',dec='.')
LL<-sf::st_as_sf(opp_sightings,coords=c('lon','lat'),crs=WGS84)
LL<-sf::st_transform(LL,IBCSO)
opp_sightings$x<-sf::st_coordinates(LL)[,1]
opp_sightings$y<-sf::st_coordinates(LL)[,2]

opp_sightings<-rbind(opp_sightings,survey_agg)
opp_sightings<-opp_sightings[order(opp_sightings$datetime),]
opp_sightings$orig_id<-opp_sightings$id
opp_sightings$id<-1:nrow(opp_sightings)

sigs<-subset(opp_sightings, species!='aggr')
sigs$data_source<-'opp_sightings'
agg<-subset(opp_sightings, species=='aggr')
agg$data_source<-'aggregations'

sightings<-rbind(agg, survey, sigs)
sightings<-sightings[order(as.Date(sightings$date)),]
sightings$id<-1:nrow(sightings)

sightings$date<-gsub('-3-','-03-', fixed=T,sightings$date)
sightings$date<-sub("-([0-9])$", "-0\\1", sightings$date)

sightings_table<-sqldf::sqldf('select data_source, species, count(*) as G, sum(best) as I from sightings group by data_source, species')
openxlsx::write.xlsx(sightings_table,file=MSM115_CET_SIGHTINGS)

SIGS<-sf::st_as_sf(sightings, coords=c('x','y'),crs=IBCSO, remove=F)

#---------------------Spatial data--------------------#
SIGS<-subset(SIGS,species %in% c('bphy','aggr'))

sigs<-SIGS

sigs_data<-sf::st_drop_geometry(sigs)

#### distance data set
STRATA<-sf::st_read(MSM115_STRATA_FILE,layer='simplified_strata')
STRATA$area_km2<-as.numeric(units::set_units(sf::st_area(STRATA),km^2))

surveyData<-ds_survey

data<-sqldf::sqldf('select MSM115_id as object, lat, lon, x,y,transect, track_length_km as effort_km, seastate, swell, glare, sub_left, sub_right, cond, side, perp_dist_m as distance, best_number as size, species, "" as stratum from surveyData')

data$subj<-NA
data$subj[data$side=='L']<-data$sub_left[data$side=='L']
data$subj[data$side=='R']<-data$sub_left[data$side=='R']

LL<-sf::st_as_sf(data, coords=c('x','y'), crs=IBCSO)
strata<-sf::st_intersects(LL,STRATA, sparse=T)
data$stratum<-as.numeric(strata)
data$RegionLabel<-STRATA$name[data$stratum]
data$SampleLabel<-data$transect
data$Area<-STRATA$area_km2[data$stratum]

region_table<-sqldf::sqldf('select RegionLabel as "Region.Label" , Area from data group by RegionLabel, Area')
sample_table<-sqldf::sqldf('select RegionLabel as "Region.Label", transect as "Sample.Label", sum(effort_km) as Effort from data group by stratum, transect')
obs_table<-sqldf::sqldf('select object, RegionLabel as "Region.Label", transect as "Sample.Label" from data where species = "bphy"')
ds_data<-sqldf::sqldf('select object, distance, size,  seastate, swell, glare, subj from data where species ="bphy"')
ds_data$seastate<-as.factor(ds_data$seastate)
ds_data$subj<-as.factor(ds_data$subj)
ds_data$swell<-as.factor(ds_data$swell)
ds_data$glare<-as.factor(ds_data$glare)

MSMS115_CET=list(distance_data = list(ds_data=ds_data, region_table = region_table, sample_table=sample_table, obs_table=obs_table),
                 sig_data=sigs_data)
save(MSMS115_CET, file=MSM115_CET_ANALYSIS_DATA, compress='gzip')
