rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #only works in rstudio / posit
source(file.path(getwd(),'_SRC','_startup.R'))

#### ----- SOAP parameter creation ---- #####
require(ggplot2)

sa_islands<-sf::st_read(MSM115_SPATIAL_DATA,layer='MSM115_islands') #only the islands
sa_islands <- sf::st_simplify(sa_islands, preserveTopology = FALSE, dTolerance = 5000) #simplify, otherwise it will be too complex
sa_islands = sf::st_cast(sa_islands,"POLYGON") #split up the multipolygon into individual polygons per island
sa_islands$area_km2<-as.numeric(sf::st_area(sa_islands))/1000^2 #recalculate the area of each poly, is in m², we care more about km²
sa_islands$id<-1:nrow(sa_islands) # create an id for each island
sa_islands<-subset(sa_islands,area_km2 >= 100) # remove all polys (islands) that are < 100km² (10x10km is very small for boundaries)

# next - the survey area
sa<-sf::st_read(MSM115_SPATIAL_DATA,layer='MSM115_survey_area') #islands still included - this will be the outside boundary
sa_coords<-data.frame(sf::st_coordinates(sa)) # get coordinates as data.frame - x column is X (capital x), y column is Y (capital y)

# this is the required format for boundaries - a list of lists for each boundary with x, y and f (the boundary condition, usually 0)
# the outer (survey area boundary)
bound<-list(
  list(x = sa_coords$X, y = sa_coords$Y, f = rep(0, nrow(sa_coords)))
)

# now we'll add the islands as new boundaries
# we'll use seq_along, which simply counts from 1 to the end for each loop
for (i in seq_along(unique(sa_islands$id))){ # go through each island poly
  island<-unique(sa_islands$id)[i] # get that island id
  poly<-subset(sa_islands, id==island) # get that poly out of all the islands
  coords<-data.frame(sf::st_coordinates(poly)) # get the coordinates, see above
  new_list<-list(x = coords$X, y = coords$Y, f = rep(0, nrow(coords))) # create the new list for the boundary, see above
  bound[[i+1]]<-new_list # add the new boundary to the list of previous boundaries (that's where seq_along comes in handy)
}

#### create knots ####
# there's a big discussion or nuisance in soap smoothing:
# knots are required to fall inside boundaries, but the way you construct them regularly and check this condition
# is different from the check that the actual gam model runs and often leads to errors saying knots are not within boundaries.
# that's why we'll cont´struct the knots from a spatial sampling process using buffered (shrinked and extended polygons) that include
# a safety buffer to construct our knot coordinates...

N <- 6 # number of knots (i.e. the max available number of knots, will be squared) N cannot be >= length(segments)!

safety_margin <- 10000 #our safety buffer
sa_shrink<-sf::st_buffer(sa,-safety_margin) # first, shrink the survey area by our safety margin
sai_shrink<-sf::st_buffer(sa_islands,safety_margin) # then extend the islands by the same safety margin

# sample N squared points in a regular pattern within the shrunk survey area 
GP<-sf::st_as_sf(sf::st_sample(sa_shrink, size=N**2, type='regular')) 
GP$id<-1:nrow(GP) # assign an id to each point
excl<-sf::st_filter(GP,sai_shrink) # check which points fall into the extended islands polygon
KNOTS<-GP # assume none were discarded
if (length(excl) > 0){ # but if some were discarded, do this instead
  KNOTS<-subset(GP, !(id %in% excl$id)) # remove all points based on their id
}

# summary plot - red points go into the model as knots, blue points (if any) are discarded
p<-ggplot()+geom_sf(data=sa)
p<-p+geom_sf(data=sa_shrink,col='red')
p<-p+geom_sf(data=sa_islands)
p<-p+geom_sf(data=sai_shrink,col='red')
p<-p+geom_sf(data=GP,col='blue')
p<-p+geom_sf(data=KNOTS,col='red')
p<-p+labs(title='MSM115 SOAP parameter setup')
p+MAINTHEME

knots<-data.frame(sf::st_coordinates(KNOTS)) # as before, get the coordinates of the points
knots<-data.frame(x=knots$X,y=knots$Y) # make it consistent with what gam wants

MSM115_SOAP_PARAMETERS<-list(knots=knots, boundary=bound)
save(MSM115_SOAP_PARAMETERS,file=SOAP_PARAMETERS, compress='gzip')

#### dump modified spatial data into geopackage ####
sf::st_write(sa_islands[,-which(names(sa_islands)=='ID')],dsn=SOAP_SPATIAL_DATA,layer='islands_simplified',driver='GPKG',append=F, delete_dsn=T)
sf::st_write(sf::st_read(MSM115_SPATIAL_DATA,layer='MSM115_islands'),dsn=SOAP_SPATIAL_DATA,layer='islands_original',append=T)
sf::st_write(sai_shrink[,-which(names(sai_shrink)=='ID')],dsn=SOAP_SPATIAL_DATA,layer='islands_simplified_buffered', append=T)
sf::st_write(sa_shrink,dsn=SOAP_SPATIAL_DATA,layer='survey_area_simplified_shrunk',append=T)
sf::st_write(sa,dsn=SOAP_SPATIAL_DATA,layer='survey_area_simplified', append=T)
sf::st_write(sf::st_read(MSM115_SPATIAL_DATA,layer='MSM115_survey_area'),dsn=SOAP_SPATIAL_DATA,layer='survey_area_original',append=T)
