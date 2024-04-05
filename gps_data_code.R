##### DATA ASSEMBLY AND CLEANING CODE FOR THE GPS PAPER 
### Versioned by MTJ on 9th Jan 2023.

### Necessary packages and working directory.
setwd('~/Desktop/rock doves/')

library(sf)
library(tidyverse)
library(ggmap)
library(lubridate)
library(raster)
library(tmap)


############ PATCH TOGETHER THE RAW DATA ====

pigeon_names <- data.frame(id = c("4039","5ff1","62d0","7e51","ac8e"),
                           name = c('EA49571', 'EA49502', 'EA49570', 'EA49503', 'EA49568'))

gps_data <- list.files('Raw_GPS_Data/', full.names = T)

read_csv(gps_data) %>% mutate(id = str_sub(UUID, start = -4)) %>% left_join(pigeon_names, by = 'id') %>%
  write_csv('gps_folder/gps_dataset.csv')


############ DATA FILTERING AND ORGANISING THE DATASET ====

## Main filtering stage to remove dodgy observations; keeps 4,043 out of 5,937 observations.
gps_dataset <- read_csv('gps_folder/gps_dataset.csv') %>% # load in the file created above.
  filter(hdop < 4) %>% # horizontal precision; as per nightjar paper; removes 68.
  filter(vdop < 4) %>% # same but vertical - think altitude? not used in nightjar; removes 251.
  filter(satused > 2) %>% # number of satellites used to pinpoint, as per nightjar paper; removes 678.
  st_as_sf(coords = c('lon', 'lat'), crs = st_crs(4326)) %>% # make a geo object; set standard CRS for now.
  dplyr::select(!horaccuracy) %>% dplyr::select(!veraccuracy) %>% dplyr::select(!posmode) %>% # not useful as far as I can tell.
  st_crop(xmin = -7.569580, ymin = 57.485571,xmax = -6.955719, ymax = 57.754007) %>% # remove travel through Uk with active tags.
  mutate(datetime = paste0(date, " ", time), # recreate a date-time object.
         datetime = lubridate::ymd_hms(datetime)) %>% # make it date-time
  dplyr::select(-date, -time) %>%
  arrange(datetime) %>% # order the dataset in chronological order as it may matter for calculations.
  filter(datetime <= as.Date('2023-08-08')) %>% # remove empty observations of one fallen tag.
  filter(datetime >= as.Date('2023-07-10')) # remove capture evening.

gps_dataset$datetime <- force_tz(gps_dataset$datetime, tzone = 'Europe/Warsaw')
gps_dataset$datetime <- with_tz(gps_dataset$datetime, tzone = 'Europe/London')
write_sf(gps_dataset, 'bin/filtered_gps_dataset.gpkg')


############ GENERATING MAP DATA ====

### Create the clachan base and habitat shape files to be used in the analysis.
clachan <- st_read('~/Desktop/GBR_adm/GBR_adm0.shp') %>% # from ESRI free database.
  st_transform(crs = 27700) %>% # transform to British National Grid
  st_crop(xmin = 85206.7663, ymin = 875224.6196, xmax = 89558.8570, ymax = 878304.2766) %>% # limit to Clachan.
  st_write('bin/clachan_base.gpkg') 
clachan_hab <- st_read('~/Downloads/HABMOS_SCOTLAND_SHP_27700/HABMOS_SCOTLAND.shp') %>% # from Nature Scot
  st_crop(xmin = 85206.7663, ymin = 875224.6196, xmax = 89558.8570, ymax = 878304.2766) %>% # limit to Clachan
  group_by(HABITAT_NA) %>% summarise(geometry = st_union(geometry)) %>% # unionised different habitats.
  st_write('bin/clachan_habitat.gpkg')

### Obtain satellite base map (NB: this uses horrible ggmap with no proper crsing; difficult to work with).
library(ggmap)
register_google(key = 'AIzaSyA0XK5mlbxJpGsM0Bj55CRq-emtNKAESpw') # needs a paid-for API key
lon <- c(-7.280598, -7.207275) # define longitude of the region
lat <- c(57.65647, 57.67809) # define latitude of the region
points <- st_coordinates(gps) # make them coordinates
map <- get_map(location = c(lon = mean(lon), lat = mean(lat)), zoom = 13, # obtain satellite image.
               maptype = "satellite", source = "google")
save(map, file = 'sat_map_clachan.RData') # save for re-use; NB: crs is 4326 or some other; not British Grid.

############ IDENTIFYING ROOST SITES ====

## Find the roost sites - assign roost as a habitat type.
roost <- gps_dataset %>% 
  st_transform(crs = st_crs(27700)) %>%
  filter(hour(datetime) > 22) %>% # based on time of day in July
  rbind(gps_dataset %>% st_transform(crs = st_crs(27700)) %>% filter(hour(datetime) < 4))

# Plot the roost sites
tm_shape(clachan) + tm_fill() + # base map
  tm_shape(roost) + tm_dots(fill = 'name')

# Identify the four roost sites.
adj <- st_distance(roost)
adj <- matrix(as.numeric(as.numeric(adj)) < 50, nrow = nrow(adj))
g <- igraph::graph_from_adjacency_matrix(adj)
roost <- roost %>% mutate(roost_site = igraph::components(g)$membership)

# Check if roost sites have been assigned correctly.
tm_shape(clachan) + tm_fill() + # base map
  tm_shape(roost) + tm_dots(fill = 'roost_site') # success.

# Draw polygons around the roost sites.
roost_sites <- roost %>% group_by(roost_site) %>% 
  summarise(geometry = st_union(geometry)) %>%
  st_convex_hull() %>%
  st_buffer(10) %>%
  mutate(across(roost_site, as.character)) %>%
  write_sf('bin/roost_sites.gpkg') # add buffer of 10m to the edges.

# Check the identified roost sites
tm_shape(clachan) + tm_fill() + # base map
  tm_shape(roost_sites) + tm_dots(fill = 'roost_site') # success.


############ CREATE AN ANALYSIS DATASET CONTAINING GPS DATA, HABITATS AND ROOST SITES ====

### Because 'arable land' is a big overlapping category in object clachan_hab it has to be modified.
## Modification here seperate all polygons with other habitat type and cuts them out of arable.
## The remaining arable is then stitched back.

clachan_arable <- clachan_hab %>% filter(HABITAT_NA == 'Arable land and market gardens') %>%
  summarise(geometry = st_union(geometry))
clachan_noarable <- clachan_hab %>% filter(HABITAT_NA != 'Arable land and market gardens') %>%
  summarise(geometry = st_union(geometry))
clachan_arable <- st_difference(clachan_arable, clachan_noarable)
clachan_new <- clachan_hab %>% dplyr::select(HABITAT_NA) %>%
  rename(habitat = HABITAT_NA) %>% filter(habitat != 'Arable land and market gardens') %>%
  group_by(habitat) %>% summarise(geometry = st_union(geometry)) %>%
  rbind(clachan_arable %>% mutate(habitat = 'Arable land'))


## Compare the two shapefiles.
t1 <- clachan_new %>% summarise(geometry = st_union(geometry))
t2 <- clachan_hab %>% summarise(geometry = st_union(geometry))

ggplot() + geom_sf(data = t1, fill = 'red') + geom_sf(data = t2, fill = 'green', alpha = 0.5) # all brown; complete overlap.

ggplot() + geom_sf(data = clachan_new, aes(fill = habitat)) # look at the final map to confirm.

clachan_new %>% st_write('bin/clachan_new_habitat.gpkg') # overwrite previous for maping.

## Create the dataset -- notably one point can be classified as more than one habitat
gps_hab <- gps_dataset %>%
  st_transform(crs = st_crs(27700)) %>% # transform to British grid to work with habitat data.
  st_intersection(clachan_new) %>% # obtain habitat information
  st_join(roost_sites %>% st_transform(crs = st_crs(27700))) %>% # add roost site data
  mutate(inroost = ifelse(is.na(roost_site), 'N', 'Y')) %>%
  arrange(datetime) %>% # arrange according to time
  distinct(datetime, geometry, .keep_all = T) # remove duplicate observations; this keeps dunes due to alphabet.
  
## gps_hab has lost some observations due to lack of matched habitat, we will want to keep them.
missing <- anti_join(gps_dataset %>% st_drop_geometry(), gps_hab %>% st_drop_geometry())
missing <- left_join(missing, gps_dataset) %>% st_as_sf() %>%
  st_transform(crs = st_crs(27700)) %>%
  st_join(roost_sites %>% st_transform(crs = st_crs(27700))) %>% # add roost site data
  mutate(inroost = ifelse(is.na(roost_site), 'N', 'Y')) %>%
  mutate(habitat = 'Unknown') %>%
  relocate(habitat, .after = datetime) %>%
  rename("alt.ellipsoid." = "alt(ellipsoid)")

gps_hab <- rbind(gps_hab, missing)
write_sf(gps_hab, 'bin/filtered_gps_habitat_dataset.gpkg')

## Extract lon/lat coordinates for plotting purposes.
gps_dataset_coordinates <- gps_hab %>%
  st_coordinates() %>%
  as.data.frame() %>%
  mutate(name = gps_hab$name, datetime = gps_hab$datetime) %>% 
  write_csv('bin/no_geometry_gps_dataset.csv')

