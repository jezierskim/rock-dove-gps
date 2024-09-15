##### DATA ASSEMBLY AND CLEANING CODE FOR THE GPS PAPER 
### Versioned by MTJ on 9th Jan 2023.

### Necessary packages and working directory.
setwd('~/Desktop/gps revision/code/bin/')

library(sf)
library(tidyverse)
library(ggmap)
library(lubridate)
library(raster)

############ DATA FILTERING AND ORGANISING THE DATASET ====

## Main filtering stage to remove dodgy observations; keeps 4,043 out of 5,937 observations.
gps_dataset <- read_csv('gps_dataset.csv') %>% # load in the file created above.
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
  filter(datetime >= as.Date('2023-07-10')) %>% # remove capture evening.
  distinct(datetime, name, .keep_all = T) # remove 20 duplicated timestamps.

gps_dataset$datetime <- force_tz(gps_dataset$datetime, tzone = 'Europe/Warsaw')
gps_dataset$datetime <- with_tz(gps_dataset$datetime, tzone = 'Europe/London')
write_sf(gps_dataset, 'filtered_gps_dataset.gpkg')


############ GENERATING MAP DATA ====

### Create the clachan base and habitat shape files to be used in the analysis.
clachan <- st_read('~/Desktop/biogeography data/OLD UK DIVA/GBR_adm0.shp') %>% # from ESRI free database.
  st_transform(crs = 27700) %>% # transform to British National Grid
  st_crop(xmin = 85206.7663, ymin = 875224.6196, xmax = 89558.8570, ymax = 878304.2766) %>% # limit to Clachan.
  write_sf('clachan_base.gpkg') 
clachan_hab <- st_read('Habitat_Map_of_Scotland_UIST.gpkg') %>% # from Nature Scot
  st_crop(xmin = 85206.7663, ymin = 875224.6196, xmax = 89558.8570, ymax = 878304.2766) %>% # limit to Clachan
  group_by(HABITAT_NAME) %>% summarise(SHAPE = st_union(SHAPE)) %>% # unionised different habitats.
  write_sf('clachan_habitat.gpkg')

### Obtain satellite base map (NB: this uses horrible ggmap with no proper crsing; difficult to work with).
library(ggmap)
register_google() # needs a paid-for API key
lon <- c(-7.280598, -7.207275) # define longitude of the region
lat <- c(57.65647, 57.67809) # define latitude of the region
points <- st_coordinates(gps) # make them coordinates
map <- get_map(location = c(lon = mean(lon), lat = mean(lat)), zoom = 13, # obtain satellite image.
               maptype = "satellite", source = "google")
save(map, file = 'sat_map_clachan.RData') # save for re-use; NB: crs is 4326 or some other; not British Grid.

############ IDENTIFYING ROOST SITES ====

## Find the roost sites - assign roost as a habitat type.
roost <- gps_dataset %>% 
  filter(hour(datetime) > 22) %>% # based on time of day in July
  rbind(gps_dataset) %>% filter(hour(datetime) < 4)

# Plot the roost sites
ggplot() +
  geom_sf(data = clachan) +
  geom_sf(data = roost)

# Identify the four roost sites.
adj <- st_distance(roost)
adj <- matrix(as.numeric(as.numeric(adj)) < 50, nrow = nrow(adj))
g <- igraph::graph_from_adjacency_matrix(adj)
roost <- roost %>% mutate(roost_site = igraph::components(g)$membership)

# Check if roost sites have been assigned correctly.
ggplot() +
  geom_sf(data = clachan) +
  geom_sf(data = roost, aes(col = as.character(roost_site)))

# Draw polygons around the roost sites.
roost_sites <- roost %>% group_by(roost_site) %>% 
  st_transform(27700) %>%
  summarise(geometry = st_union(geometry)) %>%
  st_convex_hull() %>%
  st_buffer(10) %>%
  mutate(across(roost_site, as.character)) %>%
  st_transform(4326) %>%
  write_sf('roost_sites.gpkg') # add buffer of 10m to the edges.

# Check the identified roost sites
ggplot() +
  geom_sf(data = clachan) +
  geom_sf(data = roost_sites, aes(col = as.character(roost_site)))


############ CREATE AN ANALYSIS DATASET CONTAINING GPS DATA, HABITATS AND ROOST SITES ====

### Because both 'arable land' and 'machair' are big overlapping categories in object clachan_hab it has to be modified.
## Modification here seperate all polygons with other habitat type and cuts them out of arable and machair.
## The remaining arable is then stitched back.

clachan_arable <- clachan_hab %>% filter(HABITAT_NAME == 'Arable land and market gardens') %>%
  summarise(SHAPE = st_union(SHAPE)) %>% st_transform(4326)

ggmap(clachan_sat) +
  geom_sf(data = clachan_arable, inherit.aes = F) +
  coord_sf()

clachan_machair <- clachan_hab %>% filter(HABITAT_NAME == 'H21A0 - Machairs') %>%
  summarise(SHAPE = st_union(SHAPE)) %>% st_transform(4326)

ggmap(clachan_sat) +
  geom_sf(data = clachan_machair, inherit.aes = F) +
  coord_sf()

clachan_noarable <- clachan_hab %>% filter(HABITAT_NAME != 'Arable land and market gardens') %>% 
  filter(HABITAT_NAME != 'H21A0 - Machairs') %>% summarise(SHAPE = st_union(SHAPE)) %>% st_transform(4326)

ggmap(clachan_sat) +
  geom_sf(data = clachan_noarable, inherit.aes = F) +
  coord_sf()

clachan_arable <- st_difference(clachan_arable, clachan_noarable)
clachan_machair <- st_difference(clachan_machair, clachan_noarable)

clachan_arable <- st_difference(clachan_arable, clachan_machair)
clachan_machair <- st_difference(clachan_machair, clachan_arable)

clachan_new <- clachan_hab %>% dplyr::select(HABITAT_NAME) %>% 
  st_transform(4326) %>%
  filter(HABITAT_NAME != 'Arable land and market gardens') %>%
  filter(HABITAT_NAME != 'H21A0 - Machairs') %>%
  group_by(HABITAT_NAME) %>% summarise(SHAPE = st_union(SHAPE)) %>%
  rbind(clachan_arable %>% mutate(HABITAT_NAME = 'Arable land')) %>%
  rbind(clachan_machair %>% mutate(HABITAT_NAME = 'Machair')) %>%
  rename(geometry = SHAPE, habitat = HABITAT_NAME)


## Compare the two shapefiles.
t1 <- clachan_new %>% summarise(geometry = st_union(geometry))
t2 <- clachan_hab %>% summarise(SHAPE = st_union(SHAPE))

ggplot() + geom_sf(data = t1, fill = 'red') + geom_sf(data = t2, fill = 'green', alpha = 0.5) # all brown; complete overlap.

ggplot() + geom_sf(data = clachan_new, aes(fill = habitat)) # look at the final map to confirm.

clachan_new %>% st_write('clachan_new_habitat.gpkg') # overwrite previous for maping.

## Create the dataset -- notably one point can be classified as more than one habitat
gps_hab <- gps_dataset %>%
  st_intersection(clachan_new) %>% # obtain habitat information
  st_join(roost_sites) %>% # add roost site data
  mutate(inroost = ifelse(is.na(roost_site), 'N', 'Y')) %>%
  arrange(datetime) %>% # arrange according to time
  distinct(datetime, geometry, .keep_all = T) # remove duplicate observations; this keeps dunes due to alphabet.
  
## gps_hab has lost some observations due to lack of matched habitat, we will want to keep them.
missing <- anti_join(gps_dataset %>% st_drop_geometry(), gps_hab %>% st_drop_geometry())
missing <- left_join(missing, gps_dataset) %>% st_as_sf() %>%
  st_join(roost_sites) %>% # add roost site data
  mutate(inroost = ifelse(is.na(roost_site), 'N', 'Y')) %>%
  mutate(habitat = 'Unknown') %>%
  relocate(habitat, .after = datetime) %>%
  rename("alt.ellipsoid." = "alt(ellipsoid)")

ggmap(clachan_sat) +
  geom_sf(data = missing, inherit.aes = F)

gps_hab <- rbind(gps_hab, missing)
write_sf(gps_hab, 'filtered_gps_habitat_dataset.gpkg')

## Extract lon/lat coordinates for plotting purposes.
gps_dataset_coordinates <- gps_hab %>%
  st_coordinates() %>%
  as.data.frame() %>%
  mutate(name = gps_hab$name, datetime = gps_hab$datetime) %>% 
  write_csv('no_geometry_gps_dataset.csv')

