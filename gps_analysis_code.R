##### GPS ANALYSIS CODE #####
### To accompany the manuscript: Use of anthropogenic landscapes in a wild Rock Dove Columba livia population
### William J. Smith, Stephen Portugal, Michał T. Jezierski.
## Versioned by MTJ 4th April 2024.

### 1 == SET UP THE DIRECTORY, LOAD PACKAGES AND DATASETS REQUIRED. ####################################
#### SET WORKING DIRECTORY. ##
setwd('~/Desktop/rock doves/gps_folder/')

#### PACKAGES NEEDED. ##
library(tidyverse)
library(sf)
library(lubridate)
library(ggmap)
library(cowplot)
library(adehabitatHR)
library(cowplot)
library(units)

sf_use_s2(F) # also turn of spherical geometry for sf.

#### LOAD IN THE DATASET. ##
## Data as prepared in the file 'gps_data_code.R'.
gps <- st_read('bin/filtered_gps_habitat_dataset.gpkg') %>% # the filtered dataset; see file above and METHODS.
  mutate(date = as.Date((format(as.POSIXct(datetime), "%Y-%m-%d")))) # reformat date so it can be used more easily.
clachan_base <- st_read('bin/clachan_base.gpkg') # base map of the study area.
clachan_habitat <- st_read('bin/clachan_new_habitat.gpkg') # habitat map of the study area.
load('sat_map_clachan.RData') # loading satellite map of the study area; execute all three lines.
clachan_sat <- map
rm(map) # map loaded correctly by this point.




### 2 == GENERAL STATISTICS; DISTANCE TRAVELLED ####################################

#### HOW MANY TRACKING DAYS ACROSS THE DATASET? ##
gps %>% st_drop_geometry() %>% # take the dataset; remove geometry (not needed for calculations).
  group_by(name) %>% # group across all individual rock doves.
  summarise(start_date = min(date), end_date = max(date)) %>% # summarise the minimum and maximum date.
  mutate(date_range = (end_date - start_date)+1) %>% # calculate the date range.
  print() %>% # print it for each pigeon name.
  summarise(mean = mean(date_range), sd = sd(date_range)) # and across all rock dove name, give mean and SD.

### HOW FAR DID EACH ROCK DOVE TRAVEL PER DAY? ##
### FIRST: calculate the total daily movement of each Rock Dove.
gps_dist <- gps %>% 
  group_by(name, date) %>% # per pigeon and date
  mutate(
    lead = geom[row_number() + 1], # create a lead variable for calculation
    dist = st_distance(geom, lead, by_element = T)) %>% # calculate the distance between the current point and the next one. 
  drop_units() %>% # remove units (obstruct calculations)
  st_drop_geometry() %>% # drop geometry (not needed)
  drop_na(dist) %>% # drop empty rows which result from adding lead to the last observation per grouping.
  summarise(total = sum(dist)) # summarise distance 

### SECOND: calculate the mean and standard deviation for each Rock Dove
gps_dist %>% group_by(name) %>% summarise(mean_dist = mean(total), sd_dist = sd(total)) %>% print()

figS1 <- ggplot(data = gps_dist, aes(x = name, y = total)) + ggdist::stat_halfeye(justification = -.2, point_colour = NA, .width = 0) +
  geom_point(size = 1, position = position_nudge(x = -.2)) + 
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar", color="red", width = 0.1) +
  stat_summary(fun.y=mean, geom="point", color="red", size = 2) +
  theme_bw() +
  scale_y_continuous(name = 'Total distance travelled per day (m)') +
  scale_x_discrete(name = 'Rock Dove ID') +
  theme(panel.grid = element_blank())

ggsave(figS1, path = 'figs/', filename = 'FIGS1.png', dpi = 300, width = 7, height = 4.5, units = 'in')


### 3 == HOME RANGES ####################################
### FIRST: Estimate Rock Dove home ranges. 
## Using adehabitatHR core functions; there are two calculations: one based on fixed-kernel utilisation distribution
## the other is based on minimum convex polygon
spatial_gps <- gps %>% dplyr::select(name) %>% as_Spatial() # turn the observations into spatial (other data format)
spatial_gps$name <- as.factor(spatial_gps$name) # make names as factors as required for the analysis.

## Calculate fixed kernel utilisation distribution.
kud <- kernelUD(spatial_gps, h="href") # estimates of Kernel Home-Range
(ver_95 <- getverticeshr(kud, 95)) # extract the 95% contour of home range
(ver_50 <- getverticeshr(kud, 50)) # extract the 50% contour of home range

## Calculate minimum convex polygon.
(mcp_100 <- mcp(spatial_gps, 100))

## Summary of pigeon tracking data, including the different home ranges.
options(pillar.sigfig = 2) # forces significant figures.
gps %>% st_drop_geometry() %>% # drop geometry
  group_by(name) %>% # across all pigeon names
  summarise(start_date = min(date), end_date = max(date), fix = n()) %>% # recreate the dates above + no. of fixes.
  mutate(date_range = end_date - start_date) %>% 
  mutate(sex = c('Male', 'Male', 'Male', 'Female (?)', 'Female'), # add bonus information
         age = c('Adult', 'Adult', 'Adult', 'Adult', 'Adult'),
         home_ud_50_ha = num(c(5.485856, 34.322273, 14.104215, 78.748625, 6.729580), digits = 3), # add data from HR estimation.
         home_ud_95_ha = num(c(68.37009, 134.09548, 62.03079, 308.29462, 38.66206), digits = 3),
         home_mcp_100_ha = num(c(218.68409, 129.18426, 34.62694, 206.07421, 56.46548), digits = 3)) %>%
  print() %>%
  summarise(mean_ud_50 = mean(home_ud_50_ha), sd_ud_50 = sd(home_ud_50_ha), mean_ud_95 = mean(home_ud_95_ha), sd_ud_95 = sd(home_ud_95_ha),
            mean_mcp_100 = mean(home_mcp_100_ha), sd_mcp_100 = sd(home_mcp_100_ha))


## SECOND: compare Rock Dove ranges to data from Rose et al. 2006 (MCP) and Carlson et al. 2011 (UD)
### DATA FROM ROSE ET AL. 2006:
rose_data <- c(2.9, 2.9, 3.2, 5.1, 5.2,12.5, 14.8, 16.3, 24.2, 46.7, 144.1, 150.6)
mean(rose_data)
sd(rose_data)
wilcox.test(rose_data, mcp_100$area)

### DATA FROM CARLSON ET AL. 2011
carlson_data <- c(9.263, 0.015, 0.194, 0.148, 0.619, 1.385, 0.055, 0.037, 1.233, 0.153, 14.926, 16.830, 8.748, 3.323)
carlson_data <- carlson_data*100 # convert to hectares from km2.
mean(carlson_data)
sd(carlson_data)
wilcox.test(carlson_data, mcp_100$area)

### 4 == FIGURE 2 - MOVEMENT OF ROCK DOVES THROUGH LANDSCAPE ####################################
### GPS dataset needs to be fully compatible with WGS coordinate system.
gps_wgs_coords <- gps %>% # create a dataset with X/Y rather than geometry coordinates.
  st_transform(crs = st_crs(4326)) %>% # must recreate coordinates for lon/lat projections
  st_coordinates() %>% # # extract coordinates
  as.data.frame() %>% # as data frame
  mutate(name = gps$name, datetime = gps$datetime) # add name and datetime columns

### Load inset North Uist map.
inset_nu <- st_read('bin/GBR_adm0.shp') %>% # shapefile from https://www.diva-gis.org/gdata
  st_crop(xmin = -7.577820, ymin = 57.464158, xmax = -7.009277, ymax = 57.747412) # crop to NUist box

### Create the map of North Uist.
(inset_nu_map <- ggplot() + 
  geom_sf(data = inset_nu) + # add sf object called inse nu
  geom_rect(aes(xmin = -7.280598, 
                ymin = 57.65647,
                xmax = -7.207275,
                ymax = 57.68704), col = 'black', fill = NA) + # rectangle around Clachan.
  theme_void() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()))

### Create Figure 1A to which the inset will be added.
(fig2A <- ggmap(clachan_sat) + 
  geom_point(data = gps_wgs_coords, aes(x = X, y = Y, col = name), size = 0.3) + # add points
  geom_path(data = gps_wgs_coords, aes(x = X, y = Y, col = name), size = 0.3) + # add paths
  xlab('Longitude (°)') + ylab('Latitude (°)') + # plot labels
  scale_y_continuous(limits = c(57.65,57.68408)) + # cropping scale; satellite is larger than clachan
  scale_color_discrete(name = 'Dove ID') + # here you can add colours you like using values = c('colou1', 'colour2')
  ggsn::scalebar(y.min = 57.652, x.min = -7.298783, y.max = 57.68, x.max = -7.198, # adds the scale bar
          transform = T, dist = 500, dist_unit = 'm', location = 'topright',
          box.fill = c("yellow", "white"), height = 0.02, st.dist = 0.05, st.color = 'white', st.size = 10/.pt) +
  geom_segment(aes(x = -7.29, y = 57.68, xend = -7.29, yend = 57.683), linewidth = 1, col = 'white',
                arrow = arrow(length = unit(0.25, "cm"))) + # adds the north arrow
  geom_text(aes(x = -7.29, y = 57.684), label = 'N', col = 'white', size = 14/.pt) + # add the north label
  coord_cartesian())

(fig2 <- ggdraw() + # code joining the two figures
  draw_plot(fig1A) +
  draw_plot(inset_nu_map, x = 0.79, y = 0.1, width = 0.25, height = 0.25))

ggsave(fig2, path = 'figs/', filename = 'FIG2.png', dpi = 300, width = 7, height = 4.5, units = 'in')

### 5 == FIGURE S2-S6 - MOVEMENT OF EACH PIGEON AND HOME RANGE. #################################### 
## Need to obtain polygons of home ranges first.
kud_names <- names(kud) # define names from the estimation.
d_95 <- lapply(kud, function(x) try(getverticeshr(x, 95))) # extract home ranges at 95%
sapply(1:length(d_95), function(i) {
  row.names(d_95[[i]]) <- kud_names[i] # assign them names
})
sdf_d_95<- Reduce(rbind, d_95)
sdf_d_95 <- spTransform(sdf_d_95, CRS("+init=epsg:4326")) # projection in degrees.
sdf_d_95 <- fortify(sdf_d_95) # turn to data frame so ggplot2 can use.

## Repeat for 50
kud_names <- names(kud) # define names
d_50 <- lapply(kud, function(x) try(getverticeshr(x, 50))) # extract home ranges at 50%
sapply(1:length(d_50), function(i) {
  row.names(d_50[[i]]) <- kud_names[i] # assign them names
})
sdf_d_50 <- Reduce(rbind, d_50)
sdf_d_50 <- spTransform(sdf_d_50, CRS("+init=epsg:4326"))
sdf_d_50 <- fortify(sdf_d_50)

#### FIGURES1 - EA49502.
## Plot as above repeated with each polygon recreated.
(figS2A <- ggmap(clachan_sat) + 
    geom_point(data = gps_wgs_coords %>% filter(name == 'EA49502'), aes(x = X, y = Y, col = name), size = 0.3) + 
    geom_path(data = gps_wgs_coords %>% filter(name == 'EA49502'), aes(x = X, y = Y, col = name), size = 0.3) +
    geom_polygon(data = sdf_d_50 %>% filter(id == 'homerange'), aes(x = long, y = lat, fill = id, group = group), col = 'red', alpha = 0.6) +
    xlab('Longitude (°)') + ylab('Latitude (°)') + 
    scale_y_continuous(limits = c(57.65,57.68408)) +
    scale_color_discrete(name = 'Dove ID') +
    ggsn::scalebar(y.min = 57.652, x.min = -7.298783, y.max = 57.68, x.max = -7.198,
             transform = T, dist = 500, dist_unit = 'm', location = 'topright',
             box.fill = c("yellow", "white"), height = 0.02, st.dist = 0.02, st.color = 'white', st.size = 8/.pt) +
    geom_segment(aes(x = -7.29, y = 57.68, xend = -7.29, yend = 57.683), linewidth = 1, col = 'white',
                 arrow = arrow(length = unit(0.25, "cm"))) + 
    geom_text(aes(x = -7.29, y = 57.684), label = 'N', col = 'white', size = 14/.pt) +
    geom_label(aes(y = 57.65, x = -7.27), label = 'EA49502; 50% HR', size = 14/.pt) +
    coord_cartesian() +
    theme(
          axis.title = element_blank(),
          legend.position = 'none'))
(figS2B <- ggmap(clachan_sat) + 
    geom_point(data = gps_wgs_coords %>% filter(name == 'EA49502'), aes(x = X, y = Y, col = name), size = 0.3) + 
    geom_path(data = gps_wgs_coords %>% filter(name == 'EA49502'), aes(x = X, y = Y, col = name), size = 0.3) +
    geom_polygon(data = sdf_d_95 %>% filter(id == 'homerange'), aes(x = long, y = lat, fill = id, group = group), col = 'red', alpha = 0.6) +
    xlab('Longitude (°)') + ylab('Latitude (°)') + 
    scale_y_continuous(limits = c(57.65,57.68408)) +
    scale_color_discrete(name = 'Dove ID') +
    ggsn::scalebar(y.min = 57.652, x.min = -7.298783, y.max = 57.68, x.max = -7.198,
             transform = T, dist = 500, dist_unit = 'm', location = 'topright',
             box.fill = c("yellow", "white"), height = 0.02, st.dist = 0.02, st.color = 'white', st.size = 8/.pt) +
    geom_segment(aes(x = -7.29, y = 57.68, xend = -7.29, yend = 57.683), linewidth = 1, col = 'white',
                 arrow = arrow(length = unit(0.25, "cm"))) + 
    geom_text(aes(x = -7.29, y = 57.684), label = 'N', col = 'white', size = 14/.pt) +
    geom_label(aes(y = 57.65, x = -7.27), label = 'EA49502; 95% HR', size = 14/.pt) +
    coord_cartesian() +
    theme(
          axis.title = element_blank(),
          legend.position = 'none'))

(figS2 <- plot_grid(figS2A, figS2B, labels = c('A)', 'B)'), label_fontfamily = 'serif', label_size = 10))
ggsave(figS2, path = 'figs/', filename = 'figS2.png', dpi = 300, width = 9, height = 4.5, units = 'in')

#### FIGURE S2 - EA49503.

(figS3A <- ggmap(clachan_sat) + 
    geom_point(data = gps_wgs_coords %>% filter(name == 'EA49503'), aes(x = X, y = Y), col = '#98993C', size = 0.3) + 
    geom_path(data = gps_wgs_coords %>% filter(name == 'EA49503'), aes(x = X, y = Y), col = '#98993C', size = 0.3) +
    geom_polygon(data = sdf_d_50 %>% filter(id == 'homerange1'), aes(x = long, y = lat, group = group), fill = '#98993C', col = '#98993C', alpha = 0.6) +
    xlab('Longitude (°)') + ylab('Latitude (°)') + 
    scale_y_continuous(limits = c(57.65,57.68408)) +
    scale_color_discrete(name = 'Dove ID') +
    ggsn::scalebar(y.min = 57.652, x.min = -7.298783, y.max = 57.68, x.max = -7.198,
             transform = T, dist = 500, dist_unit = 'm', location = 'topright',
             box.fill = c("yellow", "white"), height = 0.02, st.dist = 0.02, st.color = 'white', st.size = 8/.pt) +
    geom_segment(aes(x = -7.29, y = 57.68, xend = -7.29, yend = 57.683), linewidth = 1, col = 'white',
                 arrow = arrow(length = unit(0.25, "cm"))) + 
    geom_text(aes(x = -7.29, y = 57.684), label = 'N', col = 'white', size = 14/.pt) +
    geom_label(aes(y = 57.65, x = -7.27), label = 'EA49503; 50% HR', size = 14/.pt) +
    coord_cartesian() +
    theme(
          axis.title = element_blank(),
          legend.position = 'none'))
(figS3B <- ggmap(clachan_sat) + 
    geom_point(data = gps_wgs_coords %>% filter(name == 'EA49503'), aes(x = X, y = Y), col = '#98993C', size = 0.3) + 
    geom_path(data = gps_wgs_coords %>% filter(name == 'EA49503'), aes(x = X, y = Y), col = '#98993C', size = 0.3) +
    geom_polygon(data = sdf_d_95 %>% filter(id == 'homerange1'), aes(x = long, y = lat, group = group), fill = '#98993C', col = '#98993C', alpha = 0.6) +
    xlab('Longitude (°)') + ylab('Latitude (°)') + 
    scale_y_continuous(limits = c(57.65,57.68408)) +
    scale_color_discrete(name = 'Dove ID') +
    ggsn::scalebar(y.min = 57.652, x.min = -7.298783, y.max = 57.68, x.max = -7.198,
             transform = T, dist = 500, dist_unit = 'm', location = 'topright',
             box.fill = c("yellow", "white"), height = 0.02, st.dist = 0.02, st.color = 'white', st.size = 8/.pt) +
    geom_segment(aes(x = -7.29, y = 57.68, xend = -7.29, yend = 57.683), linewidth = 1, col = 'white',
                 arrow = arrow(length = unit(0.25, "cm"))) + 
    geom_text(aes(x = -7.29, y = 57.684), label = 'N', col = 'white', size = 14/.pt) +
    geom_label(aes(y = 57.65, x = -7.27), label = 'EA49503; 95% HR', size = 14/.pt) +
    coord_cartesian() +
    theme(
          axis.title = element_blank(),
          legend.position = 'none'))

(figS3 <- plot_grid(figS3A, figS3B, labels = c('A)', 'B)'), label_fontfamily = 'serif', label_size = 10))
ggsave(figS3, path = 'figs/', filename = 'figS3.png', dpi = 300, width = 9, height = 4.5, units = 'in')

#### FIGURE S3 - EA49568.

(figS4A <- ggmap(clachan_sat) + 
    geom_point(data = gps_wgs_coords %>% filter(name == 'EA49568'), aes(x = X, y = Y), col = '#39A783', size = 0.3) + 
    geom_path(data = gps_wgs_coords %>% filter(name == 'EA49568'), aes(x = X, y = Y), col = '#39A783', size = 0.3) +
    geom_polygon(data = sdf_d_50 %>% filter(id == 'homerange2'), aes(x = long, y = lat, group = group), fill = '#39A783', col = '#39A783', alpha = 0.6) +
    xlab('Longitude (°)') + ylab('Latitude (°)') + 
    scale_y_continuous(limits = c(57.65,57.68408)) +
    scale_color_discrete(name = 'Dove ID') +
    ggsn::scalebar(y.min = 57.652, x.min = -7.298783, y.max = 57.68, x.max = -7.198,
             transform = T, dist = 500, dist_unit = 'm', location = 'topright',
             box.fill = c("yellow", "white"), height = 0.02, st.dist = 0.02, st.color = 'white', st.size = 8/.pt) +
    geom_segment(aes(x = -7.29, y = 57.68, xend = -7.29, yend = 57.683), linewidth = 1, col = 'white',
                 arrow = arrow(length = unit(0.25, "cm"))) + 
    geom_text(aes(x = -7.29, y = 57.684), label = 'N', col = 'white', size = 14/.pt) +
    geom_label(aes(y = 57.65, x = -7.27), label = 'EA49568; 50% HR', size = 14/.pt) +
    coord_cartesian() +
    theme(
          axis.title = element_blank(),
          legend.position = 'none'))
(figS4B <- ggmap(clachan_sat) + 
    geom_point(data = gps_wgs_coords %>% filter(name == 'EA49568'), aes(x = X, y = Y), col = '#39A783', size = 0.3) + 
    geom_path(data = gps_wgs_coords %>% filter(name == 'EA49568'), aes(x = X, y = Y), col = '#39A783', size = 0.3) +
    geom_polygon(data = sdf_d_95 %>% filter(id == 'homerange2'), aes(x = long, y = lat, group = group), fill = '#39A783', col = '#39A783', alpha = 0.6) +
    xlab('Longitude (°)') + ylab('Latitude (°)') + 
    scale_y_continuous(limits = c(57.65,57.68408)) +
    scale_color_discrete(name = 'Dove ID') +
    ggsn::scalebar(y.min = 57.652, x.min = -7.298783, y.max = 57.68, x.max = -7.198,
             transform = T, dist = 500, dist_unit = 'm', location = 'topright',
             box.fill = c("yellow", "white"), height = 0.02, st.dist = 0.02, st.color = 'white', st.size = 8/.pt) +
    geom_segment(aes(x = -7.29, y = 57.68, xend = -7.29, yend = 57.683), linewidth = 1, col = 'white',
                 arrow = arrow(length = unit(0.25, "cm"))) + 
    geom_text(aes(x = -7.29, y = 57.684), label = 'N', col = 'white', size = 14/.pt) +
    geom_label(aes(y = 57.65, x = -7.27), label = 'EA49568; 95% HR', size = 14/.pt) +
    coord_cartesian() +
    theme(
          axis.title = element_blank(),
          legend.position = 'none'))

(figS4 <- plot_grid(figS4A, figS4B, labels = c('A)', 'B)'), label_fontfamily = 'serif', label_size = 10))
ggsave(figS4, path = 'figs/', filename = 'figS4.png', dpi = 300, width = 9, height = 4.5, units = 'in')

#### FIGURE S4 - EA49570.

(figS5A <- ggmap(clachan_sat) + 
    geom_point(data = gps_wgs_coords %>% filter(name == 'EA49570'), aes(x = X, y = Y), col = '#3CABD3', size = 0.3) + 
    geom_path(data = gps_wgs_coords %>% filter(name == 'EA49570'), aes(x = X, y = Y), col = '#3CABD3', size = 0.3) +
    geom_polygon(data = sdf_d_50 %>% filter(id == 'homerange3'), aes(x = long, y = lat, group = group), fill = '#3CABD3', col = '#3CABD3', alpha = 0.6) +
    xlab('Longitude (°)') + ylab('Latitude (°)') + 
    scale_y_continuous(limits = c(57.65,57.68408)) +
    scale_color_discrete(name = 'Dove ID') +
    ggsn::scalebar(y.min = 57.652, x.min = -7.298783, y.max = 57.68, x.max = -7.198,
             transform = T, dist = 500, dist_unit = 'm', location = 'topright',
             box.fill = c("yellow", "white"), height = 0.02, st.dist = 0.02, st.color = 'white', st.size = 8/.pt) +
    geom_segment(aes(x = -7.29, y = 57.68, xend = -7.29, yend = 57.683), linewidth = 1, col = 'white',
                 arrow = arrow(length = unit(0.25, "cm"))) + 
    geom_text(aes(x = -7.29, y = 57.684), label = 'N', col = 'white', size = 14/.pt) +
    geom_label(aes(y = 57.65, x = -7.27), label = 'EA49570; 50% HR', size = 14/.pt) +
    coord_cartesian() +
    theme(
          axis.title = element_blank(),
          legend.position = 'none'))
(figS5B <- ggmap(clachan_sat) + 
    geom_point(data = gps_wgs_coords %>% filter(name == 'EA49570'), aes(x = X, y = Y), col = '#3CABD3', size = 0.3) + 
    geom_path(data = gps_wgs_coords %>% filter(name == 'EA49570'), aes(x = X, y = Y), col = '#3CABD3', size = 0.3) +
    geom_polygon(data = sdf_d_95 %>% filter(id == 'homerange3'), aes(x = long, y = lat, group = group), fill = '#3CABD3', col = '#3CABD3', alpha = 0.6) +
    xlab('Longitude (°)') + ylab('Latitude (°)') + 
    scale_y_continuous(limits = c(57.65,57.68408)) +
    scale_color_discrete(name = 'Dove ID') +
    ggsn::scalebar(y.min = 57.652, x.min = -7.298783, y.max = 57.68, x.max = -7.198,
             transform = T, dist = 500, dist_unit = 'm', location = 'topright',
             box.fill = c("yellow", "white"), height = 0.02, st.dist = 0.02, st.color = 'white', st.size = 8/.pt) +
    geom_segment(aes(x = -7.29, y = 57.68, xend = -7.29, yend = 57.683), linewidth = 1, col = 'white',
                 arrow = arrow(length = unit(0.25, "cm"))) + 
    geom_text(aes(x = -7.29, y = 57.684), label = 'N', col = 'white', size = 14/.pt) +
    geom_label(aes(y = 57.65, x = -7.27), label = 'EA49570; 95% HR', size = 14/.pt) +
    coord_cartesian() +
    theme(
          axis.title = element_blank(),
          legend.position = 'none'))

(figS5 <- plot_grid(figS5A, figS5B, labels = c('A)', 'B)'), label_fontfamily = 'serif', label_size = 10))
ggsave(figS5, path = 'figs/', filename = 'figS5.png', dpi = 300, width = 9, height = 4.5, units = 'in')

#### FIGURE S5 - EA49571.

(figS6A <- ggmap(clachan_sat) + 
    geom_point(data = gps_wgs_coords %>% filter(name == 'EA49571'), aes(x = X, y = Y), col = '#E769F1', size = 0.3) + 
    geom_path(data = gps_wgs_coords %>% filter(name == 'EA49571'), aes(x = X, y = Y), col = '#E769F1', size = 0.3) +
    geom_polygon(data = sdf_d_50 %>% filter(id == 'homerange4'), aes(x = long, y = lat, group = group), fill = '#E769F1', col = '#E769F1', alpha = 0.6) +
    xlab('Longitude (°)') + ylab('Latitude (°)') + 
    scale_y_continuous(limits = c(57.65,57.68408)) +
    scale_color_discrete(name = 'Dove ID') +
    ggsn::scalebar(y.min = 57.652, x.min = -7.298783, y.max = 57.68, x.max = -7.198,
             transform = T, dist = 500, dist_unit = 'm', location = 'topright',
             box.fill = c("yellow", "white"), height = 0.02, st.dist = 0.02, st.color = 'white', st.size = 8/.pt) +
    geom_segment(aes(x = -7.29, y = 57.68, xend = -7.29, yend = 57.683), linewidth = 1, col = 'white',
                 arrow = arrow(length = unit(0.25, "cm"))) + 
    geom_text(aes(x = -7.29, y = 57.684), label = 'N', col = 'white', size = 14/.pt) +
    geom_label(aes(y = 57.65, x = -7.27), label = 'EA49571; 50% HR', size = 14/.pt) +
    coord_cartesian() +
    theme(
          axis.title = element_blank(),
          legend.position = 'none'))
(figS6B <- ggmap(clachan_sat) + 
    geom_point(data = gps_wgs_coords %>% filter(name == 'EA49571'), aes(x = X, y = Y), col = '#E769F1', size = 0.3) + 
    geom_path(data = gps_wgs_coords %>% filter(name == 'EA49571'), aes(x = X, y = Y), col = '#E769F1', size = 0.3) +
    geom_polygon(data = sdf_d_95 %>% filter(id == 'homerange4'), aes(x = long, y = lat, group = group), fill = '#E769F1', col = '#E769F1', alpha = 0.6) +
    xlab('Longitude (°)') + ylab('Latitude (°)') + 
    scale_y_continuous(limits = c(57.65,57.68408)) +
    scale_color_discrete(name = 'Dove ID') +
    ggsn::scalebar(y.min = 57.652, x.min = -7.298783, y.max = 57.68, x.max = -7.198,
             transform = T, dist = 500, dist_unit = 'm', location = 'topright',
             box.fill = c("yellow", "white"), height = 0.02, st.dist = 0.02, st.color = 'white', st.size = 8/.pt) +
    geom_segment(aes(x = -7.29, y = 57.68, xend = -7.29, yend = 57.683), linewidth = 1, col = 'white',
                 arrow = arrow(length = unit(0.25, "cm"))) + 
    geom_text(aes(x = -7.29, y = 57.684), label = 'N', col = 'white', size = 14/.pt) +
    geom_label(aes(y = 57.65, x = -7.27), label = 'EA49571; 95% HR', size = 14/.pt) +
    coord_cartesian() +
    theme(
          axis.title = element_blank(),
          legend.position = 'none'))

(figS6 <- plot_grid(figS6A, figS6B, labels = c('A)', 'B)'), label_fontfamily = 'serif', label_size = 10))
ggsave(figS6, path = 'figs/', filename = 'figS6.png', dpi = 300, width = 9, height = 4.5, units = 'in')

### 6 == HABITAT USE BY ROCK DOVES ####################################
## Approach taken here: obtain a dataset of 0.5h periods across all monitoring dates; then check which
## GPS points happened during which period; for each dove, keep the one it spent the most time in.
# Create a sequence of half an hour slots
full_h <- as.data.frame(with(expand.grid(date = as.Date(unique(gps$date)), hour = paste0(0:23, ":00:00")), 
     sort(as.POSIXct(paste(date, hour)))))  # create a dataset of full hours
colnames(full_h) <- 'start' # call it start as it will be the first value
half_h <- as.data.frame(with(expand.grid(date = as.Date(unique(gps$date)), hour = paste0(0:23, ":29:59")), 
                             sort(as.POSIXct(paste(date, hour))))) # create a dataset of half-hours meant for end.
colnames(half_h) <- 'end' # call it end.
date_df1 <- cbind(full_h, half_h) # bind them together; gives you every XX:00:00 - XX:29:59 period.
full_h1 <- as.data.frame(with(expand.grid(date = as.Date(unique(gps$date)), hour = paste0(0:23, ":59:59")), 
                              sort(as.POSIXct(paste(date, hour))))) # repeat for ends of fulls hours.
colnames(full_h1) <- 'end'
half_h1 <- as.data.frame(with(expand.grid(date = as.Date(unique(gps$date)), hour = paste0(0:23, ":30:00")), 
                              sort(as.POSIXct(paste(date, hour))))) # repeat for starts of half hours.
colnames(half_h1) <- 'start'
date_df2 <- cbind(half_h1, full_h1) # bind them together; gives you every XX:30:00 - XX:59:59 period.

date_df <- rbind(date_df1, date_df2) %>% arrange(start) %>% # merge all together and arrange
  mutate(int = interval(start, end)) %>% # create a time interval column.
  mutate(id = rep(seq(1, 48, 1), 48, 1392)) # give each 0.5h period ID that is the same across days.

# Now need to check if a given pigeon observation falls within that time interval.
gps_int <- date_df %>%
  left_join(gps %>% filter(inroost == 'N'), by = character()) %>% # cross-join while removing in-roost observations.
  filter(datetime %within% int) # verify if a timepoint is within an interval.

# Simplify habitats:
gps_int <- gps_int %>% mutate(habitat = recode(habitat, "H21A0 - Machairs" = 'Machair', "H2130 - Fixed dunes" = "Dunes",
                                               "Atlantic Cynosurus-Centaurea pastures" = "Atlantic pasture",
                                               "H7140 - Transition mires" = "Unknown or other",
                                               "Unknown" = "Unknown or other",
                                               "H2120 - Shifting dunes" = "Dunes",
                                               "Agriculturally-improved, re-seeded and heavily fertilised grassland, including sports fields and grass lawns" = "Arable land",
                                               "Constructed, industrial and other artificial habitats" = "Unknown or other",
                                               "Unvegetated sand beaches above the driftline" = "Unknown or other",
                                               "H2190 - Humid dune slacks" = 'Dunes',
                                               "Flood swards and related communities" = "Unknown or other"))

# Time spent in habitat.
time_in_habitat <- gps_int %>%
  group_by(name, date, id, habitat) %>% # for each pigeon, on each day, across each habitat and interval.
  count() %>% # count how many times each pigeon on a given day used which habitats in which interval.
  ungroup() %>% # remove grouping
  group_by(name, date, id) %>% # now for each pigeon on a given day and time of day.
  slice_max(n) %>% # keep only the habitat with the highest number of points in it; keeps one obs per habitat at that date and time.
  ungroup %>% # remove grouping
  dplyr::select(-n) %>% # remove counter
  group_by(id, habitat) %>% # now by each mix of time and habitat
  count() # count how many there are!

time_in_habitat$id <- as.numeric(time_in_habitat$id) # modify id of time interval to numeric
hours <- date_df %>% dplyr::select(start, id) %>% # create actual hours object to add it as a column
  mutate(hour = paste(format(as.POSIXct(start), format = "%H:%M:%S"))) %>%
  dplyr::select(hour, id) %>%
  unique()
time_in_habitat <- left_join(time_in_habitat, hours, by = 'id') # merge the two together.

### Proportions of time spent in habitat per hour

time_in_habitat <- time_in_habitat %>%
  group_by(hour) %>%
  mutate(prop = n/sum(n))


### 7 == FIGURE 3  ###################################
gps_int_roost <- date_df %>%
  left_join(gps %>% 
              mutate(habitat = ifelse(inroost == 'Y', 'Roost site', gps$habitat)), by = character()) %>% # cross-join add creating roost as a new habitat.
  filter(datetime %within% int) # verify if a timepoint is within an interval.

gps_int_roost <- gps_int_roost %>% mutate(habitat = recode(habitat, "H21A0 - Machairs" = 'Machair', "H2130 - Fixed dunes" = "Dunes",
                                                           "Atlantic Cynosurus-Centaurea pastures" = "Atlantic pasture",
                                                           "H7140 - Transition mires" = "Unknown or other",
                                                           "Unknown" = "Unknown or other",
                                                           "H2120 - Shifting dunes" = "Dunes",
                                                           "Agriculturally-improved, re-seeded and heavily fertilised grassland, including sports fields and grass lawns" = "Arable land",
                                                           "Constructed, industrial and other artificial habitats" = "Unknown or other",
                                                           "Unvegetated sand beaches above the driftline" = "Unknown or other",
                                                           "H2190 - Humid dune slacks" = 'Dunes',
                                                           "Flood swards and related communities" = "Unknown or other"))

unique(gps_int_roost$habitat)

# Time spent in habitat; NB: this is a clever pipe; NB: slice_max dunno what it does if it finds equal number.
time_in_habitat_roost <- gps_int_roost %>%
  group_by(name, date, id, habitat) %>% # for each pigeon, on each day, across each habitat and interval.
  count() %>% # count how many times each pigeon on a given day used which habitats in which interval.
  ungroup() %>% # remove grouping
  group_by(name, date, id) %>% # now for each pigeon on a given day and time of day.
  slice_max(n) %>% # keep only the habitat with the highest number of points in it; keeps one obs per habitat at that date and time.
  ungroup %>% # remove grouping
  dplyr::select(-n) %>% # remove counter
  group_by(id, habitat) %>% # now by each mix of time and habitat
  count() # count how many there are!

time_in_habitat_roost$id <- as.numeric(time_in_habitat_roost$id)
time_in_habitat_roost <- left_join(time_in_habitat_roost, hours, by = 'id') %>%
  filter(habitat != 'unknown')

time_in_habitat_roost <- time_in_habitat_roost %>%
  group_by(hour) %>%
  mutate(prop = n/sum(n))


### Actual figure code.
(fig3 <- ggplot(data = time_in_habitat_roost, aes(x = as.POSIXct(hour, format = "%H:%M:%S"), y = n, fill = habitat)) + 
    geom_bar(position = 'fill', stat = 'identity') + theme_bw() +
    xlab('Time of day (half-hour intervals)') + ylab('Proportion of time in habitat') +
    scale_x_datetime(labels = scales::date_format("%H:%M:%S")) +
    scale_fill_manual(name = '', values = c("Arable land" = "#DDCC77", "Atlantic pasture" = '#CC6677', "Roost site" = 'grey',
                                            "Machair" = '#117733', "Dunes" = '#6699CC', 'Unknown or other' = '#332288'))  +
    theme(legend.position = 'bottom',
          axis.text.x = element_text(size = 10, angle = 90),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 10, face = 'bold'),
          axis.title.x = element_text(margin = margin(t = 10)),
          axis.title.y = element_text(margin = margin(t = 10)),
          legend.text = element_text(size = 9),
          legend.title = element_text(size = 10, face = 'bold'),
          plot.margin = margin(1,0,0,0, 'cm')) +
    guides(fill = guide_legend(nrow = 2)))


(fig3 <- plot_grid(fig3, labels = c('Habitat use including roost site'), hjust = 0))

ggsave(fig3, path = 'figs/', filename = 'FIG3.png', width = 6, height = 8, units = 'in', dpi = 300)
### 8 == PIGEON ROOSTING BEHAVIOUR ####################################
## Describe pigeon roosting behaviour.

## FIRST: identify GPS tracks present at a given roost site.
where_roost <- gps %>% 
  filter(inroost == 'Y') %>% # keep in roost observations
  filter(hour(datetime) > 22 | hour(datetime) < 4) %>% # in night time.
  group_by(date, name, roost_site) %>% # for each combination of date, name and roost time
  slice(1) %>% # keep only one
  dplyr::select(name, date, roost_site)

### 9 == FIGURE 4 -- PIGEON ROOSTING BEHAVIOUR ####################################
## SECOND: make a map to help identify roosting sites.
roost_sites <- st_read('bin/roost_sites.gpkg')
(fig4A <- ggmap(clachan_sat) + xlab('Longitude (°)') + ylab('Latitude (°)') +
    geom_sf(data = roost_sites %>% st_buffer(40) %>% st_transform(crs = st_crs(4326)), aes(fill = roost_site), inherit.aes = F) +
    scale_y_continuous(limits = c(57.65,57.68408)) +
    scale_fill_discrete(name = 'Roost site ID', labels = c('(1) Derelict building', '(2) - Barn',
                                                           '(3) Derelict building', '(4) Cattle byre')) +
    coord_sf() +
    theme(axis.title = element_blank(),
      legend.position = 'top') + 
    guides(fill = guide_legend(nrow = 2)))

## THIRD: make a plot of use of each roost site by pigeons.
(fig4B <- ggplot(data = where_roost, aes(y = as.character(roost_site), x = date)) + 
    geom_point() + facet_wrap(~name) + theme_bw() +
    ylab('Roost sites') + xlab('Date') + 
    theme(panel.spacing = unit(0.5, 'cm')))

(fig4 <- plot_grid(fig4A, fig4B, nrow = 2, labels = c('A)', 'B)')))
ggsave(fig4, path = 'figs/', filename = 'FIG4.png', dpi = 300)



### 9 == DISTANCE FROM ROOST ACROSS DAYS ####################################
## This section calculates distances from roost; including same roost, total morning roost, total night roost.
### FIRST: same roost calculations; distance from roost for birds that have used the same roost in the morning and at night. 
## Identify the gps fixes on days on which those conditions are met. 
same_roost <- gps %>% 
  dplyr::select(name, datetime, date, inroost, roost_site) %>% # keep only useful columns
  arrange(datetime) %>% # arrange by time.
  group_by(name, date) %>% # for each name and date.
  slice(c(1,n())) %>% # keep top observation
  filter(inroost == 'Y') %>% 
  mutate(test = roost_site == lag(roost_site, default = NA)) %>%
  ungroup() %>%
  filter(test == TRUE) %>%
  dplyr::select(name,date,roost_site)

same_roost %>% st_drop_geometry() %>% group_by(name) %>% count() # how much data is there?

## Calculate the distance from roost; same roost situation
pigeon_names <- unique(same_roost$name) # select pigeon names.
roost_sites <- st_read('bin/roost_sites.gpkg') # read in the shapefile of roost sites
dist_to_roost_df <- data.frame() # data frame to store results.

for (i in 1:length(pigeon_names)) { # select one pigeon from the list of their names.
  pigeon <- pigeon_names[i] # save its name.
  same_roost_days <- same_roost[same_roost$name == pigeon,]$date # keep unique days that pigeons has been in the same roost.
  for (j in 1:length(same_roost_days)) { # then for each of those days.
    day <- same_roost_days[j] # take one day at a time
    roost_id <- same_roost[same_roost$date == day & same_roost$name == pigeon,]$roost_site[1] # identify which roost was used
    tmp1 <- gps %>% filter(date == day) %>% filter(name == pigeon) %>% # prepare a dataset of gps points
      dplyr::select(name, date, datetime, habitat, inroost, roost_site) # for that one pigeon.
    roost <- roost_sites %>% filter(roost_site == roost_id) # keep only the roost site it was using
    tmp2 <- data.frame(dist = st_distance(tmp1, roost)) # calculate the distance between that roost site and gps
    tmp1 <- cbind(tmp1, tmp2 %>% st_drop_geometry()) # bind them together.
    dist_to_roost_df <- rbind(dist_to_roost_df, tmp1) # and bind to the dataframe of interest.
  }
}

# Now need to merge it back with date data to know where they are at times of day.
dist_time <- date_df %>%
  left_join(dist_to_roost_df, by = character()) %>% # cross-join
  filter(datetime %within% int) %>% # keep only observations within intervals.
  dplyr::select(name, datetime, id, habitat, inroost, dist)

# Add hours to it.
dist_time$id <- as.numeric(dist_time$id)
hours <- date_df %>% dplyr::select(start, id) %>%
  mutate(hour = paste(format(as.POSIXct(start), format = "%H:%M:%S"))) %>%
  dplyr::select(hour, id) %>%
  unique()
dist_time <- left_join(dist_time, hours, by = 'id')

# Average the distances
avg_dist_to_roost <- dist_time %>%
  group_by(hour, name) %>%
  summarise(mean_dist = mean(dist), sd_dist = sd(dist)) %>%
  mutate(mean_dist = drop_units(mean_dist), lower = mean_dist-sd_dist, upper = mean_dist+sd_dist) %>%
  group_by(hour) %>%
  mutate(total_mean_dist = mean(mean_dist))

### 10 == FIGURE 5 - Average distance from night roost (m); same roost on morning and evening. ####################################
(fig5 <- ggplot(data = avg_dist_to_roost, aes(x = as.POSIXct(hour, format = "%H:%M:%S", tz = 'GMT'), y = mean_dist, col = name)) + 
   geom_line(alpha = 0.7, linetype = 'dashed') + geom_point(alpha = 0.7, size = 1) +
   geom_line(data = avg_dist_to_roost,
             aes(x = as.POSIXct(hour, format = "%H:%M:%S", tz = 'GMT'), y = total_mean_dist), col = 'black', linewidth = 1.5, inherit.aes = F) +
   scale_x_datetime(labels = scales::date_format("%H:%M:%S", tz = 'GMT', locale = '')) +
   scale_color_discrete(name = 'Dove ID') + theme_bw() +
   xlab('Time of day (half-hour intervals)') + ylab('Mean distance from night roost (m)') +
   theme(legend.position = 'bottom'))

ggsave(fig5, path = 'figs/', filename = 'FIG5.png', dpi = 300)


#### 11 == PIGEON ASSOCIATIONS ####################################
## Explore whether pigeons associate with individuals from the same roost sites
## We need 15 minute windows for this analysis so back to the crude hour methods.
# Create a sequence of 15min slots slots; this is some of the crudest code ever but works.
full_h1 <- as.data.frame(with(expand.grid(date = as.Date(unique(gps$date)), hour = paste0(0:23, ":00:00")), 
                             sort(as.POSIXct(paste(date, hour)))))  # create a dataset of full hours
colnames(full_h1) <- 'start' # call it start as it will be the first value
half_h1 <- as.data.frame(with(expand.grid(date = as.Date(unique(gps$date)), hour = paste0(0:23, ":14:59")), 
                             sort(as.POSIXct(paste(date, hour))))) # create a dataset of half-hours meant for end.
colnames(half_h1) <- 'end' # call it end.
date_df1 <- cbind(full_h1, half_h1) # bind them together; gives you every XX:00:00 - XX:29:59 period.
full_h2 <- as.data.frame(with(expand.grid(date = as.Date(unique(gps$date)), hour = paste0(0:23, ":29:59")), 
                              sort(as.POSIXct(paste(date, hour))))) # repeat for ends of fulls hours.
colnames(full_h2) <- 'end'
half_h2 <- as.data.frame(with(expand.grid(date = as.Date(unique(gps$date)), hour = paste0(0:23, ":15:00")), 
                              sort(as.POSIXct(paste(date, hour))))) # repeat for starts of half hours.
colnames(half_h2) <- 'start'
date_df2 <- cbind(half_h2, full_h2) # bind them together; gives you every XX:30:00 - XX:59:59 period.


full_h3 <- as.data.frame(with(expand.grid(date = as.Date(unique(gps$date)), hour = paste0(0:23, ":30:00")), 
                              sort(as.POSIXct(paste(date, hour)))))  # create a dataset of full hours
colnames(full_h3) <- 'start' # call it start as it will be the first value
half_h3 <- as.data.frame(with(expand.grid(date = as.Date(unique(gps$date)), hour = paste0(0:23, ":44:59")), 
                              sort(as.POSIXct(paste(date, hour))))) # create a dataset of half-hours meant for end.
colnames(half_h3) <- 'end' # call it end.
date_df3 <- cbind(full_h3, half_h3) # bind them together; gives you every XX:00:00 - XX:29:59 period.
full_h4 <- as.data.frame(with(expand.grid(date = as.Date(unique(gps$date)), hour = paste0(0:23, ":59:59")), 
                              sort(as.POSIXct(paste(date, hour))))) # repeat for ends of fulls hours.
colnames(full_h4) <- 'end'
half_h4 <- as.data.frame(with(expand.grid(date = as.Date(unique(gps$date)), hour = paste0(0:23, ":45:00")), 
                              sort(as.POSIXct(paste(date, hour))))) # repeat for starts of half hours.
colnames(half_h4) <- 'start'
date_df4 <- cbind(half_h4, full_h4) # bind them together; gives you every XX:30:00 - XX:59:59 period.

date_15_df <- rbind(date_df1, date_df2, date_df3, date_df4) %>%
  mutate(int = interval(start, end)) %>% # create a time interval column.
  arrange(start) %>%
  mutate(id = rep(seq(1, 96, 1), 96, 2784))
rm(date_df1, date_df2, date_df3, date_df4, full_h1, full_h2, full_h3, full_h4, half_h1, half_h2, half_h3, half_h4)

### How many gps fixes in each 15 minute window?
gps_fix_15 <- date_15_df %>%
  left_join(gps %>% filter(inroost == 'N'), by = character()) %>% # cross-join while removing in-roost observations.
  filter(datetime %within% int) %>% st_as_sf()

## Create a new 'hours object' for plotting.
hours_15 <- date_15_df %>% dplyr::select(start, id) %>%
  mutate(hour = paste(format(as.POSIXct(start), format = "%H:%M:%S"))) %>%
  dplyr::select(hour, id) %>%
  unique()

## Number of fixes per time interval and date of each pigeon.
gps_fix_15_summ <- gps_fix_15 %>%
  group_by(id, name, date) %>%
  count() %>%
  left_join(hours_15, by = 'id')

## Exploratory plot of the fixes above.
ggplot(data = gps_fix_15_summ, aes(x = as.POSIXct(hour, format = "%H:%M:%S"), y = n, fill = name)) + geom_bar(stat = 'identity') +
  scale_x_datetime(labels = scales::date_format("%H:%M:%S")) + facet_wrap(~date)

### How have pigeons associated within those 15 minute intervals?
asso <- as.data.frame(st_intersects(gps_fix_15 %>% # Check where does gps_fix_15 intersect with itself
                                      st_buffer(dist = 20), # with a buffer of 20 units (here meters)
                                    gps_fix_15 %>% st_buffer(dist = 20)))
t1 <- gps_fix_15 %>% # create object t1 to translate row.id into observations
  dplyr::select(name, datetime, int, id, habitat) %>%
  mutate(row.id = row_number())
asso <- left_join(asso, t1, by = 'row.id') %>% dplyr::select(-row.id) %>% rename(row.id = col.id) # translate row.id into intervals
t2 <- gps_fix_15 %>% 
  dplyr::select(name, datetime) %>% # Identify the intersections.
  rename(other = name,
         datetime2 = datetime,
         geom2 = geom) %>% mutate(row.id = row_number())
asso <- left_join(asso, t2, by = 'row.id') %>% filter(datetime2 %within% int) %>%
  mutate(other = ifelse(other == name, NA, other)) %>% drop_na() 

## Summary of the association dataset.
asso_summ <- asso %>% st_drop_geometry() %>% dplyr::select(name, other, int) %>% distinct() %>%
  group_by(name, other) %>% count()

### Add total obs per pigeon
total_obs <- gps_fix_15 %>% st_drop_geometry() %>% dplyr::select(name, int) %>% distinct() %>%
  group_by(name) %>% count() %>% rename(total_obs = n)

asso_totals <- asso_summ %>% group_by(name) %>% summarise(sum = sum(n))
total_obs <- total_obs %>% left_join(asso_totals, by = 'name') %>% mutate(alone = total_obs - sum) %>%
  dplyr::select(name, alone) %>% pivot_longer(alone, names_to = 'other', values_to = 'n')

asso_summ <- rbind(asso_summ, total_obs)

#### 12 == FIGURE S7 - Pigeon associations ####################################

(fig4 <- ggplot(data = asso_summ, aes(x = name, y = n, fill = other)) +
    geom_bar(position="fill", stat="identity") +
    xlab('Dove ID') + ylab('Proportion of time intervals') + theme_bw() + 
    scale_fill_manual(name = 'Associating with',values = c('alone' = 'grey', 'EA49571' = '#E769F1', 'EA49570' = '#3CABD3', 'EA49568' = '#39A783', 'EA49503' = '#98993C', 'EA49502' = '#f8756d'),
                      labels = c('Not-associated', 'EA49502', 'EA49503', 'EA49568', 'EA49570', 'EA49571')) +
    ggtitle('Proportion of 15-minute time intervals \nby association between tagged Rock Doves') +
    theme(legend.position = 'bottom'))
ggsave(fig4, path = 'figs/', filename = 'FIG4.png', width = 5, unit = 'in',  dpi = 300)


### Number of 15 min intervals with more than 2 pigeons anywhere mins anywhere.
f <- gps_fix_15 %>% st_drop_geometry %>% 
  group_by(int) %>% count() %>% filter(n >= 2) ## the number of rows here is the answer.
nrow(f)
### Number of intervals with pigeons being within 20m of each other in those 15 minute interval.
g <- asso %>% st_drop_geometry() %>% 
  group_by(int) %>% count() %>% filter(n >= 2) ## number of rows here is the answer.
nrow(g)

