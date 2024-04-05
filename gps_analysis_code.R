##### GPS ANALYSIS CODE #####
### To accompany the manuscript: Use of anthropogenic landscapes in a wild Rock Dove Columba livia population
### William J. Smith, Stephen Portugal, Michał T. Jezierski.
## Versioned by MTJ 4th April 2024.

### Working directory and necessary packages. ####
library(tidyverse)
library(sf)
library(lubridate)
library(ggmap)
library(cowplot)
library(adehabitatHR)
library(cowplot)
library(units)
sf_use_s2(F)

setwd('~/Desktop/rock doves/gps_folder/')


#### LOAD IN DATASETS PREPARED IN gps_data_code.R ====
## TO WILL: The directory in the shared folder is 'bin/'.
gps <- st_read('bin/filtered_gps_habitat_dataset.gpkg') %>%
  mutate(date = as.Date((format(as.POSIXct(datetime), "%Y-%m-%d")))) # the filtered gps dataset.
clachan_base <- st_read('bin/clachan_base.gpkg') # base grey map for Clachan.
clachan_habitat <- st_read('bin/clachan_new_habitat.gpkg') # habitat shapefile.
load('sat_map_clachan.RData') # Satellite map, all three lines below.
clachan_sat <- map
rm(map)

#### RESULTS SECTION 1 - GENERAL STATISTICS REGARDING POINTS ====

## How many tracking days do we have across pigeons?
gps %>% st_drop_geometry() %>% # remove geometry
  group_by(name) %>% # across all pigeon dates.
  summarise(start_date = min(date), end_date = max(date)) %>% # summarise the minimum and maximum date.
  mutate(date_range = end_date - start_date) %>% # calculate the date range.
  print() %>% # print it for each pigeon name.
  summarise(mean = mean(date_range), sd = sd(date_range)) # and across all name, give mean and SD.

## Estimate home range
# Using adehabitatHR core functions.
spatial_gps <- gps %>% dplyr::select(name) %>% as_Spatial() # turn the observations into spatial (other data format)
spatial_gps$name <- as.factor(spatial_gps$name) # make names as factors as required for the analysis
kud <- kernelUD(spatial_gps, h="href") # estimates of Kernel Home-Range
(ver_95 <- getverticeshr(kud, 95)) # extract the 95% contour of home range
(ver_50 <- getverticeshr(kud, 50)) # extract the 50% contour of home range
(mcp_50 <- mcp(spatial_gps, 50))
(mcp_95 <- mcp(spatial_gps, 95))


## Summary of pigeon tracking data 
options(pillar.sigfig = 2) # forces significant figures.
gps %>% st_drop_geometry() %>% # drop geometry
  group_by(name) %>% # across all pigeon names
  summarise(start_date = min(date), end_date = max(date), fix = n()) %>% # recreate the dates above + no. of fixes.
  mutate(date_range = end_date - start_date) %>% 
  mutate(sex = c('Male', 'Male', 'Male', 'Female (?)', 'Female'), # add bonus information
         age = c('Adult', 'Adult', 'Adult', 'Adult', 'Adult'),
         home_50_ha = num(c(5.485856, 34.322273, 14.104215, 78.748625, 6.729580), digits = 3), # add data from HR estimation.
         home_95_ha = num(c(68.37009, 134.09548, 62.03079, 308.29462, 38.66206), digits = 3)) %>%
  print() %>%
  summarise(mean_50 = mean(home_50_ha), sd_50 = sd(home_50_ha), mean_95 = mean(home_95_ha), sd_95 = sd(home_95_ha))

## Compare pigeon home ranges but focusing on convex polygons.
# What are the means and sds of the data?
mean(mcp_50$area)
sd(mcp_50$area)
mean(mcp_95$area)
sd(mcp_95$area)

# Data from other papers
rose_data <- c(2.9, 2.9, 3.2, 5.1, 5.2,12.5, 14.8, 16.3, 24.2, 46.7, 144.1, 150.6)
mean(rose_data)
sd(rose_data)
wilcox.test(rose_data, mcp_95$area)

# Get populations level mcp as well?


##### FIGURE 1 - OVERALL MOVEMENT OF PIGEONS THROUGH LANDSCAPE.
gps_wgs_coords <- gps %>% # create a dataset with X/Y rather than geometry coordinates.
  st_transform(crs = st_crs(4326)) %>% # must recreate coordinates for lon/lat projections
  st_coordinates() %>% # # extract coordinates
  as.data.frame() %>% # as data frame
  mutate(name = gps$name, datetime = gps$datetime) # add name and datetime columns

inset_nu <- st_read('bin/GBR_adm0.shp') %>% # shapefile from https://www.diva-gis.org/gdata
  st_crop(xmin = -7.577820, ymin = 57.464158, xmax = -7.009277, ymax = 57.747412) # crop to NUist box

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

(fig1A <- ggmap(clachan_sat) + 
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

(fig1 <- ggdraw() + # code joining the two figures
  draw_plot(fig1A) +
  draw_plot(inset_nu_map, x = 0.79, y = 0.1, width = 0.25, height = 0.25))

ggsave(fig1, path = 'figs/', filename = 'FIG1.png', dpi = 300, width = 7, height = 4.5, units = 'in')

##### FIGURE S1-S5 - MOVEMENT OF EACH PIGEON AND HOME RANGE.
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
(figS1A <- ggmap(clachan_sat) + 
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
(figS1B <- ggmap(clachan_sat) + 
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

(figS1 <- plot_grid(figS1A, figS1B, labels = c('A)', 'B)'), label_fontfamily = 'serif', label_size = 10))
ggsave(figS1, path = 'figs/', filename = 'FIGS1.png', dpi = 300, width = 9, height = 4.5, units = 'in')

#### FIGURE S2 - EA49503.

(figS2A <- ggmap(clachan_sat) + 
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
(figS2B <- ggmap(clachan_sat) + 
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

(figS2 <- plot_grid(figS2A, figS2B, labels = c('A)', 'B)'), label_fontfamily = 'serif', label_size = 10))
ggsave(figS2, path = 'figs/', filename = 'FIGS2.png', dpi = 300, width = 9, height = 4.5, units = 'in')

#### FIGURE S3 - EA49568.

(figS3A <- ggmap(clachan_sat) + 
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
(figS3B <- ggmap(clachan_sat) + 
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

(figS3 <- plot_grid(figS3A, figS3B, labels = c('A)', 'B)'), label_fontfamily = 'serif', label_size = 10))
ggsave(figS3, path = 'figs/', filename = 'FIGS3.png', dpi = 300, width = 9, height = 4.5, units = 'in')

#### FIGURE S4 - EA49570.

(figS4A <- ggmap(clachan_sat) + 
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
(figS4B <- ggmap(clachan_sat) + 
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

(figS4 <- plot_grid(figS4A, figS4B, labels = c('A)', 'B)'), label_fontfamily = 'serif', label_size = 10))
ggsave(figS4, path = 'figs/', filename = 'FIGS4.png', dpi = 300, width = 9, height = 4.5, units = 'in')

#### FIGURE S5 - EA49571.

(figS5A <- ggmap(clachan_sat) + 
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
(figS5B <- ggmap(clachan_sat) + 
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

(figS5 <- plot_grid(figS5A, figS5B, labels = c('A)', 'B)'), label_fontfamily = 'serif', label_size = 10))
ggsave(figS5, path = 'figs/', filename = 'FIGS5.png', dpi = 300, width = 9, height = 4.5, units = 'in')

#### RESULTS SECTION 2 - TIME SPENT ACROSS HABITATS ====
## Approach taken here: obtain a dataset of 0.5h periods across all monitoring dates; then check which
## GPS points happened during which period; for each pigeon, keep the one it spent the most time in.
# Create a sequence of half an hour slots; this is some of the crudest code ever but works.
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

# Time spent in habitat; NB: this is a clever pipe.
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


##### FIGURE 2 - Time spent in each habitat 
## Subplot A, excluding roost sites (also only day hours).
(fig2a <- ggplot(data = time_in_habitat, aes(x = as.POSIXct(hour, format = "%H:%M:%S"), y = n, fill = habitat)) + 
  geom_bar(position = 'fill', stat = 'identity') + theme_bw() +
  xlab('Time of day (half-hour intervals)') + ylab('Proportion of time in habitat') + # note I cut some blocks out below
  scale_x_datetime(labels = scales::date_format("%H:%M:%S"), limits = c(as.POSIXct('08:30:00', format = "%H:%M:%S"), as.POSIXct('20:00:00', format = "%H:%M:%S"))) +
    scale_fill_manual(name = '', values = c("Arable land" = "#DDCC77", "Atlantic pasture" = '#CC6677',
                                 "Machair" = '#117733', "Dunes" = '#6699CC', 'Unknown or other' = '#332288')) +
  theme(legend.position = 'none',
        axis.text.x = element_text(size = 10, angle = 90),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 10, face = 'bold'),
        axis.title.x = element_blank(),
        axis.title.y = element_text(margin = margin(t = 10)),
        plot.margin = margin(1,0,0,0, 'cm')) +
    guides(fill = guide_legend(nrow = 2)))

## Subplot B, including roost sites (all day).
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
(fig2b <- ggplot(data = time_in_habitat_roost, aes(x = as.POSIXct(hour, format = "%H:%M:%S"), y = n, fill = habitat)) + 
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


(fig2 <- plot_grid(fig2a, fig2b, rel_heights = c(0.4, 0.6), nrow = 2, labels = c('A) Habitat use without roost site', 
                                                                                  'B) Habitat use including roost site'), hjust = 0))

ggsave(fig2, path = 'figs/', filename = 'FIG2.png', width = 6, height = 8, units = 'in', dpi = 300)

#### RESULTS SECTION 3 - DISTANCE FROM ROOST SITE ACROSS DAYS ====
## Will wants birds that have roosted in the same site on those days.
same_roost <- gps %>% 
  dplyr::select(name, datetime, date, inroost, roost_site) %>% # keep only useful columns
  arrange(datetime) %>% # arrange by time.
  group_by(name, date) %>% # for each name and date.
  slice(c(1,n())) %>% # keep top observation
  filter(inroost == 'Y') %>% 
  mutate(test = roost_site == lag(roost_site, default = NA)) %>%
  ungroup() %>%
  filter(test == TRUE) %>% # I will be honest; i forgot what each step does.
  dplyr::select(name,date,roost_site)

same_roost %>% st_drop_geometry() %>% group_by(name) %>% count() # how much data is there?

# Also interesting - not same roost.
## Code as above but keeping the birds that went to a different place.
not_same_roost <- gps %>% dplyr::select(name, datetime, date, inroost, roost_site) %>%
  arrange(datetime) %>% 
  group_by(name, date) %>% 
  slice(c(1,n())) %>%
  filter(inroost == 'Y') %>% 
  mutate(test = roost_site == lag(roost_site, default = NA)) %>%
  ungroup() %>%
  filter(test == FALSE) %>% 
  dplyr::select(name,date)

not_same_roost %>% st_drop_geometry() %>% group_by(name) %>% count() 

# Calculate distance from roost.
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

(fig3 <- ggplot(data = avg_dist_to_roost %>% filter(name != 'EA49570'), aes(x = as.POSIXct(hour, format = "%H:%M:%S"), y = mean_dist, col = name)) + 
  geom_line() + geom_point() +
  geom_line(data = avg_dist_to_roost %>% filter(name != 'EA49570'),
            aes(x = as.POSIXct(hour, format = "%H:%M:%S"), y = total_mean_dist), col = 'black', linewidth = 1.5, alpha = 0.5, inherit.aes = F) +
  scale_x_datetime(labels = scales::date_format("%H:%M:%S")) +
  scale_color_discrete(name = 'Dove ID') + theme_bw() +
  xlab('Time of day (half-hour intervals)') + ylab('Mean distance from night roost (m)') +
  theme(legend.position = 'bottom'))

ggsave(fig3, path = 'figs/', filename = 'FIG3.png', dpi = 300)
  
##### PIGEON ROOSTING BEHAVIOURS ====
### Where do pigeons spend the night?
where_roost <- gps %>% 
  filter(inroost == 'Y') %>% # keep in roost observations
  filter(hour(datetime) > 22 | hour(datetime) < 4) %>% # in night time.
  group_by(date, name, roost_site) %>% # for each combination of date, name and roost time
  slice(1) %>% # keep only one
  dplyr::select(name, date, roost_site)

### Make a map to help locate these.
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

## Make the plot of use of roosting sites
(fig4B <- ggplot(data = where_roost, aes(y = as.character(roost_site), x = date)) + 
    geom_point() + facet_wrap(~name) + theme_bw() +
    ylab('Roost sites') + xlab('Date') + 
    theme(panel.spacing = unit(0.5, 'cm')))

(fig4 <- plot_grid(fig4A, fig4B, nrow = 2, labels = c('A)', 'B)')))
ggsave(fig4, path = 'figs/', filename = 'FIG4.png', dpi = 300)

#### PIGEON ASSOCIATIONS =====
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

### FIGURE 5 - how do pigeons associate with each other?

(fig5alt <- ggplot(data = asso_summ, aes(x = name, y = n, fill = other)) +
    geom_bar(position="fill", stat="identity") +
    xlab('Dove ID') + ylab('Proportion of time intervals') + theme_bw() + 
    scale_fill_manual(name = 'Associating with',values = c('alone' = 'grey', 'EA49571' = '#E769F1', 'EA49570' = '#3CABD3', 'EA49568' = '#39A783', 'EA49503' = '#98993C', 'EA49502' = '#f8756d'),
                      labels = c('Not-associated', 'EA49502', 'EA49503', 'EA49568', 'EA49570', 'EA49571')) +
    ggtitle('Proportion of 15-minute time intervals \nby association between tagged Rock Doves') +
    theme(legend.position = 'bottom'))
ggsave(fig5alt, path = 'figs/', filename = 'FIG5.png', width = 5, unit = 'in',  dpi = 300)


(fig5 <- ggplot(data = asso, aes(x = name, y = other)) + geom_count() +
  xlab('Dove ID') + ylab('Dove ID') + theme_bw() + 
    ggtitle('No. of fixes where two doves were\nwithin 20m of each other a in 15-minutes interval') +
  theme(legend.position = 'bottom'))

ggsave(fig5, path = 'figs/', filename = 'FIG56.png', dpi = 300)
  

### Number of 15 min intervals with more than 2 pigeons anywhere mins anywhere.
f <- gps_fix_15 %>% st_drop_geometry %>% 
  group_by(int) %>% count() %>% filter(n >= 2) ## the number of rows here is the answer.
nrow(f)
### Number of intervals with pigeons being within 20m of each other in those 15 minute interval.
g <- asso %>% st_drop_geometry() %>% 
  group_by(int) %>% count() %>% filter(n >= 2) ## number of rows here is the answer.
nrow(g)

################################ DANGER ZONE CODE BELOW IS OBSOLETE BUY MAY BE USEFUL ====================

pigeon_df <- gps_fix_15 %>% filter(name == pigeon) %>% dplyr::select(name, datetime, int, id, habitat) %>%
  mutate()


for (i in 1:length(names)) {
  pigeon <- names[i]
  pigeon_ints <- gps_fix_15[gps_fix_15$name == pigeon,]$int
  for (j in 1:length(pigeon_ints)) {
    int1 <- pigeon_ints[i]
    tmp1 <- gps_fix_15 %>% filter(name == pigeon & datetime %within% int1) %>% 
      dplyr::select(name, datetime, id, habitat) %>% st_buffer(dist = 20)
    tmp2 <- gps_fix_15 %>% filter(name != pigeon & datetime %within% int1) %>%
      dplyr::select(name, datetime, id, habitat) %>% st_buffer(dist = 20)
    st1 <- st_intersection(tmp1, tmp2) %>% rename('ot_name' = 'name.1', 'ot_datetime' = 'datetime.1') %>%
      dplyr::select(name, datetime, id, habitat, ot_name, ot_datetime)
    asso_fix_15 <- rbind(asso_fix_15, st1)
  }
}











asso <- gps_int %>% left_join(hours, by = 'id')

t <- asso %>% group_by(hour, name) %>% count()
ggplot(data = t, aes(x = as.POSIXct(hour, format = "%H:%M:%S"), y = n, col = name)) + geom_line() +
  scale_x_datetime(labels = scales::date_format("%H:%M:%S"))

### Calculate distance from roost.
roost_dist <- as.data.frame(st_distance(gps_hab, roost_sites))
roost_dist <- drop_units(roost_dist)
colnames(roost_dist) <- c('roost1', 'roost2', 'roost3', 'roost4')
roost_dist_df <- cbind(gps_hab, roost_dist) %>% pivot_longer(roost1:roost4, names_to = 'roost_id', values_to = 'roost_dist')
library(units)
ggplot(data = roost_dist_df, aes(x = hour(time), y = roost_dist, col = roost_id)) + geom_line() + facet_wrap(~name)

### Distance from morning roost.

day_roost <- gps %>%
  filter(inroost == "Y", hour(datetime) < 8) %>%
  group_by(date, name) %>%
  arrange(desc(datetime)) %>%
  slice(1) %>%
  ungroup() %>%
  select(name, date, time, datetime, roost_site)

test_df <- gps_hab 
date_df <- test_df %>% select(date) %>% st_drop_geometry() %>% distinct()
name_vec <- gps_hab %>% st_drop_geometry() %>% select(name) %>% distinct()

dist_df <- data.frame()

for (i in 1:length(date_df$date)) {
  for (j in 1:length(name_vec$name)) {
    dat <- date_df$date[i]
    nam <- name_vec$name[j]
    postroost <- test_df %>% filter(date == dat, name == nam)
    roost <- day_roost %>% filter(date == dat, name == nam) %>% select(roost_site) %>% st_drop_geometry() 
    if (nrow(roost) == 0 || nrow(postroost) == 0) {
      next
    }
    site <- roost_sites %>% filter(roost_site == roost$roost_site)
    dist_to_roost <- as.data.frame(st_distance(postroost, site))
    tmp <- data.frame(name = nam,
                      date = dat,
                      time = postroost$time,
                      dist = dist_to_roost) %>% distinct()
    
    colnames(tmp) <- c('name', 'date', 'time', 'dist')
    dist_df <- rbind(dist_df, tmp)
  }
}

dist_df$dist <- drop_units(dist_df$dist)
dist_df_avg <- dist_df %>% mutate(hour = hour(time)) %>% group_by(hour, name) %>%
  summarise(mean_dist = mean(dist))
dist_df_avg_all <- dist_df %>% mutate(hour = hour(time)) %>% group_by(hour) %>% summarise(mean_dist = mean(dist))

ggplot() + geom_line(data = dist_df_avg, aes(x = hour, y = mean_dist, col = name)) +
  geom_point(data = dist_df_avg, aes(x = hour, y = mean_dist, col = name)) +
  geom_line(data = dist_df_avg_all, aes(x = hour, y = mean_dist), col = 'black', size = 2) +
  geom_point(data = dist_df_avg_all, aes(x = hour, y = mean_dist), col = 'black', size = 2)

time_in_habitat_roost$id <- as.numeric(time_in_habitat_roost$id)
hours <- date_df %>% select(start, id) %>%
  mutate(hour = paste(format(as.POSIXct(start), format = "%H:%M:%S"))) %>%
  select(hour, id) %>%
  unique()
time_in_habitat_roost <- left_join(time_in_habitat_roost, hours, by = 'id')


# Calculate proportion of time spent in a given habitat category, with no roost sites.
time_spent_noroost <- gps %>% filter(inroost == 'N') %>% 
  mutate(hour = hour(datetime)) %>% arrange(name, datetime)

hours <- unique(time_spent_noroost$hour)
tmp1 <- time_spent_noroost[1:100,] %>% select(name, habitat, datetime) %>% mutate(hour = hour(datetime)) 
td <- tmp1 %>% mutate(diff = datetime - lag(datetime)) %>% replace(is.na(.), 0)
complete_hours <- expand.grid(datetime = seq(min(tmp1$datetime), max(tmp1$datetime), by = "hour"))

# Merge the complete sequence with your data frame
merged_data <- complete_hours %>%
  left_join(tmp1, by = "datetime") %>%
  arrange(datetime)

# Fill in missing values in 'name' and 'habitat' columns
merged_data <- merged_data %>%
  fill(name, habitat, .direction = "down")

# Fill in missing values in 'hour' column
merged_data$hour <- ifelse(is.na(merged_data$hour), hour(merged_data$datetime), merged_data$hour)

# Calculate the time difference ('diff') column
merged_data <- merged_data %>%
  mutate(diff = difftime(lead(datetime), datetime, units = "mins"))

# Remove the last row as it will have NA values
merged_data <- merged_data[1:(nrow(merged_data) - 1), ]


std <- td %>% group_by(habitat) %>% summarise(sum_diff = sum(diff))

unique(tmp1$habitat)

mutate(habitat = case_when(
  habitat == "H21A0 - Machairs" ~ 1,
  habitat == "H2130 - Fixed dunes" ~ 2,
  habitat == "H2120 - Shifting dunes" ~ 3,
  habitat == "Arable land" ~ 4,
  habitat == "Atlantic Cynosurus-Centaurea pastures" ~ 5
))


  
tmp1 %>% select(name, hour, datetime, habitat) %>% st_drop_geometry() %>% print(1:20) 
  
  df %>%
  arrange(category, randtime) %>%
  group_by(category) %>%
  mutate(diff = randtime - lag(randtime),
         diff_secs = as.numeric(diff, units = 'secs'))

test <- gps %>% mutate(hour = hour(datetime))

rmhours <- gps %>% mutate(hour = hour(datetime)) %>% filter(inroost == 'N') %>%
  group_by(hour) %>% summarise(n = n()) %>% dplyr::select(hour)
hours <- as.integer(hours$hour)


ggplot(data = time_spent_noroost, aes(x = hour, y = time_spent, fill = habitat)) + 
  geom_bar(stat = 'identity', position = 'fill') + theme_bw() +
  scale_fill_discrete(labels = c('Improved grassland', 'Arable land', expression('Atlantic'~italic(Cynosurus-Centaurea)~'pastures'),
                                 'Artificial sites', 'Flood sward', 'Shifting dunes', 'Fixed dunes', 'Humid dune slacks',
                                 'Machair', 'Transition mires', 'Unvegetated beaches')) +
         theme(legend.position = 'bottom') 


#### RESULTS SECTION 1 - GENERAL STATISTICS REGARDING POINTS ====

## How many tracking days do we have across pigeons?
gps %>% st_drop_geometry() %>% group_by(name) %>%
  summarize(start_date = min(date), end_date = max(date)) %>%
  mutate(date_range = end_date - start_date) %>%
  print() %>%
  summarise(mean = mean(date_range), sd = sd(date_range))





### Will wants birds that go to same site morng and evening


### Return to any roost site

ggplot(data = gps_hab, aes(x = hour(datetime), fill = name)) + geom_bar()

test_any_roost <- gps_hab %>% select(name, date, time, datetime) %>%
  st_distance(roost_sites) %>% as.data.frame() %>%
  cbind(gps_hab %>% select(name, date, time, datetime) %>% st_drop_geometry()) %>%
  pivot_longer(V1:V4, names_to = 'roost_site', values_to = 'dist') %>%
  mutate(dist = drop_units(dist)) %>% mutate(hour = hour(time)) %>%
  group_by(name, roost_site, hour) %>% summarise(mean_dist = mean(dist))

ggplot(data = test_any_roost, aes(x = hour, y = mean_dist, col = name)) + geom_point() + geom_line() + facet_wrap(~roost_site)

# Plot the bar graph
ggplot(test_prob, aes(x = name, y = total_time_spent, fill = HABITAT_NA)) +
  geom_bar(stat = "identity", position = 'fill') +
  labs(title = "Proportion of Time Spent in Each Habitat",
       x = "Habitat",
       y = "Proportion of Time Spent") +
  theme_minimal()

ggplot(data = test_prob, aes(x = HABITAT_NA, fill = name)) + geom_bar() +
  theme(axis.text.x = element_text(angle = 90))

tm_shape(clachan) + tm_fill() + tm_shape(gps) + tm_dots() 

## How much data do we have per individual?
data_count <- gps %>% group_by(date, name) %>% summarise(n = n()) %>% st_drop_geometry()
ggplot(data = data_count, aes(x = date, y = n, col = name)) + geom_line()

### Map of of behaviour per bird
library(rnaturalearth)
north_uist <- ne_download(scale = 10, type = 'land', category = 'physical') %>% st_as_sf(crs = st_crs(4326)) %>%
  st_crop(xmin = -7.580566, ymin = 57.498855, xmax = -7.014771, ymax = 57.724686)

### Need to classify roost sites etc.