####### ANALYSIS CODE =======================================================================================
### Code contains analytical steps taken on data from GPS-tracking of Rock Doves in Jul/Aug 2023. 
### Manuscript: Use of anthropogenic landscapes in a wild Rock Dove Columba livia population
### In review in: Ornithology
### Authors: William J. Smith, Stephen J. Portugal, Michał T. Jezierski.
### Versioned by MTJ 12th September 2024.

### 1 == SET UP THE DIRECTORY, LOAD PACKAGES AND DATASETS REQUIRED. ####################################
#### SET WORKING DIRECTORY. ##
setwd('~/Desktop/gps revision/code/')

#### PACKAGES NEEDED. ##
library(tidyverse)
library(sf)
library(lubridate)
library(ggmap)
library(cowplot)
library(adehabitatHR)
library(cowplot)
library(units)
library(ctmm)

sf_use_s2(F) # also turn of spherical geometry for sf.

#### LOAD IN THE DATASET. ##
## Data as prepared in the file 'gps_data_code.R'.
gps <- st_read('bin/filtered_gps_habitat_dataset.gpkg') %>% # the filtered dataset; see file above and METHODS.
  mutate(date = as.Date((format(as.POSIXct(datetime), "%Y-%m-%d")))) # reformat date so it can be used more easily.
clachan_base <- st_read('bin/clachan_base.gpkg') # base map of the study area.
clachan_habitat <- st_read('bin/clachan_new_habitat.gpkg') # habitat map of the study area.
load('bin/sat_map_clachan.RData') # loading satellite map of the study area; execute all three lines.
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

### And the total across the dataset.
gps_dist %>% ungroup() %>% summarise(mean_dist = mean(total), sd_dist = sd(total)) %>% print()

### THIRD: check the impact of 10/07 (day of tagging on this)

gps_dist_no_tag <- gps_dist %>% filter(date != as.Date('2023-07-10'))

wilcox.test(gps_dist$total, gps_dist_no_tag$total) # no difference for all Rock Doves

## Individual tests.
ea49502 <- gps_dist %>% filter(name == 'EA49502')
ea49502_nt <- gps_dist_no_tag %>% filter(name == 'EA49502')
wilcox.test(ea49502$total, ea49502_nt$total) 

ea49503 <- gps_dist %>% filter(name == 'EA49503')
ea49503_nt <- gps_dist_no_tag %>% filter(name == 'EA49503')
wilcox.test(ea49502$total, ea49502_nt$total) 

ea49568 <- gps_dist %>% filter(name == 'EA49568')
ea49568_nt <- gps_dist_no_tag %>% filter(name == 'EA49568')
wilcox.test(ea49502$total, ea49502_nt$total) 

ea49570 <- gps_dist %>% filter(name == 'EA49570')
ea49570_nt <- gps_dist_no_tag %>% filter(name == 'EA49570')
wilcox.test(ea49502$total, ea49502_nt$total) 

ea49571 <- gps_dist %>% filter(name == 'EA49571')
ea49571_nt <- gps_dist_no_tag %>% filter(name == 'EA49571')
wilcox.test(ea49502$total, ea49502_nt$total) 

### 3 == HOME RANGES ####################################
### FIRST: Estimate Rock Dove home ranges. 
## Using adehabitatHR and ctmm estimate home ranges of Rock Doves:
## 1) Using ctmm calculate Rock Dove home ranges for reporting in Table 1.
## 2) Using adehabitatHR calculate maximum convex polygon ranges for comparison with data on Feral Pigeons.

### 1) - Calculate AKDE ranges using ctmm.
## Load in the dataset as telemtry object (standardised using Movebank).
ctmm_gps <- as.telemetry('bin/Rock Dove, Jezierski&Smith, Outer Hebrides(1).csv')
names <- c('EA49502', 'EA49503', 'EA49568', 'EA49570', 'EA49571')

for (i in 1:5) {
  guess <- ctmm.guess(ctmm_gps[[i]], interactive = F)
  fit <- ctmm.fit(ctmm_gps[[i]], guess)
  ud <- akde(ctmm_gps[[i]], fit)
  print(names[i])
  print(summary(ud))
}

## Average AKDE range and standard deviation
mean(c(103.1283, 147.2843, 66.60396, 321.1663, 53.36784))
sd(c(103.1283, 147.2843, 66.60396, 321.1663, 53.36784))

### 2) - Calculate MCP ranges using adehabitatHR.
## Convert the GPS data into Spatial for use with adehabitat HR.
spatial_gps <- gps %>% dplyr::select(name) %>% as_Spatial() # turn the observations into spatial (other data format)
spatial_gps$name <- as.factor(spatial_gps$name) # make names as factors as required for the analysis.

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
         home_AKDE_est = num(c(103.1283, 147.2843, 66.60396, 321.1663, 53.36784), digits = 2),
         home_AKDE_min = num(c(90.68896, 118.6105, 51.17934, 229.7206, 46.40701), digits = 2),
         home_AKDE_max = num(c(116.3555, 179.0158, 84.02922, 427.6948, 60.80797), digits = 2),
         home_mcp_100_ha = num(c(218.68409, 129.18426, 34.62694, 206.07421, 56.46548), digits = 2)) %>%
  print() %>%
  summarise(mean_AKDE = mean(home_AKDE_est), sd_AKDE = sd(home_AKDE_est),
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

## THIRD: to examine up to modern standards, also examine Rock Dove ranges using AKDE methods from ctmm package.

### 4 == FIGURE 2 - MOVEMENT OF ROCK DOVES THROUGH LANDSCAPE ####################################
### GPS dataset needs to be fully compatible with WGS coordinate system.
gps_wgs_coords <- gps %>% # create a dataset with X/Y rather than geometry coordinates.
  st_transform(crs = st_crs(4326)) %>% # must recreate coordinates for lon/lat projections
  st_coordinates() %>% # # extract coordinates
  as.data.frame() %>% # as data frame
  mutate(name = gps$name, datetime = gps$datetime) # add name and datetime columns

### Load inset North Uist map.
inset_nu <- st_read('bin/UIST.gpkg') %>% # shapefile from https://www.diva-gis.org/gdata
  st_crop(xmin = -7.577820, ymin = 57.464158, xmax = -7.009277, ymax = 57.747412) # crop to NUist box

### Create the map of North Uist with study area.
(inset_nu_map <- ggplot() + 
    geom_sf(data = inset_nu) + # add sf object called inse nu
    geom_rect(aes(xmin = -7.280598, 
                  ymin = 57.65647,
                  xmax = -7.207275,
                  ymax = 57.68704), col = 'purple', fill = NA) + # rectangle around Clachan.
    theme_void() +
    theme(panel.grid = element_blank(),
          panel.background = element_rect(fill = 'white'),
          panel.border = element_rect(colour = 'black', linewidth = 1, fill = NA),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank()))

ggsave(inset_nu_map, path = 'figs/FIG1/', filename = 'NU.png', width = 1, height = 1, units = 'in')

### Create the map of the UK with North Uist in it.
inset_uk <- st_read('bin/British_Isles_map.gpkg')

inset_uk_map <- ggplot() +
    geom_sf(data = inset_uk) +
  geom_rect(aes(xmin = -7.577820,
                ymin = 57.464158,
                xmax = -7.009277,
                ymax = 57.747412), col = 'purple', fill = NA, linewidth = 0.2) +
  theme_void() +
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(colour = 'black', linewidth = 1, fill = NA),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())

ggsave(inset_uk_map, path = 'figs/FIG1/', filename = 'UK.png', width = 1, height = 1.5, units = 'in')

### Load in roost sites for inclusion on the map.
roost_sites <- st_read('bin/roost_sites.gpkg') %>%
  st_centroid()

### Create Figure 1A to which the inset will be added; depict movement of all Rock Doves.
(fig1A <- ggmap(clachan_sat) +
    geom_point(data = gps_wgs_coords, aes(x = X, y = Y, col = name), size = 0.3) + # add points
    geom_path(data = gps_wgs_coords, aes(x = X, y = Y, col = name), size = 0.3) + # add paths
    geom_sf(data = roost_sites, inherit.aes = F, size = 2, shape = 25, col = 'red', fill = 'black') +
    coord_sf(crs = st_crs(4326)) +
    xlab('Longitude (°)') + ylab('Latitude (°)') + # plot labels
    scale_y_continuous(limits = c(57.65,57.68408)) + # cropping scale; satellite is larger than clachan
    scale_color_discrete(name = 'Dove ID') + # here you can add colours you like using values = c('colou1', 'colour2')
    ggsn::scalebar(y.min = 57.652, x.min = -7.298783, y.max = 57.68, x.max = -7.198, # adds the scale bar
                   transform = T, dist = 500, dist_unit = 'm', location = 'bottomright',
                   box.fill = c("yellow", "white"), height = 0.02, st.dist = 0.05, st.color = 'white', st.size = 10/.pt) +
    geom_segment(aes(x = -7.29, y = 57.68, xend = -7.29, yend = 57.683), linewidth = 1, col = 'white',
                 arrow = arrow(length = unit(0.25, "cm"))) + # adds the north arrow
    geom_text(aes(x = -7.29, y = 57.684), label = 'N', col = 'white', size = 14/.pt) + # add the north label
    theme(legend.position = 'top'))

ggsave(fig1A, path = 'figs/FIG1/', filename = 'Fig1A.png', dpi = 300, width = 7, height = 4.5, units = 'in')

### Create Figure 1C; a day in life of EA49568.
ea49568 <- gps_wgs_coords %>% filter(name == 'EA49568')
ea49568 <- ea49568[1:67,] %>% mutate(hour=format(as.POSIXct(datetime), format = "%H")) 
ea49568$hour <- as.numeric(ea49568$hour)

(fig1C <- ggmap(clachan_sat) +
    geom_path(data = ea49568, aes(x = X, y = Y, col = hour), size = 1.5) + # add paths
    geom_sf(data = roost_sites, inherit.aes = F, size = 4, shape = 25, col = 'red', fill = 'black') +
    coord_sf(crs = st_crs(4326)) +
    xlab('Longitude (°)') + ylab('Latitude (°)') + # plot labels
    scale_y_continuous(limits = c(57.668, 57.676)) + # cropping scale; satellite is larger than clachan
    scale_x_continuous(limits = c(-7.25, -7.22)) +
    scale_color_viridis_c(name = 'Hour of day (h)', breaks = c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23), option = 'plasma') +
    ggsn::scalebar(y.min = 57.668, x.min = -7.25, y.max = 57.676, x.max = -7.22, # adds the scale bar
                   transform = T, dist = 200, dist_unit = 'm', location = 'bottomright',
                   box.fill = c("yellow", "white"), st.color = 'white', height = 0.05, st.dist = 0.05, st.size = 0.5, st.bottom = F) +
    theme(legend.position = 'top',
          legend.key.width = unit(1, 'cm')))

ggsave(fig1C, path = 'figs/FIG1/', filename = 'Fig1C.png', width = 4.5, height = 3.5, units = 'in')


### 5 == HABITAT USE BY ROCK DOVES ####################################
### Divide the day into half-hour periods. Across all dove days (each Dove has a day), on how many days was a
### Dove present in that habitat in that time?

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

rm(full_h, full_h1, half_h, half_h1)

### 7 == FIGURE 3 - USE OF HABITATS BY ROCK DOVES ###################################
### Merge GPS data with time intervals. Add a habitat type for roost site.
gps_int_roost <- date_df %>%
  left_join(gps %>% 
              mutate(habitat = ifelse(inroost == 'Y', 'Roost site', gps$habitat)), by = character()) %>% # cross-join add creating roost as a new habitat.
  filter(datetime %within% int) %>% # verify if a timepoint is within an interval.
  filter(date > as.Date('2023-07-10')) # add to obtain no first day dataset.

gps_int_roost <- gps_int_roost %>% # simplify habitat names; shift minor habitats into unknown/other
  mutate(habitat = recode(habitat, "H21A0 - Machairs" = 'Machair', "H2130 - Fixed dunes" = "Dunes",
        "Atlantic Cynosurus-Centaurea pastures " = "Atlantic pasture",
        "H7140 - Transition mires" = "Unknown or other",
        "Unknown" = "Unknown or other",
        "[Carex dioica], [Carex pulicaris] and [Carex flava] fens" = 'Unknown or other',
        "Buildings of cities, towns and villages" = 'Unknown or other',
        "H2120 - Shifting dunes" = "Dunes",
        "Agriculturally-improved, re-seeded and heavily fertilised grassland, including sports fields and grass lawns" = "Arable land",
        "Constructed, industrial and other artificial habitats" = "Unknown or other",
        "Unvegetated sand beaches above the driftline" = "Unknown or other",
        "H2190 - Humid dune slacks" = 'Dunes',
        "Flood swards and related communities" = "Unknown or other"))

# Time spent in habitat; per each Dove day was the Dove present in that habitat or not? What is the proportion of presences
# in that habitat at that time by Rock Doves?
time_in_habitat_roost <- gps_int_roost %>%
  group_by(name, date, id, habitat) %>% # for each pigeon, on each day, across each habitat and interval.
  count() %>% # count how many times each pigeon on a given day used which habitats in which interval.
  ungroup() %>% # remove grouping
  group_by(id, habitat) %>% # now by each mix of interval and habitat
  count() # count how many there are.

# The result is: sum of unique combinations of name/date/interval/habitat (or being present in a habitat on a Dove day
# at that interval) across the entire sampling period.

time_in_habitat_roost$id <- as.numeric(time_in_habitat_roost$id) # change interval id to numeric

hours <- date_df %>% dplyr::select(start, id) %>% # create actual hours object to add it as a column
  mutate(hour = paste(format(as.POSIXct(start), format = "%H:%M:%S"))) %>%
  dplyr::select(hour, id) %>%
  unique()
time_in_habitat_roost <- left_join(time_in_habitat_roost, hours, by = 'id') %>%
  filter(habitat != 'unknown')

# Calculate proportional presence in habitat per each hour of day.
time_in_habitat_roost <- time_in_habitat_roost %>%
  group_by(hour) %>%
  mutate(prop = n/sum(n))

## Figure 2A === Create a plot of proportional occurence in habitat.
(fig2A <- ggplot(data = time_in_habitat_roost, aes(x = as.POSIXct(hour, format = "%H:%M:%S"), y = n, fill = habitat)) + 
    geom_bar(position = 'fill', stat = 'identity') + theme_bw() +
    xlab('BST Time of day (half-hour intervals)') + ylab('Proportion of presences in habitat') +
    scale_x_datetime(date_breaks = "1 hour", date_labels = "%H:00", expand = c(0, 0)) +
    scale_fill_manual(name = '', values = c("Arable land" = "#DDCC77", "Atlantic pasture" = '#CC6677', "Roost site" = 'grey',
                                            "Machair" = '#117733', "Dunes" = '#6699CC', 'Unknown or other' = '#332288'))  +
    theme(legend.position = 'top',
          axis.text.x = element_text(size = 10, angle = 90),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 10),
          axis.title.x = element_text(margin = margin(t = 10)),
          axis.title.y = element_text(margin = margin(t = 10)),
          legend.text = element_text(size = 9),
          legend.title = element_text(size = 10),
          plot.margin = margin(1,0,0,0, 'cm')) +
    guides(fill = guide_legend(nrow = 2)))

ggsave(fig2A, path = 'figs/FIG2/', filename = 'Fig2A.png', dpi = 300, width = 7, height = 4.5, units = 'in')

## Figure 2B === Plot the study area and North Uist with distribution of habitats of interest.

## Load in and recode habitat names in a Uist habitat map:
uist_habitats <- st_read('bin/clachan_new_habitat.gpkg') %>% st_transform(crs = st_crs(4326)) %>%
  st_crop(ymin = 57.69, ymax = 57.64, xmin = -7.29, xmax = -7.19) %>%
  mutate(habitat = ifelse(habitat == "Machair", 'Machair',
                          ifelse(habitat == "H2130 - Fixed dunes", "Dunes", 
                                 ifelse(habitat == "Atlantic Cynosurus-Centaurea pastures ", "Atlantic pasture",
                                        ifelse(habitat == "H2120 - Shifting dunes", "Dunes",
                                               ifelse(habitat == "Arable land", 'Arable land',
                                               ifelse(habitat == "Agriculturally-improved, re-seeded and heavily fertilised grassland, including sports fields and grass lawns", "Arable land",
                                                      ifelse(habitat == "H2190 - Humid dune slacks", 'Dunes',
                                                             ifelse(habitat == "H2150 - Decalcified fixed dunes", 'Dunes',
                                                                    ifelse(habitat == "H2110 - Embryonic shifting dunes", 'Dunes',
                                                                           ifelse(habitat =="Shifting coastal dunes", 'Dunes', 'Unknown or other')))))))))))


(fig2B <- ggmap(clachan_sat) +
  geom_sf(data = uist_habitats, aes(fill = habitat), inherit.aes = F, alpha = 0.5) +
  geom_point(data = gps_wgs_coords, aes(x = X, y = Y), size = 0.5, alpha = 0.5) + # add points
  geom_path(data = gps_wgs_coords, aes(x = X, y = Y), size = 0.3, alpha = 0.5, col = 'grey') + # add paths
  geom_sf(data = roost_sites, inherit.aes = F, size = 2, shape = 25, col = 'red', fill = 'black') +
  coord_sf(crs = st_crs(4326)) +
  xlab('Longitude (°)') + ylab('Latitude (°)') + # plot labels
  scale_y_continuous(limits = c(57.65,57.68408)) +
  scale_fill_manual(name = '', values = c("Arable land" = "#DDCC77", "Atlantic pasture" = '#CC6677',
                                          "Machair" = '#117733', "Dunes" = '#6699CC', 'Unknown or other' = '#332288')) +
  ggsn::scalebar(y.min = 57.652, x.min = -7.298783, y.max = 57.68, x.max = -7.198, # adds the scale bar
                   transform = T, dist = 500, dist_unit = 'm', location = 'bottomright',
                   box.fill = c("yellow", "white"), height = 0.02, st.dist = 0.05, st.color = 'white', st.size = 10/.pt) +
  geom_segment(aes(x = -7.29, y = 57.68, xend = -7.29, yend = 57.683), linewidth = 1, col = 'white',
                 arrow = arrow(length = unit(0.25, "cm"))) + # adds the north arrow
  geom_text(aes(x = -7.29, y = 57.684), label = 'N', col = 'white', size = 14/.pt) + 
  theme(legend.position = 'none'))

ggsave(fig2B, path = 'figs/FIG2/', filename = 'Fig2B.png', dpi = 300, width = 7, height = 4.5, units = 'in')
