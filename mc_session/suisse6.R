#library(tidyverse) # wrangling tabular data and plotting

library(sf) # processing spatial vector data
library(sp) # another vector data package necessary for continuity
# And a lot of different packages to test their interpolation functions
library(gstat)  # inverse distance weighted, Kriging
library(fields) # Thin Plate Spline
library(interp) # Triangulation
library(mgcv)   # Spatial GAM
library(automap)# Automatic approach to Kriging
# Finally, some packages to make pretty plots
#library(patchwork)
library(viridis)
library(spatstat)
library(sf)
library(sfheaders)
library(ggplot2)
#library(dismo)
library(terra)
library(raster)
library(tmap)
library(tmaptools)
library(rspat)
library(geostatsp)
library(ggthemes)
library(dplyr)
# get the data from geospatstat

area=sf::st_read("mc_session/plot.shp")
stations=sf::st_read("mc_session/stations_prelim.gpkg")
statdat = openxlsx::read.xlsx("~/Downloads/all_GW1000A-WIFIFC29(202308290000-202308291731).xlsx")
area_sp=as(area,"Spatial")
plotborder= sp::spTransform(area_sp,CRSobj = "EPSG:32633")

statdat$
# get the extent of switzerland
extent.harzplot <- raster::extent(plotborder)




# get the extent of switzerland
extent.plot<- raster::extent(extent.harzplot)

# gcreate a point data grid  1000m*1000m +border distance
plot_grid <- expand.grid(x = seq(from = round(extent.plot@xmin),
to = round(extent.plot@xmax),
by = 700),
y = seq(from = round(extent.plot@ymin),
to = round(extent.plot@ymax),
by = 700))

# convert data to terra vvector format
plotborder_vect <- terra::vect(plotborder)

# calculate voronoi
t_swiss_rain_voronoi <- terra::voronoi(plotborder_vect )
# reconvert it
sf_swiss_rain_voronoi=st_as_sf(t_swiss_rain_voronoi)

# map it
tmap_mode("view")
tmap::qtm(swissRain) + 
tmap::qtm(sf_swiss_rain_voronoi,fill = NULL) +
tmap::qtm(swissBorder, fill = NULL,borders =  "red")

## plot it
ggplot() +
  geom_point(data = rain_df,
             mapping = aes(x = x, y = y, color = rain), size = 0.8) +
  scale_fill_viridis_d()+
  coord_sf(crs = st_crs(crs))+   theme_bw()
  
  
####

rain_sf = sf::st_as_sf(swissRain)
rain_df = sfheaders::sf_to_df(rain_sf,fill = TRUE)
rain_df = dplyr::select(rain_df,rain,x,y)

bbox = extent(swissRain)
bbox@ymin = bbox@ymin -  30000
bbox@xmin = bbox@xmin -  30000
bbox@ymax = bbox@ymax +  30000
bbox@xmax = bbox@xmax +  30000
# # no 1
# grd_template <- expand.grid(
# X = seq(from = bbox@xmin, to = bbox@xmax, by = 1000),
# Y = seq(from = bbox@ymin, to = bbox@ymax, by = 1000) # 20 m resolution
# )

# spatial points
alt_grd_template_sf <- bbox %>% 
  st_bbox(bbox) %>% 
  st_as_sfc() %>% 
  st_make_grid(
    cellsize = c(1000, 1000),
    what = "centers"
  ) %>%
  st_as_sf() %>%
  cbind(., st_coordinates(.)) %>% 
  st_drop_geometry() %>% 
  mutate(Z = 0)

# raster
alt_grd_template_raster <- alt_grd_template_sf %>% 
  raster::rasterFromXYZ(
    crs = crs)

# Nearest Neighbor
fit_NN <- gstat::gstat( # using package {gstat} 
  formula = rain ~ 1,    # The column `NH4` is what we are interested in
  data = as(rain, "Spatial"), # using {sf} and converting to {sp}, which is expected
  nmax = 10, nmin = 3 # Number of neighboring observations used for the fit
)

# Inverse Distance Weighting
fit_IDW <- gstat::gstat( # The setup here is quite similar to NN
  formula = rain ~ 1,
  data = as(rain_sf, "Spatial"),
  nmax = 10, nmin = 3,
  set = list(idp = 0.5) # inverse distance power
)

# Thin Plate Spline Regression
mtx = cbind(swissRain$x, swissRain$y ,swissRain$rain)
colnames(mtx) <- c("x","y","rain")
fit_TPS <- fields::Tps( # using {fields}
  x = as.matrix(mtx[, c("x", "y")]), # accepts points but expects them as matrix
  Y = rain_df$rain,  # the dependent variable
  miles = FALSE     # EPSG 25833 is based in meters
)

# Generalized Additive Model
fit_GAM <- mgcv::gam( # using {mgcv}
  rain ~ s(x, y),      # here come our X/Y/Z data - straightforward enough
  data = rain_df      # specify in which object the data is stored
)

# Next we use a couple of functions that have a slightly different modus
# operandi as they in fact already return interpolated Z values.

# Triangular Irregular Surface
fit_TIN <- interp::interp( # using {interp}
  x = rain_df$x,           # the function actually accepts coordinate vectors
  y = rain_df$y,
  z = rain_df$rain,
  xo = s_rain_grid$x,     # here we already define the target grid
  yo = s_rain_grid$y,
  output = "points"
) %>% bind_cols()

# Automatized Kriging  
fit_KRIG <- automap::autoKrige( # using {automap}
  formula = rain ~ 1,           # The interface is similar to {gstat} but
  input_data = swissRain        # {automap} makes a lot of assumptions for you
) %>% 
  .$krige_output %>%            # the function returns a complex object with lot's of metainfo
  as.data.frame() %>%           # we keep only the data we are interested in
  dplyr::select(X = x1, Y = x2, Z = var1.pred) 

######

# Triangular Irregular Surface
interp_TIN <- mask(raster::rasterFromXYZ(fit_TIN, crs = crs),swissBorder)
# Automatized Kriging
interp_KRIG <- mask(raster::rasterFromXYZ(fit_KRIG, crs = crs),swissBorder)

# Nearest Neighbor
interp_NN <-mask(raster:: interpolate(alt_grd_template_raster, fit_NN),swissBorder)
# Inverse Distance Weighting
interp_IDW <- mask(raster::interpolate(alt_grd_template_raster, fit_IDW),swissBorder)
# Thin Plate Spline Regression
interp_TPS <- mask(raster::interpolate(alt_grd_template_raster, fit_TPS),swissBorder)

# Generalized Additive Model
names(alt_grd_template_sf) = c("x","y","z")
interp_GAM <- alt_grd_template_sf %>%
mutate(z = predict(fit_GAM, .)) %>%
raster::rasterFromXYZ(crs = crs)
interp_GAM = mask(interp_GAM,swissBorder)

vector=st_intersection(vector,st_as_sf(swissBorder,"Spatial"))


names(interp_NN) <- "Nearest Neighbor"
names(interp_IDW) <- "Inverse Distance Weighted"
names(interp_KRIG) <- "Kriging"
names(interp_TPS) <- "Thin Plate Spline Regression"
names(interp_TIN) <-"Triangular Irregular Surface"
names(interp_GAM) <-"Generalized Additive Model"


plot_my_rasters <- function(raster_object){
  tmap_mode("view")
  tm_shape(raster_object) +
    tm_raster(breaks = c(0,10,20,30,40,50,60),style = "fixed",title = "precp. mm", palette=get_brewer_pal("Blues", n =6, plot=TRUE)) +
    #tm_add_legend(type = "fill")+
    tm_shape(vector) +
    tm_polygons(alpha = 0.0)+
    tm_legend(outside = TRUE)+
    tm_shape(swissRain) +
    tm_symbols(col = "blue",size = "rain",scale =0.2) +
    tm_shape(st_as_sf(swissBorder)) +
    tm_polygons(alpha = 0.0)+
    tm_graticules()+
    tm_layout(title = names(raster_object), main.title = names(raster_object), main.title.size = 2, main.title.position = "top")
  
}

raster_object =  list(
"Nearest Neighbor" = interp_NN,
"Inverse Distance Weighted" = interp_IDW,
"Kriging" = interp_KRIG,
"Thin Plate Spline Regression" = interp_TPS,
"Triangular Irregular Surface" = interp_TIN,
"Generalized Additive Model" = interp_GAM
)

plotlist <-  lapply(raster_object, plot_my_rasters)


m2=tmap_arrange(plotlist,ncol = 2,nrow = 3,sync = TRUE)

