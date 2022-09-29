is.element("sf", installed.packages())

install.packages("sf", dep=TRUE)
library(INLA)         # Integrated Nested Laplace Approximation package
library(dplyr)        # A package for data manipulation
library(sf)           # Simple feature for R
library(spdep)        # Functions and tests for evaluating spatial patterns 
# and autocorrelation
library(SpatialEpi)   # Methods and data for spatial epidemiology

# Packages used for visualization 
library(RColorBrewer) # A package providing colour palettes for shading maps 
# and other plots
library(tmap)         # A package for static and interactive maps
library(ggplot2)      # A package that implements the grammar of graphics, which is a term used to
# break up graphs into semantic components, such as geometries and layers.
library(mapview)      # A package for interactive maps
library(cowplot)      # Add-on to ggplot. It provides features that help us with creating
# publication-quality figures

LTLA <- st_read("Practical2a/LTLA_shp.shp")

class(LTLA)

# Check geometry
st_geometry_type(LTLA) 

# Check what CRS this file data is in
st_crs(LTLA) # the data are encoded using a Transverse Mercator Projection. 
# The Airy ellipsoid is being used (+ellps=airy) and the units are meters (+units=m)

# Check the spatial extent of the shapefile (i.e. the geographic "edge" or location that is the furthest north, south east and west) 
st_bbox(LTLA)

# View all of the metadata and attributes for this shapefile object
LTLA

# the default plot of an sf object is a multi-plot of all attributes

plot(LTLA)          # plot all the attributes
plot(LTLA$geometry) # plot only the boundaries

## Static Map with red and white
ggplot() + 
  geom_sf(data = LTLA, color = "red", fill = "white") + 
  ggtitle("Map of LTLAs in England") + 
  coord_sf() +    #axis limits and CRS
  labs(x = "Longitude", y = "Latitude", fill = "") +
  theme_bw() +    # dark-on-light theme
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

## Static map for LTLA
tm_shape(LTLA) +
  tm_fill("tomato") +
  tm_borders(lty="dashed", col="gold") +
  tm_style("natural", bg.color="grey90") +
  tm_layout(title="Map of LTLAs in England")

# as seen before, but with specified presentation mode
tmap_mode("plot")
tm_shape(LTLA) +
  tm_fill("grey90") +
  tm_borders(lty="solid", col="skyblue1") +
  tm_style("albatross") +
  tm_layout(title="Static Map of LTLAs in England")

# interactive mode
tmap_mode("view")
tm_shape(LTLA) +
  tm_fill("grey90") +
  tm_borders(lty="solid", col="black") +
  tm_style("natural") +
  tm_layout(title="Interactive Map of LTLAs in England")

mapview(LTLA)


## Covid deaths expected and SMR

COVID19Deaths<-read.csv("Practical2a/COVID19Deaths.csv")

COVID19Deaths$SMR <- COVID19Deaths$O/COVID19Deaths$E

England_SMR <- left_join(LTLA, COVID19Deaths, by = c("LTLA" = "LTLA"))

plot(England_SMR["O"])

plot(England_SMR["SMR"], breaks = c(0,0.5,1,1.5,2))

ggplot() + 
  geom_sf(data = England_SMR, aes(fill = SMR), col = NA) +
  scale_fill_viridis_c(limits = c(0,2), option = "turbo") +
  # colorspace::scale_fill_continuous_divergingx("RdYlBu", limits = c(0,2))+
  theme_bw() -> MapSMR1

MapSMR1

tmap_mode("plot")
tm_shape(England_SMR) + 
  tm_polygons("SMR",palette="RdYlGn", style="cont", n=8) +
  tm_borders(lty="solid", col="black") +
  tm_style("natural") +
  tm_layout(title="Map of SMRs at LTLA level in England")

tmap_mode("plot")
tm_shape(England_SMR) + 
  tm_polygons("IMD",style="pretty", n=8, alpha=0.5) +
  tm_borders(lty="solid", col="black") +
  tm_style("natural") +
  tm_layout(title="Map of IMD at LTLA level in England")

breaks =  c(0,0.5,1,1.5,2)
England_SMR <- mutate(England_SMR, SMR_cat = cut(SMR, breaks, include.lowest = TRUE)) # mutate() adds new variables 
# and preserves existing ones

ggplot() + 
  geom_sf(data = England_SMR, aes(fill = SMR_cat), col = NA) +
  theme_bw() + 
  scale_fill_brewer(palette = "OrRd") + 
  guides(fill=guide_legend(title="SMR"))

breaks =  c(0,0.5,1,1.5,2)
England_SMR <- mutate(England_SMR, SMR_cat = cut(SMR, breaks, include.lowest = TRUE)) # mutate() adds new variables 
# and preserves existing ones

ggplot() + 
  geom_sf(data = England_SMR, aes(fill = SMR_cat), col = NA) +
  theme_bw() + 
  scale_fill_brewer(palette = "OrRd") + 
  guides(fill=guide_legend(title="SMR"))


ID<- seq(1,317)
formula_iid <- O ~ f(ID, model="iid",
                     hyper=list(prec = list(
                       prior = "pc.prec",
                       param = c(0.5 / 0.31, 0.01))))  
mod_iid <- inla(formula=formula_iid, 
                family="poisson", 
                data=England_SMR, 
                E=E, 
                control.compute=list(dic=TRUE, waic=TRUE))

mod_iid$waic$waic
mod_iid$dic$dic

formula_reg <- O ~ IMD + f(ID, model="iid",
                           hyper=list(prec = list(
                             prior = "pc.prec",
                             param = c(0.5 / 0.31, 0.01))))  
mod_reg <- inla(formula=formula_reg, 
                family="poisson", 
                data=England_SMR, 
                E=E, 
                control.compute=list(dic=TRUE, waic=TRUE))
mod_reg$waic$waic


RR_COVID<-c()
for(i in 1:317){
  RR_COVID[i] <- inla.emarginal(function(x) exp(x), 
                                mod_iid$marginals.random$ID[[i]])
}

England_SMR$RR <- RR_COVID

ggplot() + 
  geom_sf(data = England_SMR, aes(fill = RR), col = NA) + 
  theme_bw() + 
  scale_fill_viridis_c(limits = c(0,2), option = "turbo") -> MapRR1

MapRR1

plot_grid(MapSMR1, MapRR1, ncol = 2, align = 'v', labels="AUTO", rel_widths = c(1, 1))
