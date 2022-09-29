library(tidyverse)
library(sf)           
library(SpatialEpi)   
library(RColorBrewer)

London <- st_read("Practical2b/GreaterLondon_ward_river.shp")

st_crs(London) 
# You can see that:
# the data are encoded using a Transverse Mercator Projection. 
# the airy ellipsoid is used  
# the units are meters 
London

# the default plot of an sf object is a multi-plot of all attributes

plot(London)          # plot all the attributes
plot(London[,c(1,2)]) # plot selected attributes

ggplot() + 
  geom_sf(data = London, color = "red", fill = "white") + 
  ggtitle("Map of London") + 
  coord_sf() +    # axis limits and CRS
  labs(x = "Longitude", y = "Latitude", fill = "") +
  theme_bw() +    # dark-on-light theme
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

HaemMal <- read_csv("Practical2b/HaemMal.csv") # the structure is a tibble

HaemByWard <- HaemMal %>% 
  group_by(STwardcode, POLY_ID) %>% 
  summarise(cases = sum(CASES), 
            population = sum(POPULATION)) %>% 
  arrange(POLY_ID)

HaemMal <- arrange(HaemMal, POLY_ID, SEX, AGE_GROUP)

HaemByWard$E <- expected(population=HaemMal$POPULATION, 
                         cases=HaemMal$CASES, 
                         n.strata=44)
HaemByWard$SMR <- HaemByWard$cases/HaemByWard$E

London_SMR <- left_join(London, HaemByWard)

# display all the palettes
display.brewer.all()

# display colorblind-friendly brewer palettes
display.brewer.all(colorblindFriendly = TRUE)

# display single RColorBrewer palette by specifying its name 
display.brewer.pal(n = 5, name = "OrRd")

# display hexadecimal color
brewer.pal(n = 5, name = "OrRd")

breaks =  c(0,0.5,1,1.5,2,2.5)
London_SMR <- mutate(London_SMR, SMR_cat = cut(SMR, breaks)) 
# mutate() adds new variable 
# and preserves existing ones

# pick a palette you like and plot the SMRs
MapSMR <- ggplot() + 
  geom_sf(data = London_SMR, aes(fill = SMR_cat)) +
  theme_bw() +     # plot white background and black gridlines
  scale_fill_brewer(palette = "OrRd") + 
  guides(fill=guide_legend(title="SMR"))

MapSMR
