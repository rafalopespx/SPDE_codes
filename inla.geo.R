library(spData); library(sf); library(tidyverse)

world <- st_read(system.file("shapes/world.gpkg", package="spData"))
plot(world["pop"])
class(world)

world_sp <- as(world, "Spatial")
class(world_sp)

world_sf <- st_as_sf(world_sp)
class(world_sf)  

world %>% select(lifeExp) %>% plot()

worlg_geo <- st_geometry(world)

world %>% st_geometry() %>% plot()
