library(tidyverse)
library(spTimer)
library(mapview) 
library(lubridate)
library(mvtsplot)
library(INLA)
library(inlabru)
library(viridis)

data(NYdata)
dim(NYdata)

head(NYdata)

stations <- cbind(unique(NYdata[,1]), unique(NYdata[,2:3]))

# set the map projection to a common projection standard such as WGS84 via the argument crs = 4326
mapview(stations, xcol = "Longitude", ycol = "Latitude", crs = 4326, grid = FALSE)

NYdata <- NYdata %>%
  mutate(date = make_date(Year, Month, Day))
glimpse(NYdata)

# select ozone data
O3 <- NYdata %>% select(s.index, o8hrmax, date)
dim(O3) 

O3_wide = O3 %>% spread(s.index, o8hrmax)
dim(O3_wide) 

O3_wide <- O3_wide[,-1] # remove date
O3_wide <- data.matrix(O3_wide)
dim(O3_wide)

colnames(O3_wide) <- unique(O3[,1])

# Daily ozone levels for 28 monitoring stations, Jul 1– Aug 31, 2006.
mvtsplot(O3_wide, 
         group = NULL, 
         xtime = NULL, 
         norm = c("global"),
         levels = 3, 
         smooth.df = NULL, 
         margin = TRUE, 
         sort =NULL,
         main = "", 
         palette = "PRGn", 
         rowstat = "median", 
         xlim,
         bottom.ylim = NULL, 
         right.xlim=NULL, 
         gcol = 3)

summary(NYdata$o8hrmax)

NYdata %>% 
  ggplot()+
  geom_histogram(aes(o8hrmax), 
                 col="orange")

library(GGally)
ggpairs(NYdata[,7:10]) # print correlations between variables

library(corrr)
tab_cor = NYdata %>%
  select(o8hrmax, cMAXTMP, WDSP, RH) %>%
  correlate() %>%
  shave(upper = TRUE) %>%
  fashion(decimals = 2, na_print = "—") 

tab_cor

NYdata$time = rep(1:n_distinct(NYdata$date),
                  n_distinct(NYdata$s.index))

NYdata = NYdata %>% filter(time <= 30)

NYdata$sqrto8hrmax = sqrt(NYdata$o8hrmax)
NYdata %>% 
  ggplot()+
  geom_histogram(aes(sqrto8hrmax), col="orange")

# validation sites
s <- c(8, 11, 12, 14, 18, 21, 24, 28) # 8 validation stations

# Training data
DataFit = NYdata %>% 
  filter(!(s.index %in% s))
dim(DataFit)

# Dim = 600 (20stations*30days) * 12

# Validation data
DataValPred =  NYdata %>% 
  filter(s.index %in% s)
dim(DataValPred)

bnd = inla.nonconvex.hull(cbind(DataFit$Longitude, DataFit$Latitude),
                          convex = 1)

mesh = inla.mesh.2d(loc = cbind(DataFit$Longitude, DataFit$Latitude),
                    max.n.strict = c(100, 20))

ggplot() +
  gg(mesh) +
  geom_point(data = DataFit, aes(Longitude, Latitude)) 

spde = inla.spde2.matern(mesh = mesh)
spde$n.spde #n. of mesh vertices

coordinates(DataFit) = c("Longitude","Latitude")
class(DataFit)

## A model with intercept, 3 linear effects of MaxTemp, 
# WindSpeed and Relative Humidty and spatial temporal field by an AR(1)
cmp  = sqrto8hrmax ~ Intercept(1) + cMAXTMP + WDSP + RH +
  SPDE(coordinates, model = spde,
       group = time, control.group = list(model = "ar1")) 

library(inlabru)
lik = like(formula = sqrto8hrmax ~ Intercept  + cMAXTMP + WDSP + RH + SPDE,
           family = "gaussian",
           data = DataFit)

#prob(sigma > sigma0)=alpha
#prob(sigma > 0.2) = 0.1
pc.prec = list(prec = list(prior = "pc.prec", param = c(0.2, 0.1)))

fit = bru(cmp, lik,
          options =  list(control.family = list(hyper = pc.prec)))
summary(fit)

names(fit$marginals.fixed)
names(fit$marginals.random)
names(fit$marginals.hyperpar)
fit$summary.fixed[,c("mean","0.025quant","0.975quant")]

int.plot <- plot(fit, "Intercept")
cMAXTMP.plot <- plot(fit, "cMAXTMP")
WDSP.plot <- plot(fit, "WDSP")
RH.plot <- plot(fit, "RH")

multiplot(int.plot, cMAXTMP.plot, WDSP.plot, RH.plot, ncols = 2)

# SPDE is the chosen name for the spatial field
spde.range <- spde.posterior(fit, "SPDE", what = "range")
spde.var <- spde.posterior(fit, "SPDE", what = "variance")
range.plot <- plot(spde.range)
var.plot <- plot(spde.var)

multiplot(range.plot, var.plot)

coordinates(DataValPred) = c("Longitude","Latitude")
class(DataValPred)

ValPred = predict(fit, 
                  DataValPred,
                  ~ (Intercept + SPDE + cMAXTMP + WDSP + RH)^2,## returning to the original scale
                  nsamples = 200)
head(ValPred@data)

as.data.frame(ValPred) %>%
  ggplot() + 
  geom_line(aes(time, o8hrmax, group = s.index), color = "red") +
  geom_line(aes(time, median, group = s.index)) +
  geom_ribbon(aes(time, ymin = q0.025, ymax = q0.975, group = s.index), alpha = 0.5) +
  facet_wrap(~ s.index)

data(NYgrid)
dim(NYgrid)

class(NYgrid)

head(NYgrid)

NYgrid %>% 
  distinct(Longitude, Latitude) %>% 
  ggplot()+
  geom_point(aes(Longitude, Latitude)) +
  gg(mesh)

NYgrid <- NYgrid %>%
  mutate(date = make_date(Year, Month, Day))


NYgrid$time = rep(1:n_distinct(NYgrid$date),
                  n_distinct(NYgrid$s.index))

glimpse(NYgrid)

NYgrid = NYgrid %>% filter(time <= 6)

coordinates(NYgrid) = c("Longitude","Latitude")
class(NYgrid)

grid = sp::SpatialPixelsDataFrame(NYgrid@coords,
                                  data = NYgrid@data)
class(grid)

dim(grid@data)

glimpse(grid@data)

GridPred = predict(fit, grid,
                   ~ (Intercept + SPDE + cMAXTMP + WDSP + RH)^2)
head(GridPred@data)

ggplot() + 
  gg(GridPred, aes(Longitude, Latitude, fill = median))  +
  facet_wrap(~ time) + 
  scale_fill_viridis() +
  coord_equal()+
  geom_point(data = stations, aes(Longitude, Latitude))
