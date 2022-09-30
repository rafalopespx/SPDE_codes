df = readRDS("./data/Piemonte_Data.rds")
class(df)

head(df)

df = df[df$time <= 50,]

library(tidyverse)
library(inlabru)

border = readRDS("./data/Piemonte_Border.rds")
class(border)

ggplot()+
  gg(border) +
  coord_equal()

library(tidyverse)
library(viridis)

df %>% 
  filter(time<=3) %>%
  ggplot() + 
  geom_point(aes(UTMX, UTMY, color = PM10), size = 2)+
  facet_wrap(.~ Date, ncol = 2, nrow = 2) +
  scale_color_viridis() +
  coord_equal() +
  gg(border)

library(INLA)

mesh = inla.mesh.2d(loc = cbind(df$UTMX, df$UTMY),
                    offset = c(20, 40),
                    max.edge = c(30, 50))

ggplot() +
  gg(mesh) +
  geom_point(data = df, aes(UTMX, UTMY)) + 
  gg(border)

spde = inla.spde2.matern(mesh = mesh)
spde$n.spde #n. of mesh vertices

coordinates(df) = c("UTMX","UTMY")
class(df)

cmp  = logPM10 ~ Intercept(1) + 
  SPDE(coordinates, 
       model = spde,
       group = time, 
       control.group = list(model = "ar1")) +
  A + #dem(A, model = "linear") +  
  TEMP #temp(TEMP, model = "linear")

lik = like(formula = logPM10 ~ Intercept + SPDE + A + TEMP,
           family = "gaussian",
           data = df)

fit = bru(cmp, lik)
fit$summary.fixed[,c("mean","0.025quant","0.975quant")]

pred_at_station = predict(fit, 
                          df, 
                          ~ Intercept + SPDE + A + TEMP, 
                          n.samples = 1000)

sel = c(1, 7, 10, 19, 23, 24)

as.data.frame(pred_at_station) %>%
  dplyr::filter(Station.ID %in% sel) %>% 
  ggplot() + 
  geom_line(aes(time, logPM10, group = Station.ID), color = "red") +
  geom_line(aes(time, mean, group = Station.ID)) +
  geom_ribbon(aes(time, ymin = q0.025, ymax = q0.975, group = Station.ID), alpha = 0.5) +
  facet_wrap(.~Station.ID)

covariate_grid = readRDS("./data/covariate_grid.rds")
class(covariate_grid)

head(covariate_grid@data)

ggplot() + 
  gg(covariate_grid, aes(fill=A)) + 
  gg(border) +
  coord_equal() + 
  scale_fill_viridis()

ggplot() + 
  gg(covariate_grid, aes(fill=TEMP)) +
  facet_wrap(~ time) +
  gg(border) +
  coord_equal() +
  scale_fill_viridis()

pred = predict(fit, 
               covariate_grid, 
               ~ Intercept + SPDE + A + TEMP,
               seed = 2, 
               n.samples = 1000)
head(pred@data)

predprob = predict(fit, 
                   covariate_grid, 
                   ~ (Intercept + SPDE + A + TEMP) > log(50),
                   seed = 2, 
                   n.samples = 500)

ggplot() + 
     gg(predprob, aes(UTMX_km, UTMY_km, fill = mean))  +
     facet_wrap(.~ time, ncol = 2, nrow = 2) + 
     scale_fill_viridis() + coord_equal() + gg(border)
