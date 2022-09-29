library(tidyverse)
library(INLA)
library(inlabru)
library(viridis) #for colors

# Remove old objects
remove(list=ls())
# Remember to set the correct working directory
load("temperature.croatia.Rdata")
ls()
glimpse(stations_data)
glimpse(grid_cov)

ggplot()+
  geom_point(data = grid_cov, aes(Lon,Lat), col = "gray") +
  geom_point(data = stations_data, aes(Lon,Lat), col = "red", size = 2) +
  theme_bw()

stations_data %>% 
  ggplot() +
  geom_histogram(aes(MDTEMP, after_stat(density)),bins = 13, alpha=0.5) +
  geom_density(aes(MDTEMP)) 

stations_data %>% 
  ggplot() +
  geom_point(aes(Lon, Lat, col = HRdem), size = 2) +
  scale_color_viridis()

grid_cov %>% 
  ggplot() +
  geom_raster(aes(Lon, Lat, fill=HRdem)) +
  scale_fill_viridis()

bnd = inla.nonconvex.hull(cbind(stations_data$Lon,stations_data$Lat), 
                          convex=0.25)

croatia.mesh = inla.mesh.2d(loc = cbind(stations_data$Lon,
                                        stations_data$Lat),
                            boundary = bnd,
                            offset = c(1, 2),
                            max.edge = c(3, 8),
                            cutoff = 0.3)
plot(croatia.mesh,asp=1)
points(stations_data$Lon, stations_data$Lat, pch=21, cex=1.2, bg="white",col=1)

ggplot() + 
  gg(croatia.mesh) +
  geom_point(data = stations_data, aes(Lon, Lat))

spde = inla.spde2.matern(mesh = croatia.mesh)
spde$n.spde

A.est = inla.spde.make.A(croatia.mesh,
                         loc = as.matrix(cbind(stations_data$Lon,stations_data$Lat)))
dim(A.est)

stack.est = inla.stack(data = list(temp = stations_data$MDTEMP),
                       A = list(A.est, 1, 1),
                       effects = list(spatial.index = 1:spde$n.spde,
                                      Intercept = rep(1, nrow(stations_data)),
                                      HRdem = stations_data$HRdem), 
                       tag="est")

A.pred = inla.spde.make.A(croatia.mesh,
                          loc = cbind(grid_cov$Lon,grid_cov$Lat))
dim(A.pred)

stack.pred = inla.stack(data = list(temp = NA), 
                        A = list(A.pred, 1, 1),
                        effects = list(spatial.index = 1:spde$n.spde,
                                       Intercept = rep(1, nrow(grid_cov)),
                                       HRdem = grid_cov$HRdem),
                        tag = "pred")

fullstack = inla.stack(stack.est, stack.pred)

formula = temp ~ -1 + Intercept + HRdem + f(spatial.index, model = spde)

output = inla(formula,
              data = inla.stack.data(fullstack, spde = spde),
              family = "gaussian",
              control.predictor = list(A = inla.stack.A(fullstack), compute = TRUE),
              control.compute = list(dic = TRUE))
output$summary.fixed

inla.smarginal(output$marginals.fixed$Intercept) %>% 
  dplyr::bind_rows() %>%  #from list to data frame
  ggplot() +
  geom_line(aes(x,y)) +
  ggtitle("Intercept")


inla.smarginal(output$marginals.fixed$HRdem) %>% 
  dplyr::bind_rows() %>%  #from list to data frame
  ggplot() +
  geom_line(aes(x,y)) +
  ggtitle("Elevation")

sigma2e_marg =
  inla.tmarginal(function(x) 1/x,
                 output$marginals.hyperpar$"Precision for the Gaussian observations")

inla.zmarginal(sigma2e_marg)

output.field = inla.spde2.result(inla = output, 
                                 name = "spatial.index",
                                 spde = spde,
                                 do.transf = TRUE)
var.nom.marg = output.field$marginals.variance.nominal[[1]]
inla.zmarginal(var.nom.marg)

range.nom.marg = output.field$marginals.range.nominal[[1]]
inla.zmarginal(range.nom.marg)

output$dic$dic

coordinates(stations_data) = c("Lon","Lat")
class(stations_data)

cmp1 = MDTEMP ~ Intercept(1) + HRdem + s.field(main = coordinates, model = spde) 

lik1 = like(formula = MDTEMP ~ Intercept + HRdem + s.field,
            family = "gaussian",
            data = stations_data)
outputbru <- bru(cmp1, lik1)

cbind(
  rbind(output$summary.fixed[,c("mean", "0.5quant", "sd")],
        outputbru$summary.fixed[,c("mean", "0.5quant","sd")]),
  method = rep(c("inla.stack", "inlabru"), each=2)
)

cbind(
  rbind(output$summary.hyperpar[,c("mean", "0.5quant", "sd")],
        outputbru$summary.hyperpar[,c("mean", "0.5quant","sd")]),
  method = rep(c("inla.stack", "inlabru"), each=3)
)

index.pred = inla.stack.index(stack = fullstack, "pred")$data

grid_cov$post.mean.pred = output$summary.linear.predictor[index.pred, "mean"]
grid_cov$post.sd.pred = output$summary.linear.predictor[index.pred, "sd"]

grid_cov %>% 
  ggplot() +
  geom_raster(aes(Lon, Lat, fill = post.mean.pred)) +
  scale_fill_viridis() +
  ggtitle("Posterior mean - inla.stack") 

grid_cov %>% 
  ggplot() +
  geom_raster(aes(Lon, Lat, fill = post.sd.pred)) +
  scale_fill_viridis() +
  ggtitle("Posterior SD - inla.stack") 

coordinates(grid_cov) = c("Lon", "Lat")
gridded(grid_cov) = TRUE
class(grid_cov)

pred = predict(outputbru, grid_cov,
               ~ Intercept + HRdem + s.field,
               seed = 1,
               n.samples = 500)
head(pred@data)

ggplot() + 
  gg(pred, aes(Lon, Lat, fill = mean)) +
  ggtitle("Posterior mean - inlabru") + 
  scale_fill_viridis()

cmp2 = MDTEMP ~ 1 + HRdem + HRdsea + s.field(main = coordinates, model = spde) 
lik2 = like(formula = MDTEMP ~ Intercept + HRdem + HRdsea + s.field,
            family = "gaussian",
            data = stations_data)
outputbru2 = bru(cmp2, lik2)

outputbru$dic$dic
outputbru2$dic$dic

output$dic$dic
