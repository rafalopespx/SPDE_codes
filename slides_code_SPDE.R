library(INLA)
data(SPDEtoy)
dim(SPDEtoy)
head(SPDEtoy, n=3)
summary(SPDEtoy$y)

library(tidyverse)

SPDEtoy %>% 
  ggplot() +
  geom_point(aes(s1,s2, col=y),size=2)

coords <- as.matrix(SPDEtoy[,1:2])
mesh0 <- inla.mesh.2d(loc = coords,max.edge = 0.1)
plot(mesh0)
points(coords)

library(inlabru)

ggplot() + 
  gg(mesh0) +
  geom_point(data = data.frame(coords),
             aes(s1, s2))

mesh1 <- inla.mesh.2d(loc = coords,
                      max.edge = c(0.1, 0.1))

ggplot() + 
  gg(mesh1) +
  geom_point(data = data.frame(coords),
             aes(s1, s2))
mesh2 <- inla.mesh.2d(loc = coords,max.edge = c(0.1, 0.2))

ggplot() +
  gg(mesh2) +
  geom_point(data = data.frame(coords),aes(s1, s2))

mesh3 <- inla.mesh.2d(loc = coords,max.edge = c(0.1, 0.2),offset = c(0.4,0.1))

ggplot()+
  gg(mesh3)+
  geom_point(data = data.frame(coords),aes(s1, s2))

mesh4 <- inla.mesh.2d(loc = coords,
                      max.edge = c(0.1, 0.2),
                      offset = c(0.1,0.4))
ggplot() + 
  gg(mesh4) +
  geom_point(data = data.frame(coords),aes(s1, s2))


domain <- matrix(cbind(c(0,1,1,0.7,0),c(0,0,0.7,1,1)), ncol=2)
mesh5domain <- inla.mesh.2d(loc.domain = domain,
                            max.edge = c(0.04, 0.2), 
                            offset = c(0.1, 0.4))

ggplot() + 
  gg(mesh5domain) +
  geom_point(data = data.frame(domain), aes(X1, X2)) +
  geom_point(data = data.frame(coords),
             aes(s1, s2), col = "red", alpha = 0.5)

mesh5 <- inla.mesh.2d(loc.domain = domain, 
                      max.edge = c(0.04, 0.2), 
                      cutoff = 0.5,
                      offset = c(0.1, 0.4))

ggplot() + 
  gg(mesh5) +
  geom_point(data = data.frame(domain),aes(X1, X2)) +
  geom_point(data = data.frame(coords),
             aes(s1, s2), col = "red", alpha = 0.5)

mesh6 <- inla.mesh.2d(loc.domain = domain,
                      max.edge = c(0.04, 0.2),
                      cutoff = 0.05,
                      offset = c(0.1, 0.4))
ggplot() + 
  gg(mesh6) +
  geom_point(data = data.frame(domain),
             aes(X1, X2)) +
  geom_point(data = data.frame(coords),
             aes(s1, s2), col = "red", alpha = 0.5)


set.seed(44)
loc = matrix(runif(20), 10, 2)
boundary = inla.nonconvex.hull(loc,convex=0.2)
meshNC <- inla.mesh.2d(loc = loc,boundary = boundary,max.edge = c(0.04, 0.2))

ggplot() + 
  gg(meshNC)+
  geom_point(data = data.frame(loc),
             aes(X1, X2))

A.est6 <- inla.spde.make.A(mesh = mesh6,
                           loc = coords)
dim(A.est6)
table(rowSums(A.est6))
table(rowSums(A.est6>0))
table(colSums(A.est6) > 0)

spde = inla.spde2.matern(mesh = mesh6)
spde$n.spde

formula = y ~ -1 + intercept + f(spatial.field, model = spde)

output6 <- inla(formula,
                data = list(y = SPDEtoy$y,
                            intercept = rep(1,spde$n.spde),
                            spatial.field = 1:spde$n.spde), 
                control.predictor = list(A = A.est6, compute = TRUE))
output6$summary.fixed[,c("mean","0.025quant","0.975quant")]
output6$summary.hyperpar[,c("mean","0.025quant","0.975quant")]

output6.field <- inla.spde2.result(inla = output6, 
                                   name = "spatial.field",
                                   spde = spde)
names(output6.field)

inla.emarginal(function(x) x, output6.field$marginals.variance.nominal[[1]])
inla.emarginal(function(x) x, output6.field$marginals.range.nominal[[1]])
inla.zmarginal(output6.field$marginals.range.nominal[[1]])

spde = inla.spde2.pcmatern(mesh6,
                           prior.range = c(0.01,0.1),
                           prior.sigma = c(100,0.1))

library(INLA)
data(SPDEtoy)
dim(SPDEtoy)
head(SPDEtoy, n=3)
summary(SPDEtoy$y)

stack.est <- inla.stack(
  data = list(y = SPDEtoy$y),
  A = list(1, A.est6),
  effects = list(intercept = rep(1,nrow(SPDEtoy)),
                 spatial.field = 1:spde$n.spde),
  tag="est")

formula = y ~ -1 + intercept + f(spatial.field, model = spde)

output6 <- inla(formula,
                data = inla.stack.data(stack.est),
                control.predictor = list(A = inla.stack.A(stack.est), 
                                         compute = TRUE))

output6$summary.fixed[,c("mean","0.025quant","0.975quant")]

output6$summary.hyperpar[,c("mean","0.025quant","0.975quant")]

grid.x = 20
grid.y = 20
pred.grid <- expand.grid(x = seq(0, 1, length.out = grid.x),
                         y = seq(0, 1, length.out = grid.y))
A.pred6 <- inla.spde.make.A(mesh = mesh6,
                            loc = as.matrix(pred.grid))

stack.pred <- inla.stack(data = list(y = NA),
                         A = list(1, A.pred6),
                         effects = list(intercept = rep(1, nrow(pred.grid)),
                                        spatial.field = 1:spde$n.spde),
                         tag = "pred")

join.stack <- inla.stack(stack.est, stack.pred) #full stack object

output6pred <- inla(formula,
                    data = inla.stack.data(join.stack),
                    control.compute = list(return.marginals=TRUE,
                                           return.marginals.predictor=TRUE),
                    control.predictor = list(A = inla.stack.A(join.stack), compute = TRUE))

index.pred <- inla.stack.index(join.stack, tag = "pred")$data
length(index.pred)

output6pred$summary.linear.predictor[index.pred[1:3],c("mean","sd")]

distr.point1 = output6pred$marginals.linear.predictor[index.pred[1]][[1]]
distr.point1.smooth = inla.smarginal(distr.point1)

ggplot(data.frame(distr.point1.smooth))+
  geom_line(aes(x,y))

1 - inla.pmarginal(13, distr.point1)

library(inlabru)
post.mean.pred = output6pred$summary.linear.predictor[index.pred,"mean"]
post.mean.df = SpatialPixelsDataFrame(SpatialPoints(pred.grid),
                                      data = data.frame(post.mean.pred))

library(viridis)
ggplot() +
  gg(post.mean.df, aes(x, y, post.mean.pred)) +
  ggtitle("Posterior mean") + 
  coord_fixed() + #x and y axis with the same length in the plot
  scale_fill_viridis()

post.sd.pred = output6pred$summary.linear.predictor[index.pred,"sd"]
post.sd.df = SpatialPixelsDataFrame(SpatialPoints(pred.grid),
                                    data = data.frame(post.sd.pred))
ggplot()+
  gg(post.sd.df, aes(x, y, post.sd.pred)) +
  ggtitle("Posterior standard deviation") + 
  coord_fixed() +
  scale_fill_viridis()

# Data simulation
n1 <- 200
x1 <- runif(n1)
y1 <- rnorm(n1, mean = 3 + 2 * x1 )
df1 <- data.frame(y = y1, x = x1)

library(inlabru)
# Model component
cmp1 = y ~  Intercept(1) + x

# Likelihood
lik1 = like(formula = y ~ x + Intercept,
            #formula = y ~ .,
            family = "gaussian",
            data = df1)

# Model fit
fit1 <- bru(cmp1, lik1)
fit1$summary.fixed[,c("mean","sd")]


# Model component
cmp2 = y ~  Intercept(1) + beta(x, model="linear")
# Likelihood
lik2 = like(formula = y ~ beta + Intercept,
            family = "gaussian",
            data = df1)
# Model fit
fit2 <- bru(cmp2, lik2)
fit2$summary.fixed[,c("mean","sd")]

coordinates(SPDEtoy) = c("s1","s2")
class(SPDEtoy)

domain <- matrix(cbind(c(0,1,1,0.7,0), 
                       c(0,0,0.7,1,1)),
                 ncol=2)

mesh6 <- inla.mesh.2d(loc.domain = domain,
                      max.edge = c(0.04, 0.2),
                      cutoff = 0.05,
                      offset = c(0.1, 0.4))
spde = inla.spde2.matern(mesh = mesh6)

cmp_toy <- y ~ Intercept(1) + s.field(main = coordinates, model = spde) 

like_toy <- like(formula = y ~ Intercept + s.field,
                 data = SPDEtoy,
                 family = "gaussian")
# Run inlabru!
fit <- bru(cmp_toy, like_toy)
class(fit)

fit$summary.fixed[,c("mean","sd")]

coordinates(pred.grid) = c("x","y")
gridded(pred.grid) = TRUE
class(pred.grid)

pred <- predict(fit, 
                pred.grid, 
                ~ Intercept + s.field, 
                seed = 1, 
                n.samples = 500)
class(pred)

head(pred@data)

ggplot() + 
  gg(pred, aes(x, y, fill = mean)) +
  ggtitle("Posterior mean") + 
  coord_fixed() +
  scale_fill_viridis()

ggplot() + 
  gg(pred, aes(x, y, fill = sd)) +
  ggtitle("Posterior Standard deviation") + 
  coord_fixed() +
  scale_fill_viridis()


