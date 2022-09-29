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
