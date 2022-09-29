## Lip cancer in scotland

library(SpatialEpi)

LipCancer <- read.csv("scotlip.csv")
LipCancer

formula.inla <- y ~ 1 + 
  f(RECORD_ID,model="iid", 
    hyper=list(prec=list(prior="loggamma",
                         param=c(1,0.01))))

lipcancer.poisson <- inla(formula.inla,family="poisson",
                          data=LipCancer, E=E,
                          control.predictor=list(compute=TRUE),
                          control.compute=list(config=TRUE),
                          control.fixed=list(mean.intercept=0,prec.intercept=0.00001))

sigma.v<- inla.tmarginal(function(x) sqrt(1/x),
                         lipcancer.poisson$marginals.hyperpar[[1]])

inla.qmarginal(seq(0,1,0.2),sigma.v)

joint.post <- inla.posterior.sample(100,lipcancer.poisson)
names(joint.post[[1]])

joint.post[[1]]$latent[1:3,]