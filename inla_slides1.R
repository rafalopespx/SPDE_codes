a <- 1
b <- 1
theta <- rbeta(1,a,b)
n <- 1000
y <- rbinom(1, size=n, p=theta)
a1 <- a + y
b1 <- n - y + b
sim <- rbeta(n=50000, shape1=a1, shape2=b1) 
logodds <- log(sim/(1-sim))
plot(logodds)
hist(logodds)

summary(iris)
formula <- Petal.Length ~ 1 + Petal.Width
output <- inla(formula, 
               family = "gaussian", 
               data = iris)
output$summary.fixed
output$summary.hyperpar
names(output$marginals.fixed)
beta1_post <-output$marginals.fixed[[2]]
marg <- inla.smarginal(beta1_post)
q <-inla.qmarginal(0.05,beta1_post)

plot(marg,t="l",
     ylab="",xlab="", 
     main=expression(paste("p(",beta[1], "| y)")))
polygon(c(marg$x[marg$x <= q ], q),
          c(marg$y[marg$x <= q ], 0),
          col = "slateblue1", border = 1)

inla.pmarginal(q,beta1_post)
d <-inla.dmarginal(q,beta1_post)

prec_post <-output$marginals.hyperpar[[1]]

