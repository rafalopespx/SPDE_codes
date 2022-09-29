set.seed(1234) #set the seed
n = 100
sigma = 0.1
beta.0 =  2
beta.1 = 0.5
x = runif(n)
eta = beta.0 +  beta.1*x
y = rnorm(n,eta,sigma)
plot(x,y) #plot the data
my.data = data.frame(y=y,x=x)

library(INLA)
formula = y ~ 1 + x

output = inla(formula, family="gaussian", data=my.data)
summary(output)
names(output$marginals.fixed)
class(output$marginals.fixed)
post.beta0 = output$marginals.fixed[[1]]
post.beta1 = output$marginals.fixed$x

####beta1
plot(inla.smarginal(post.beta1),
     type="l",
     xlab="",
     ylab="",
     main=expression(paste("Post. marg. of ", beta[1])))

hpd.beta1 = inla.hpdmarginal(p=0.95,post.beta1)



# Prepare everything for the hpd polygon
xvalues = seq(hpd.beta1[1], hpd.beta1[2], length.out = 256)
yvalues = inla.dmarginal(xvalues,post.beta1)
xvalues = c(hpd.beta1[1],xvalues,hpd.beta1[2])
yvalues = c(0,yvalues,0)

polygon(xvalues, 
        yvalues, 
        col = "slateblue1", 
        border = 1)


####Intercept
plot(inla.smarginal(post.beta0),type="l",xlab="",ylab="",
     main=expression(paste("Post. marg. of ", beta[0])))
hpd.beta0 = inla.hpdmarginal(p=0.95,post.beta0)

# Prepare everything for the hpd polygon
xvalues = seq(hpd.beta0[1], hpd.beta0[2], length.out = 256)
yvalues = inla.dmarginal(xvalues,post.beta0)
xvalues = c(hpd.beta0[1],xvalues,hpd.beta0[2])
yvalues = c(0,yvalues,0)

polygon(xvalues, yvalues,
        col = "slateblue1", border = 1)


post.sigma = inla.tmarginal(function(x) sqrt(1/x), 
                            output$marginals.hyperpar[[1]]) 
#We use [[1]] instead of the long name

plot(inla.smarginal(post.sigma),
     type="l",
     xlab="",
     ylab="",
     main=expression(paste("Post. marg. of ", sigma)))
hpd.sigma = inla.hpdmarginal(p=0.95,post.sigma)

# Prepare everything for the hpd polygon
xvalues = seq(hpd.sigma[1], hpd.sigma[2], length.out = 256)
yvalues = inla.dmarginal(xvalues,post.sigma)
xvalues = c(hpd.sigma[1],xvalues,hpd.sigma[2])
yvalues = c(0,yvalues,0)

polygon(xvalues, yvalues, col = "slateblue1", border = 1)


inla.zmarginal(post.sigma)


output2 = inla(formula, 
               family="gaussian", 
               data=my.data,
               ## Changing the precision of NormDist used for estimate
               control.fixed=list(mean=0,
                                  prec=1,
                                  mean.intercept=0,
                                  prec.intercept=0.0001))


plot(inla.smarginal(post.beta0),
     type="l",
     xlab="",
     ylab="",
     main=expression(paste("Post. marg. of ", 
                           beta[0])))
lines(inla.smarginal(output2$marginals.fixed[[1]]),col=2)
abline(v=beta.0)
legend("topleft",col=c(1,2),lty=c(1,1),
       legend=c("Default","Normal(0,0.0001)"),box.lty=0)


output3 = inla(formula, 
               family="gaussian", 
               data=my.data,
               ## Changing for a NormDist(0,1) on the log of precision
               control.family=list(hyper=list(
                 prec=list(prior="gaussian",
                           param=c(0,1))))
)

output4 = inla(formula, 
               family="gaussian", 
               data=my.data,
               ## Changing for a Gamma(1,0.01)
               control.family=list(hyper=list(
                 prec=list(prior="loggamma",
                           param=c(1,0.01))))
)


post.sigma3 = inla.tmarginal(function(x) sqrt(1/x), 
                             output3$marginals.hyperpar[[1]]) 
post.sigma4 = inla.tmarginal(function(x) sqrt(1/x), 
                             output4$marginals.hyperpar[[1]]) 

plot(inla.smarginal(post.sigma),type="l",xlab="",
     ylab="",main=expression(paste("Post. marg. of ", sigma)))
lines(inla.smarginal(post.sigma3),col=2)
lines(inla.smarginal(post.sigma4),col=3)
abline(v=sigma)
legend("topright",col=c(1,2,3),lty=c(1,1,1),
       legend=c("Default","Normal(0,1)","logGamma(1,0.01)"),box.lty=0)

library(tidyverse)
smarginal0<-inla.smarginal(post.sigma) %>% 
  bind_rows
smarginal2<-inla.smarginal(post.sigma3) %>% 
  bind_rows()
smarginal3<-inla.smarginal(post.sigma4) %>% 
  bind_rows()

smarginal0 %>% 
  ggplot(aes(x = x , y = y, 
             col = "marginal for post.sigma"))+
  geom_line()+
  geom_line(data = smarginal2, 
            aes(x = x, y = y, 
                col = "marginal for post.sigma3"))+
  geom_line(data = smarginal3, aes(x = x , y = y , 
                                   col = "marginal for post.sigma4"))+
  theme_bw()+
  scale_color_viridis_d(name = "Marginal Posterior", option = "viridis")


