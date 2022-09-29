library(dplyr)        # A package for data manipulation
library(sf)           # Simple feature for R
library(spdep)        # Functions and tests for evaluating spatial patterns 
# and autocorrelation
library(tidyr)

library(INLA)         # Integrated Nested Laplace Approximation package
library(ggplot2)      # A package that implements the grammar of graphics, which is a term used to
# break up graphs into semantic components, such as geometries and layers.
library(viridis)      # A package providing color palettes 
library(patchwork)

# For tables in RMarkdown
library(knitr)
library(kableExtra)

RESP_DATA <- read.csv("RESP_DATA.csv", header=TRUE)

kable(RESP_DATA %>%
        group_by(year) %>%
        summarise(observed = sum(observed), expected=sum(expected)), booktabs = T, caption = "Hospital admissions by year") %>%
  kable_styling(bootstrap_options = "striped", full_width = F, position = "center")

GGHB <- st_read("GGHB.shp")

ggplot() + 
  geom_sf(data = GGHB, color = "blue", fill = "white") + 
  coord_sf() +    #axis limits and CRS
  theme_bw() +    # dark-on-light theme
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

GGHB_nb <- poly2nb(GGHB, snap=1000, queen=TRUE)
summary(GGHB_nb)

nb2INLA("GGHB.graph",GGHB_nb)
GGHB.adj <- paste(getwd(),"/GGHB.graph",sep="")

RESP_DATA %>% group_by(SP_ID) %>% 
  summarize(observed = sum(observed), 
            expected = sum(expected)) %>% 
  dplyr::rename(O = observed, E = expected) -> RESP_DATAagg


RESP_DATAagg %>% mutate(SMR = O/E) -> RESP_DATAagg

RESP_DATAagg$SMRcat <- cut(RESP_DATAagg$SMR, 
                           breaks=c(min(RESP_DATAagg$SMR), 
                                    0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 
                                    max(RESP_DATAagg$SMR)), include.lowest = T)

map_SMR <- left_join(GGHB, RESP_DATAagg, by = c("SP_ID" = "SP_ID"))

ggplot() + geom_sf(data = map_SMR, col = NA) + aes(fill = SMRcat) +
  theme_bw() + scale_fill_viridis_d() + 
  guides(fill=guide_legend(title="SMR")) 

ID<- seq(1,271)
formula_BYM2 <- O ~ f(ID, 
                      model="bym2", 
                      graph=GGHB.adj,
                      hyper=list(prec = list(
                        prior = "pc.prec",
                        param = c(0.5 / 0.31, 0.01)),
                        phi = list(
                          prior = "pc",
                          param = c(0.5, 2 / 3)))) 


sBYM.model <- inla(formula=formula_BYM2, 
                   family="poisson", 
                   data=RESP_DATAagg, 
                   E=E, 
                   control.compute=list(dic=TRUE, waic=TRUE))

#Relative risks
RR_sBYM<-c()
for(i in 1:271){
  RR_sBYM[i] <- inla.emarginal(function(x) exp(x), 
                               sBYM.model$marginals.random$ID[[i]])
}

#Posterior probabilities
RR_sBYM_marg <- sBYM.model$marginals.random$ID[1:271]
PP_sBYM <- lapply(RR_sBYM_marg, function(x) {1-inla.pmarginal(0,x)})

resRR_PP <- data.frame(resRR=RR_sBYM, 
                       PP=unlist(PP_sBYM),
                       SP_ID=RESP_DATAagg[,1])

resRR_PP$resRRcat <- cut(resRR_PP$resRR, 
                         breaks=c(min(resRR_PP$resRR), 
                                                  0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 
                                                  max(resRR_PP$resRR)),include.lowest = T)

# breakpoints
resRR_PP$PPcat <- cut(resRR_PP$PP, c(0, 0.2, 0.8, 1.00), include.lowest = TRUE)

map_RR_PP <- left_join(GGHB, resRR_PP, by = c("SP_ID" = "SP_ID"))

ggplot() + geom_sf(data = map_RR_PP) + aes(fill = resRRcat) +
  theme_bw() + scale_fill_brewer(palette = "PuOr") + 
  guides(fill=guide_legend(title="RR")) + ggtitle("RR Spatial model") + 
  theme(text = element_text(size=15), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(), plot.title = element_text(size = 12, face = "bold")
  )-> p1

ggplot() + geom_sf(data = map_RR_PP) + aes(fill = PPcat) +
  theme_bw() +
  scale_fill_viridis(
    option = "plasma", name="PP",
    discrete = T,
    direction = -1,
    guide = guide_legend(
      title.position = 'top',
      reverse = T
    )) +  ggtitle("PP Spatial model") + theme(text = element_text(size=15), 
                                              axis.text.x = element_blank(), 
                                              axis.text.y = element_blank(), plot.title = element_text(size = 12, face = "bold")
    ) -> p2

p1|p2

sBYM.model$summary.hyperpar

#Join the data with the shapefile so the order of the shapefile is maintained.  
RESP_DATA_ST <- left_join(GGHB, RESP_DATA, by="SP_ID")
#Rename the columns of Observed and Expected as we did before
RESP_DATA_ST <- RESP_DATA_ST  %>% dplyr::rename(O = observed, E = expected)
#Create the ID for year (time)
RESP_DATA_ST$ID.time <- RESP_DATA_ST$year - 2006
#Create the ID for space
RESP_DATA_ST$ID.space <- rep(seq(1,271),each=5)

formula_ST_noint <- O ~ f(ID.space, model="bym2", graph=GGHB.adj,
                          hyper=list(prec = list(
                            prior = "pc.prec",
                            param = c(0.5 / 0.31, 0.01)),
                            phi = list(
                              prior = "pc",
                              param = c(0.5, 2 / 3)))) + f(ID.time,model="rw1", hyper=list(prec = list(
                                prior = "pc.prec",
                                param = c(0.5 / 0.31, 0.01))))

stBYM.model <- inla(formula=formula_ST_noint, family="poisson", data=RESP_DATA_ST, E=E, control.compute=list(dic=TRUE, waic=TRUE))


#Spatial Relative risks
RR_stBYM<-c()
for(i in 1:271){
  RR_stBYM[i] <- inla.emarginal(function(x) exp(x), 
                                stBYM.model$marginals.random$ID.space[[i]])
}
#Posterior probabilities (for spatial RR)
RR_stBYM_marg <- stBYM.model$marginals.random$ID.space[1:271]
tRR_stBYM_marg <- stBYM.model$marginals.random$ID.time[1:5]
PP_stBYM <- lapply(RR_stBYM_marg, function(x) {1-inla.pmarginal(0,x)})
tPP_stBYM <- lapply(tRR_stBYM_marg, function(x) {1-inla.pmarginal(0,x)})

#Temporal Relative risks and CI95, and Posterior Probability and CI95
RR_stRW_RR<-c()
RR_stRW_lo<-c()
RR_stRW_hi<-c()

# tPP_stBYM %>% 
#   unlist() %>% plot()

# PP_stRW_RR<-c()
# PP_stRW_lo<-c()
# PP_stRW_hi<-c()

for(i in 1:5){
  #Posterior mean
  RR_stRW_RR[i] <- inla.emarginal(function(x) exp(x), 
                                  stBYM.model$marginals.random$ID.time[[i]])
  # PP_stRW_RR[i] <- inla.emarginal(function(x) exp(x), 
  #                                 stBYM.model$marginals.random$ID.space[[i]])
  #2.5% quantile 
  RR_stRW_lo[i] <- inla.qmarginal(0.025,inla.tmarginal(function(x) exp(x), stBYM.model$marginals.random$ID.time[[i]]))
  # PP_stRW_lo[i] <- inla.qmarginal(0.025,inla.tmarginal(function(x) exp(x), stBYM.model$marginals.random$ID.space[[i]]))
  #97.5% quantile 
  RR_stRW_hi[i] <- inla.qmarginal(0.975, inla.tmarginal(function(x) exp(x), stBYM.model$marginals.random$ID.time[[i]]))
  # PP_stRW_hi[i] <- inla.qmarginal(0.975, inla.tmarginal(function(x) exp(x), stBYM.model$marginals.random$ID.space[[i]]))
}

RR_stRW <- data.frame(RR=RR_stRW_RR,
                      low=RR_stRW_lo,
                      high=RR_stRW_hi)

# PP_stRW <- data.frame(PP=PP_stRW_RR, 
#                       low=PP_stRW_lo, 
#                       high=PP_stRW_hi)

ggplot(RR_stRW, aes(seq(2007,2011), RR)) + 
  geom_line() + 
  ggtitle("ST model No Int") + 
  geom_ribbon(aes(ymin=low,ymax=high), alpha=0.2) + 
  labs(x="year")-> Temp1
Temp1

ggplot(PP_stRW, aes(seq(2007,2011), PP)) + 
  geom_line() + 
  ggtitle("ST model No Int") + 
  geom_ribbon(aes(ymin=low,ymax=high), alpha=0.2) + 
  labs(x="year")-> Temp2
Temp2

resRR_PP_st <- data.frame(resRR=RR_stBYM, 
                          PP=unlist(PP_stBYM),
                          SP_ID=RESP_DATAagg[,1])
# breakpoints
resRR_PP_st$resRRcat <- cut(resRR_PP_st$resRR, 
                            breaks=c(min(resRR_PP_st$resRR), 
                                                        0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 
                                                        max(resRR_PP_st$resRR)),include.lowest = T)

resRR_PP_st$PPcat <- cut(resRR_PP_st$PP, c(0, 0.2, 0.8, 1.00), include.lowest = TRUE)

map_RR_ST <- left_join(GGHB, resRR_PP_st, by = c("SP_ID" = "SP_ID"))

ggplot() + geom_sf(data = map_RR_ST) + aes(fill = resRRcat) +
  theme_bw() + scale_fill_brewer(palette = "PuOr") + 
  guides(fill=guide_legend(title="RR")) +  ggtitle("RR ST model") +
  theme(text = element_text(size=15), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(), plot.title = element_text(size = 12, face = "bold")
  ) -> p3

ggplot() + geom_sf(data = map_RR_ST) + aes(fill = PPcat) +
  theme_bw() +
  scale_fill_viridis(
    option = "plasma",
    name = "PP ST model",
    discrete = T,
    direction = -1,
    guide = guide_legend(
      title.position = 'top',
      reverse = T
    )) +  ggtitle("PP ST model") + theme(text = element_text(size=15), 
                                         axis.text.x = element_blank(), 
                                         axis.text.y = element_blank(), plot.title = element_text(size = 12, face = "bold")
    )-> p4

(p1|p2) / (p3|p4)

RESP_DATA_ST$ID.space.time <- seq(1,dim(RESP_DATA_ST)[1])
formula_ST_intI <- O ~ f(ID.space, model="bym2", graph=GGHB.adj,
                         hyper=list(prec = list(
                           prior = "pc.prec",
                           param = c(0.5 / 0.31, 0.01)),
                           phi = list(
                             prior = "pc",
                             param = c(0.5, 2 / 3)))) + 
  f(ID.time,model="rw1", hyper=list(prec = list(
    prior = "pc.prec",
    param = c(0.5 / 0.31, 0.01))))+
  f(ID.space.time,model="iid", hyper=list(prec = list(
    prior = "pc.prec",
    param = c(0.5 / 0.31, 0.01))))


stIntI.BYM.model <- inla(formula=formula_ST_intI, 
                         family="poisson", 
                         data=RESP_DATA_ST, 
                         E=E, 
                         control.compute=list(dic=TRUE, 
                                              waic=TRUE))

#Spatial Relative risks
RR_stIntI.BYM<-c()
for(i in 1:271){
  RR_stIntI.BYM[i] <- inla.emarginal(function(x) exp(x), 
                                     stIntI.BYM.model$marginals.random$ID.space[[i]])
}
#Posterior probabilities (for spatial RR)
RR_stIntI.BYM_marg <- stIntI.BYM.model$marginals.random$ID.space[1:271]
tRR_stIntI.BYM_marg <- stIntI.BYM.model$marginals.random$ID.time[1:5]
stRR_stIntI.BYM_marg <- stIntI.BYM.model$marginals.random$ID.space.time[1:1355]
PP_stIntI.BYM <- lapply(RR_stIntI.BYM_marg, function(x) {1-inla.pmarginal(0,x)})
tPP_stIntI.BYM <- lapply(tRR_stIntI.BYM_marg, function(x) {1-inla.pmarginal(0,x)})
stPP_stIntI.BYM <- lapply(stRR_stIntI.BYM_marg, function(x) {1-inla.pmarginal(0,x)})

#Temporal Relative risks and CI95
RR_stIntI.RW_RR<-c()
RR_stIntI.RW_lo<-c()
RR_stIntI.RW_hi<-c()

for(i in 1:5){
  #Posterior mean
  RR_stIntI.RW_RR[i] <- inla.emarginal(function(x) exp(x), 
                                       stIntI.BYM.model$marginals.random$ID.time[[i]])
  #2.5% quantile 
  RR_stIntI.RW_lo[i] <- inla.qmarginal(0.025,inla.tmarginal(function(x) exp(x), stIntI.BYM.model$marginals.random$ID.time[[i]]))
  #97.5% quantile 
  RR_stIntI.RW_hi[i] <- inla.qmarginal(0.975, inla.tmarginal(function(x) exp(x), stIntI.BYM.model$marginals.random$ID.time[[i]]))
}

RR_stIntI.RW<- data.frame(RR=RR_stIntI.RW_RR,
                          low=RR_stIntI.RW_lo,
                          high=RR_stIntI.RW_hi)

ggplot(RR_stIntI.RW, aes(seq(2007,2011), RR)) + geom_line() + ggtitle("ST model Int I") + geom_ribbon(aes(ymin=low,ymax=high), alpha=0.2) + labs(x="year")->Temp2
Temp1 | Temp2

resRR_PP_stIntI <- data.frame(resRR=RR_stIntI.BYM, 
                              PP=unlist(PP_stIntI.BYM),
                              SP_ID=RESP_DATAagg[,1])
# breakpoints
resRR_PP_stIntI$resRRcat <- cut(resRR_PP_stIntI$resRR, 
                                breaks=c(min(resRR_PP_stIntI$resRR), 
                                                                0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 
                                                                max(resRR_PP_stIntI$resRR)),include.lowest = T)

resRR_PP_stIntI$PPcat <- cut(resRR_PP_stIntI$PP,
                             c(0, 0.2, 0.8, 1.00), 
                             include.lowest = TRUE)

map_RR_ST.IntI <- left_join(GGHB, 
                            resRR_PP_stIntI, 
                            by = c("SP_ID" = "SP_ID"))

resRR_PP_stIntI <- data.frame(resRR=RR_stIntI.BYM, 
                                                                                                           PP=unlist(PP_stIntI.BYM),
                                                                                                           SP_ID=RESP_DATAagg[,1])
# breakpoints
resRR_PP_stIntI$resRRcat <- cut(resRR_PP_stIntI$resRR, breaks=c(min(resRR_PP_stIntI$resRR), 
                                                                0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 
                                                                max(resRR_PP_stIntI$resRR)),include.lowest = T)

resRR_PP_stIntI$PPcat <- cut(resRR_PP_stIntI$PP, c(0, 0.2, 0.8, 1.00), include.lowest = TRUE)

map_RR_ST.IntI <- left_join(GGHB, resRR_PP_stIntI, by = c("SP_ID" = "SP_ID"))

ggplot() + geom_sf(data = map_RR_ST.IntI) + aes(fill = resRRcat) +
  theme_bw() + scale_fill_brewer(palette = "PuOr") + 
  guides(fill=guide_legend(title="RR")) +  ggtitle("RR ST model Int I") +
  theme(text = element_text(size=15), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(), plot.title = element_text(size = 12, face = "bold")
  ) -> p5

ggplot() + geom_sf(data = map_RR_ST.IntI) + aes(fill = PPcat) +
  theme_bw() +
  scale_fill_viridis(
    option = "plasma",
    name = "PP ST model Int I",
    discrete = T,
    direction = -1,
    guide = guide_legend(
      title.position = 'top',
      reverse = T
    )) +  ggtitle("PP ST model Int I") + theme(text = element_text(size=15), 
                                               axis.text.x = element_blank(), 
                                               axis.text.y = element_blank(), plot.title = element_text(size = 12, face = "bold")
    )-> p6

(p1|p2) / (p3|p4) / (p5|p6)


RESP_DATA_ST$intI<-stIntI.BYM.model$summary.random$ID.space.time$mean
RESP_DATA_ST$intI_cat <- cut(RESP_DATA_ST$intI,  
                             breaks=c(-1,-0.05, -0.01, 0.01, 0.05, 1),include.lowest = T)
ggplot() +
  geom_sf(data = RESP_DATA_ST, aes(fill = intI_cat))+ theme_bw() +  scale_fill_brewer(palette = "PuOr") + 
  guides(fill=guide_legend(title=NULL)) + 
  theme(text = element_text(size=20), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank()) +
  facet_wrap(~ year, ncol = 3, labeller=labeller(ID.year=c("1"="2007","2"="2008","3"="2009","4"="2010","5"="2011"))) +
  labs("")

dat.hyper2 <- 
  round(
    data.frame(median = stIntI.BYM.model$summary.hyperpar[,4],
               LL = stIntI.BYM.model$summary.hyperpar[,3], 
               UL = stIntI.BYM.model$summary.hyperpar[,5]),
    digits = 3)

row.names(dat.hyper2) <- 
  rownames(stIntI.BYM.model$summary.hyperpar)

knitr::kable(dat.hyper2, 
             caption = "Posterior median and 95% CrI of hyperparameters.") %>%  
  kable_styling(bootstrap_options = "striped", 
                full_width = F, 
                position = "center")

dat.WAIC <- data.frame(model = c("Spatial", "SpatTemp no int", "SpatTemp typeI"), 
                       WAIC = round(c(sBYM.model$waic$waic, stBYM.model$waic$waic, stIntI.BYM.model$waic$waic))
)

row.names(dat.WAIC) <- NULL

knitr::kable(dat.WAIC, 
             caption = "WAIC of the fifferent models") %>%  
  kable_styling(bootstrap_options = "striped", 
                full_width = F, 
                position = "center")
