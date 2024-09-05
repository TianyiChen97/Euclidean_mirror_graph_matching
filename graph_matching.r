pacman::p_load(segmented, igraph, irlba, locfit, tidyverse, doParallel, broom, vegan, Matrix)
library(igraph)
library(iGraphMatch)
registerDoParallel(detectCores()-1)
set.seed(2)
df.mds=NULL
n=200
p <- 0.8
q <- 0.2

m <- tmax <- 20

stp=0.2

delta <- (1-stp)/tmax
tstar <- tmax/2

df <- doSim_shuffle(n,tmax,delta,p,q,tstar,startpoint = stp)
df <- df %>% mutate(Xhat = map(g, function(x) full.ase(x,2)$Xhat[,1,drop=F]))%>% 
  mutate (num_e=map(g, ~sum(as.matrix(.x)))) %>%
  mutate(Xhat_shuffle = map(shuffle_g, function(x) full.ase(x,2)$Xhat[,1,drop=F])) %>%
  mutate(shuffle_g_GM = map2(shuffle_g, lag(shuffle_g),~ if (!is_null(.y)) graph_mathing(.y, .x) else .x)) %>%
  mutate(Xhat_shuffle_GM = map(shuffle_g_GM, function(x) full.ase(x,2)$Xhat[,1,drop=F]))

D2 <- getD(df$Xhat) 
D2_shuffle <- getD(df$Xhat_shuffle) 
D2_shuffle_GM <- getD(df$Xhat_shuffle_GM) 

df.mds <- doMDS(D2,doplot = T)
df.mds_shuffle <- doMDS(D2_shuffle,doplot = T)
df.mds_shuffle_GM <- doMDS(D2_shuffle_GM,doplot = T)

mds <- df.mds$mds
mds_shuffle <- df.mds_shuffle$mds
mds_shuffle_GM <- df.mds_shuffle_GM$mds

df.iso <- doIso(mds, mdsd=3)
df.iso_shuffle <- doIso(mds_shuffle, mdsd=3)
df.iso_shuffle_GM <- doIso(mds_shuffle_GM, mdsd=3)

par(mfrow=c(1,3))
plot(1:tmax,df.iso$iso)
plot(1:tmax,df.iso_shuffle$iso)
plot(1:tmax,df.iso_ shuffle_GM$iso)

