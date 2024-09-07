pacman::p_load(segmented, igraph, RSpectra, locfit, tidyverse, doParallel, broom, vegan, Matrix)
library(igraph)
library(iGraphMatch)

#library('remotes')
#install_version("Matrix", version='1.6-2' , INSTALL_opts = '--no-lock')


registerDoParallel(detectCores()-1)
set.seed(2)
df.mds=NULL
n=200
p <- 0.8
q <- 0.2

m <- tmax <- 10

stp=0.2

delta <- (1-stp)/tmax
tstar <- tmax/2

df <- doSim_shuffle(n,tmax,delta,p,q,tstar,startpoint = stp)
df <- df %>% mutate(Xhat = map(g, function(x) full.ase(x,2)$Xhat[,1,drop=F]))%>% 
  mutate (num_e=map(g, ~sum(as.matrix(.x)))) %>%
  mutate(Xhat_shuffle = map(shuffle_g, function(x) full.ase(x,2)$Xhat[,1,drop=F])) %>%
  mutate(shuffle_g_GM = map2(shuffle_g, lag(shuffle_g),~ if (!is_null(.y)) graph_mathing(.y, .x,20) else .x)) %>%
  mutate(Xhat_shuffle_GM = map(shuffle_g_GM, function(x) full.ase(x,2)$Xhat[,1,drop=F]))



D2 <- getD(df$Xhat) 
D2_shuffle <- getD(df$Xhat_shuffle) 
D2_shuffle_GM <- getD(df$Xhat_shuffle_GM) 

df.mds <- doMDS(D2,doplot = F)
df.mds_shuffle <- doMDS(D2_shuffle,doplot = F)
df.mds_shuffle_GM <- doMDS(D2_shuffle_GM,doplot = F)

mds <- df.mds$mds
mds_shuffle <- df.mds_shuffle$mds
mds_shuffle_GM <- df.mds_shuffle_GM$mds

df.iso <- doIso(mds, mdsd=2)
df.iso_shuffle <- doIso(mds_shuffle, mdsd=2)
df.iso_shuffle_GM <- doIso(mds_shuffle_GM, mdsd=2)

#par(mfrow=c(1,3))
#plot(1:tmax,df.iso$iso)
#plot(1:tmax,df.iso_shuffle$iso)
#plot(1:tmax,df.iso_shuffle_GM$iso)

linf_error=function(x){
  obf=NULL
  for (nk in 2:(tmax-1)) { ## find the point which minimize the obj func Sk, that is the change point 
  obf[nk]=linf_cp(1:tmax,x,nk)[1]
  }
  ecp=min(which(obf==min(obf[-1])))
  return((ecp-tstar)/tmax)
}

linf_error(df.iso$iso)

nmc=100
mm=c(10,16,20)
n=100

p <- 0.8
q <- 0.5
m <- tmax <- 10
stp=0.2

out <- foreach (mc = 1:nmc) %do% {
  tmp1 <- tmp2 <- tmp3 <- rep(0,length(mm))
  for(i in 1:length(mm)){
    m <- tmax <- mm[i]
    delta <- (1-stp)/tmax
    tstar <- tmax/2
    df <- doSim_shuffle(n,tmax,delta,p,q,tstar,startpoint = stp)
    df <- df %>% mutate(Xhat = map(g, function(x) full.ase(x,2)$Xhat[,1,drop=F]))%>% 
    mutate (num_e=map(g, ~sum(as.matrix(.x)))) %>%
    mutate(Xhat_shuffle = map(shuffle_g, function(x) full.ase(x,2)$Xhat[,1,drop=F])) %>%
    mutate(shuffle_g_GM = map2(shuffle_g, lag(shuffle_g),~ if (!is_null(.y)) graph_mathing(.y, .x,10) else .x)) %>%
    mutate(Xhat_shuffle_GM = map(shuffle_g_GM, function(x) full.ase(x,2)$Xhat[,1,drop=F]))

    D2 <- getD(df$Xhat)
    D2_shuffle <- getD(df$Xhat_shuffle) 
    D2_shuffle_GM <- getD(df$Xhat_shuffle_GM) 

    df.mds <- doMDS(D2,doplot = F)
    df.mds_shuffle <- doMDS(D2_shuffle,doplot = F)
    df.mds_shuffle_GM <- doMDS(D2_shuffle_GM,doplot = F)

    mds <- df.mds$mds
    mds_shuffle <- df.mds_shuffle$mds
    mds_shuffle_GM <- df.mds_shuffle_GM$mds

    df.iso <- doIso(mds, mdsd=2)
    df.iso_shuffle <- doIso(mds_shuffle, mdsd=2)
    df.iso_shuffle_GM <- doIso(mds_shuffle_GM, mdsd=2)


    tmp1[i]=linf_error(df.iso$iso)
    tmp2[i]=linf_error(df.iso_shuffle$iso)
    tmp3[i]=linf_error(df.iso_shuffle_GM$iso)

    print(c(mc,m,tmp1[i],tmp2[i],tmp3[i]))
  }
  list(tmp1,tmp2,tmp3)
}

save.image(file = paste0(getwd(), "/data.RData"))
