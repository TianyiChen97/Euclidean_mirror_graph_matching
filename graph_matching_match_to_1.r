
# nolint start
load(file = paste0(getwd(), "/data.RData"))
pacman::p_load(segmented, igraph, RSpectra, locfit, tidyverse, doParallel, broom, vegan, Matrix)
library(igraph)
library(iGraphMatch)

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
  mutate(shuffle_g_GM_alltoone = map(shuffle_g, ~graph_mathing(shuffle_g[[1]],.x,100))) %>%
  mutate(shuffle_g_GM_pairwise = purrr::accumulate(shuffle_g[-1],
                                          .f = function(acc, curr) graph_mathing(acc, curr, 50),
                                          .init = shuffle_g[[1]]))%>%  # Sequential matching with correction
  mutate(Xhat_shuffle_GM_alltoone = map(shuffle_g_GM_alltoone, function(x) full.ase(x,2)$Xhat[,1,drop=F]))%>%
  mutate(Xhat_shuffle_GM_pairwise = map(shuffle_g_GM_pairwise, function(x) full.ase(x,2)$Xhat[,1,drop=F]))

D2 <- getD(df$Xhat) 
D2_shuffle <- getD(df$Xhat_shuffle)
D2_shuffle_GM_alltoone <- getD(df$Xhat_shuffle_GM_alltoone)
D2_shuffle_GM_pairwise <- getD(df$Xhat_shuffle_GM_pairwise) 

par(mfrow=c(1,3))
df.mds <- doMDS(D2,doplot = T)
df.mds_shuffle <- doMDS(D2_shuffle,doplot = F)
df.mds_shuffle_GM_alltoone <- doMDS(D2_shuffle_GM_alltoone,doplot = F)
df.mds_shuffle_GM_pairwise <- doMDS(D2_shuffle_GM_pairwise,doplot = F)

mds <- df.mds$mds
mds_shuffle <- df.mds_shuffle$mds
mds_shuffle_GM_alltoone <- df.mds_shuffle_GM_alltoone$mds
mds_shuffle_GM_pairwise <- df.mds_shuffle_GM_pairwise$mds

df.iso <- doIso(mds, mdsd=1)
df.iso_shuffle <- doIso(mds_shuffle, mdsd=1)
df.iso_shuffle_GM_alltoone <- doIso(mds_shuffle_GM_alltoone, mdsd=1)
df.iso_shuffle_GM_pairwise <- doIso(mds_shuffle_GM_pairwise, mdsd=1)

par(mfrow=c(1,4))
plot(1:tmax,df.iso$iso)
plot(1:tmax,df.iso_shuffle$iso)
plot(1:tmax,df.iso_shuffle_GM_alltoone$iso)
plot(1:tmax,df.iso_shuffle_GM_pairwise$iso)


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

    for (dd in 1:3){
    df.iso <- doIso(mds, mdsd=dd)
    df.iso_shuffle <- doIso(mds_shuffle, mdsd=dd)
    df.iso_shuffle_GM <- doIso(mds_shuffle_GM, mdsd=dd)

    tmp1[dd]=linf_error(df.iso$iso)
    tmp2[dd]=linf_error(df.iso_shuffle$iso)
    tmp3[dd]=linf_error(df.iso_shuffle_GM$iso)


    }

    print(c(mc,m,tmp1[i],tmp2[i],tmp3[i]))
  }
  list(tmp1,tmp2,tmp3)
}

save.image(file = paste0(getwd(), "/data.RData"))



list(tmp1,tmp2,tmp3)


msehat1=Reduce('cbind', lapply(out, "[[", 1)) ## this will summarize tmp1 

sm_mse1=as.data.frame(matrix(0,3,4))
sm_mse1[,2]=apply(abs(msehat1)^2, 1, mean)
sm_mse1[,1]=apply(abs(msehat1)^2, 1, mean)-apply(abs(msehat1)^2, 1, sd)*1.96/sqrt(nmc)  
sm_mse1[,3]=apply(abs(msehat1)^2, 1, mean)+apply(abs(msehat1)^2, 1, sd)*1.96/sqrt(nmc) 
sm_mse1[,4]=mm

sm_mse

msehat2=Reduce('cbind', lapply(out, "[[", 2))

sm_mse2=as.data.frame(matrix(0,3,4))
sm_mse2[,2]=apply(abs(msehat2)^2, 1, mean)
sm_mse2[,1]=apply(abs(msehat2)^2, 1, mean)-apply(abs(msehat2)^2, 1, sd)*1.96/sqrt(nmc)  
sm_mse2[,3]=apply(abs(msehat2)^2, 1, mean)+apply(abs(msehat2)^2, 1, sd)*1.96/sqrt(nmc) 
sm_mse2[,4]=mm

sm_mse2


msehat3=Reduce('cbind', lapply(out, "[[", 3))

sm_mse3=as.data.frame(matrix(0,3,4))
sm_mse3[,2]=apply(abs(msehat3)^2, 1, mean)
sm_mse3[,1]=apply(abs(msehat3)^2, 1, mean)-apply(abs(msehat3)^2, 1, sd)*1.96/sqrt(nmc)  
sm_mse3[,3]=apply(abs(msehat3)^2, 1, mean)+apply(abs(msehat3)^2, 1, sd)*1.96/sqrt(nmc) 
sm_mse3[,4]=mm

sm_mse3


sm_mse1$type <- "True 1-1"
sm_mse2$type <- "shuffled"
sm_mse3$type <- "shuffled then GM"
sm_mse_all <- rbind(sm_mse1, sm_mse2, sm_mse3)

library(ggplot2)
plottt <- ggplot(sm_mse_all, aes(x=V4, y=V2, color=type, linetype=type)) + 
  geom_line() +
  geom_errorbar(aes(ymin=V1, ymax=V3)) +
  scale_x_continuous(breaks = mm) +
  labs(y='relative MSE', x='m', color='Type', linetype='Type') +
  theme(axis.text=element_text(size=25), axis.title=element_text(size=25, face="bold"))
print(plottt)






n=300
p <- 0.4
q <- 0.3
m <- tmax <- 10
stp=0.1
delta <- (1-stp)/tmax
tstar <- tmax/2
nmc=200

out_dd <- foreach (mc = 1:nmc) %do% {
  tmp1 <- tmp2 <- tmp3 <- rep(0,3)
  df <- doSim_shuffle(n,tmax,delta,p,q,tstar,startpoint = stp)
  df <- df %>% mutate(Xhat = map(g, function(x) full.ase(x,2)$Xhat[,1,drop=F]))%>% 
    mutate (num_e=map(g, ~sum(as.matrix(.x)))) %>%
    mutate(Xhat_shuffle = map(shuffle_g, function(x) full.ase(x,2)$Xhat[,1,drop=F])) %>%
    mutate(shuffle_g_GM = map(shuffle_g, ~graph_mathing(shuffle_g[[1]],.x,100))) %>%
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

  for (dd in 1:3){
      df.iso <- doIso(mds, mdsd=dd)
      df.iso_shuffle <- doIso(mds_shuffle, mdsd=dd)
      df.iso_shuffle_GM <- doIso(mds_shuffle_GM, mdsd=dd)

      tmp1[dd]=linf_error(df.iso$iso)
      tmp2[dd]=linf_error(df.iso_shuffle$iso)
      tmp3[dd]=linf_error(df.iso_shuffle_GM$iso)
    }
  print(c(mc,m,tmp1,tmp2,tmp3))
  list(tmp1,tmp2,tmp3)
}


#load(file = paste0(getwd(), "/data_0911.RData"))

msehat1=Reduce('cbind', lapply(out_dd, "[[", 1)) ## this will summarize tmp1 

sm_mse1=as.data.frame(matrix(0,3,4))
sm_mse1[,2]=apply(abs(msehat1)^2, 1, mean)
sm_mse1[,1]=apply(abs(msehat1)^2, 1, mean)-apply(abs(msehat1)^2, 1, sd)*1.96/sqrt(nmc)  
sm_mse1[,3]=apply(abs(msehat1)^2, 1, mean)+apply(abs(msehat1)^2, 1, sd)*1.96/sqrt(nmc) 
sm_mse1[,4]=1:3

sm_mse1
msehat2=Reduce('cbind', lapply(out_dd, "[[", 2))

sm_mse2=as.data.frame(matrix(0,3,4))
sm_mse2[,2]=apply(abs(msehat2)^2, 1, mean)
sm_mse2[,1]=apply(abs(msehat2)^2, 1, mean)-apply(abs(msehat2)^2, 1, sd)*1.96/sqrt(nmc)  
sm_mse2[,3]=apply(abs(msehat2)^2, 1, mean)+apply(abs(msehat2)^2, 1, sd)*1.96/sqrt(nmc) 
sm_mse2[,4]=1:3

sm_mse2

msehat3=Reduce('cbind', lapply(out_dd, "[[", 3))

sm_mse3=as.data.frame(matrix(0,3,4))
sm_mse3[,2]=apply(abs(msehat3)^2, 1, mean)
sm_mse3[,1]=apply(abs(msehat3)^2, 1, mean)-apply(abs(msehat3)^2, 1, sd)*1.96/sqrt(nmc)  
sm_mse3[,3]=apply(abs(msehat3)^2, 1, mean)+apply(abs(msehat3)^2, 1, sd)*1.96/sqrt(nmc) 
sm_mse3[,4]=1:3
sm_mse3

# Combine the data frames and add a new column 'type' to distinguish them
sm_mse1$type <- "True 1-1"
sm_mse2$type <- "shuffled"
sm_mse3$type <- "shuffled then GM"
sm_mse_all <- rbind(sm_mse1, sm_mse2, sm_mse3)

# Plot

library(ggplot2)
plottt <- ggplot(sm_mse_all, aes(x=V4, y=V2, color=type, linetype=type)) + 
  geom_line() +
  geom_errorbar(aes(ymin=V1, ymax=V3)) +
  scale_x_continuous(breaks = 1:3) +
  labs(y='relative MSE', x='MDS embedding dim d for the (d -> 1)-iso-mirror', color='Type', linetype='Type') +
  theme(axis.text=element_text(size=25), axis.title=element_text(size=25, face="bold"))


plottt <- plottt + labs(title = paste("p =", p, "q =", q,"m=",m, "n =", n, "nmc =", nmc, 'max_iter=', 100))
plottt <- plottt + theme(
  axis.text = element_text(size = 25),
  axis.title = element_text(size = 25, face = "bold"),
  legend.text = element_text(size = 25),
  legend.title = element_text(size = 25, face = "bold"),
  plot.title = element_text(size = 20, face = "bold") 
)



print(plottt)


save.image(file = paste0(getwd(), "/data_0912.RData"))
