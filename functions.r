pacman::p_load(segmented, igraph, RSpectra, locfit, tidyverse, doParallel, broom, vegan, Matrix)
library(igraph)
library(iGraphMatch)
registerDoParallel(detectCores()-1)

##generate RDPG from given latent position X
rdpg.sample <- function(X, rdpg_scale=FALSE) {
  P <- X %*% t(X)
  if (rdpg_scale) {
    P <- scales::rescale(P)
  }
  n <-  nrow(P)
  U <- matrix(0, nrow = n, ncol = n)
  U[col(U) > row(U)] <- runif(n*(n-1)/2)
  U <- (U + t(U))
  diag(U) <- runif(n)
  A <- (U < P) + 0 ;
  diag(A) <- 0
  return(graph_from_adjacency_matrix(A,"undirected"))
}

##ASE for a network A with embedding dimension d
full.ase <- function(A, d, diagaug=TRUE, doptr=FALSE) {
  #    require(irlba)
  
  # doptr
  if (doptr) {
    g <- ptr(A)
    A <- g[]
  } else {
    A <- A[]
  }
  
  # diagaug
  if (diagaug) {
    diag(A) <- rowSums(A) / (nrow(A)-1)
  }
  
  A.svd <- svds(A,k=d)
  Xhat <- A.svd$u %*% diag(sqrt(A.svd$d))
  Xhat.R <- NULL
  
  if (!isSymmetric(A)) {
    Xhat.R <- A.svd$v %*% diag(sqrt(A.svd$d))
  }
  
  return(list(eval=A.svd$d, Xhat=Matrix(Xhat), Xhat.R=Xhat.R))
}

## Generate latent positions process from Model3 Constant drift with jump probability first-order changepoint 
doSim_shuffle <- function(n=300, tmax=40, delta=0.1, p=0.4, q=0.9, tstar=20,startpoint=0.1)
{
  glist <- NULL
  Xlist <- NULL
  Xt <- matrix(0,n,(tmax+1))
  
  for (t in 2:(tstar+1)) {
    tmp <- runif(n) < p
    Xt[,t] <- Xt[,t] + Xt[,t-1]
    Xt[tmp,t] <- Xt[tmp,t] + delta
  }
  
  for (t in (tstar+2):(tmax+1)) {
    tmp <- runif(n) < q
    Xt[,t] <- Xt[,t] + Xt[,t-1]
    Xt[tmp,t] <- Xt[tmp,t] + delta
  }
  Xt <- Xt[,-1]
  Xt <- Xt+startpoint
  
  df <- tibble(time=1:tmax) %>%
    mutate(Xt = map(time, function(x) matrix(Xt[,x],n,1)  )) %>%
    mutate(g = map(Xt, ~rdpg.sample(.))) %>%
    mutate(shuffle_g =map(g,~shuffle_graph(.)) )
  
  df
}


shuffle_graph <- function(A){
  G=as.matrix(A)
  permu_vec=sample(1:n)
  random_perm=diag(n)[permu_vec,]
  permu_G=as.matrix(random_perm%*%G%*%t(random_perm))
  G_graph=graph_from_adjacency_matrix(permu_G,mode ="undirected")
  
  return(G_graph)
}

graph_mathing <- function(stand,mess){
  
  G1=as.matrix( stand )
  G2=as.matrix( mess )
  
  gm=gm(A=stand,B=mess,start = "rds", max_iter = 20)
  perm=diag(length(gm[,2]))[gm[,2],]
  
  new_graph=as.matrix(perm%*%as.matrix( mess )%*% t(perm))
  
  print(c("Pre_Frob=",sqrt(sum((G1-G2)^2)),"Post_Frob=", sqrt(sum((G1-new_graph)^2)) ) )
  
  return( graph_from_adjacency_matrix(new_graph,mode = 'undirected') )
}


#graph_mathing(df$g[[2]],df$g[[3]])

## Procruste/i.e gain W in the estimated d_MV distance. Note in our case latent positions are 1d so this is trivial with w=1or -1
procrustes2 <- function(X, Y) {
  tmp <- t(X) %*% Y
  tmp.svd <- svd(tmp)
  W <- tmp.svd$u %*% t(tmp.svd$v)
  newX <- X %*% W
  return(list(newX = newX, error = norm(newX-Y, type="F"), W = W))
}


## Get distance matrix
getD <- function(Xlist, k=0, etype="proc") {
  m <- length(Xlist)
  if (k==0) {
    ind <- 1:n
  } else {
    ind <- which(Yhat==k)
  }
  comb <- combn(m,2)
  Dout <- foreach (k = 1:ncol(comb), .combine='rbind') %dopar% {
    i <- comb[1,k]
    j <- comb[2,k]
    #cat("i = ", i, ", j = ", j, "\n")
    
    if (etype == "proc") {
      Xhati <- Xlist[[i]][ind,] # 32277 x Khat
      Xhatj <- Xlist[[j]][ind,]
      proc <- procrustes2(as.matrix(Xhati), as.matrix(Xhatj))
      Xhati <- Xhati %*% proc$W
    } else {
      Xhati <- Xlist[[i]][ind,] # 32277 x Khat
      Xhatj <- Xlist[[j]][ind,]
    }
    
    D <- norm(Xhati - Xhatj, type="2")^2/n
    tibble(i=i, j=j, D=D)
  }
  D2 <- matrix(0,m,m)
  D2[t(comb)] <- Dout$D
  D2 <- (D2 + t(D2)) / 1
  #as.dist(D2)
  D2 <- sqrt(D2)
  D2
}

## Apply CMDS on distance matrix 
doMDS <- function(D, doplot=TRUE)
{
  tmax <- m <- nrow(D)
  mds <- cmdscale(D, m-1)
  df.mds <- tibble(ind=1:tmax, time=sprintf("%2d",1:tmax), x=mds[,1], y=mds[,2], z=mds[,3], w=mds[,4])
  
  if (doplot) {
    plot(apply(mds,2,sd), type="b", main="", xlab="dimension", ylab="column stdev")
    
    p1 <- df.mds %>% ggplot(aes(x=ind, y=x, color=time, group=1)) +
      geom_point(size=3) + geom_line() +
      geom_vline(xintercept = tstar, linetype="dashed") +
      theme(legend.position = "none") + labs(x="time",y="mds1") #+
    print(p1)
    
    p1 <- df.mds %>% ggplot(aes(x=ind, y=y, color=time, group=1)) +
      geom_point(size=3) + geom_line() +
      geom_vline(xintercept = tstar, linetype="dashed") +
      theme(legend.position = "none") + labs(x="time",y="mds2") #+
    print(p1)
    
    p1 <- df.mds %>% ggplot(aes(x=ind, y=z, color=time, group=1)) +
      geom_point(size=3) + geom_line() +
      geom_vline(xintercept = tstar, linetype="dashed") +
      theme(legend.position = "none") + labs(x="time",y="mds3") #+
    print(p1)
    
    p2 <- df.mds %>% ggplot(aes(x=x, y=y, color=time)) +
      geom_point(size=3) +
      geom_label_repel(aes(label=time), size=2) +
      theme(legend.position = "none") + labs(x="mds1",y="mds2") #+
    print(p2)
  }
  
  return(list(mds=mds, df.mds=df.mds))
}

#apply ISOMAP on the CMDS result with chosen dimension mdsd from CMDS step defaultly it always embeds to 1 
doIso <- function(mds, mdsd=2, isod=1, doplot=F)
{
  df.iso <- NULL
  dis <- vegdist(mds[,1:mdsd,drop=F], "euclidean")
  knn <- 1
  success <- FALSE
  while(!success) {
    tryCatch({
      iso = isomap(dis, k=knn, ndim=isod, path="shortest")$points
      success <- TRUE
    },
    error = function(e) {
      knn <<- knn + 1
    })
  }
  iso2 <- tibble(iso=iso[,1]) %>% mutate(i=1:nrow(mds), ind=df.mds$ind, time=df.mds$time, knn=knn)
  df.iso <- rbind(df.iso, cbind(iso2, mdsd=mdsd))
  df.iso <- df.iso %>% group_by(mdsd) %>% mutate(iso = if(iso[1] > 0) {-iso} else {iso}) %>% ungroup()
  
  if (doplot) {
    p <- df.iso %>% filter(mdsd==mdsd) %>%
      ggplot(aes(x=ind, y=iso, color=time, group=1)) +
      geom_point(size=3) + geom_line() +
      theme(legend.position = "none") +
      labs(x="time", y="isomap embedding") +
      # scale_x_date(breaks = scales::breaks_pretty(8), labels=label_date_short()) +
      theme(axis.text.x=element_text(hjust=0.7))
    # theme(axis.text.x = element_text(size = 12, angle = 90, vjust = 0.3),
    #       axis.text.y = element_text(size = 12),
    #       axis.title = element_text(size = 14, face="bold"))
    #    p <- p + scale_x_date(breaks = scales::breaks_pretty(8), labels=label_date_short())
    print(p)
    
    df.isok <- df.iso %>% filter(mdsd==mdsd) #%>% mutate(date2 = format(ymd(paste0(date,"-01")),"%m/%y"))
    row.names(df.isok) <- df.isok$time
    fit <- lm(iso ~ i, data=df.isok)
    # print(tidy(fit))
    # print(glance(fit))
    myfor <- augment(fit)
    myfor2 <- myfor %>% mutate(date=.rownames,
                               ranks = rank(.sigma),
                               mycol=sprintf("%2d",rank(.fitted)))
    p <- myfor2 %>%
      ggplot(aes(.fitted, .resid)) +
      geom_point(aes(color=mycol)) +
      geom_hline(yintercept = 0, linetype="dashed", color="grey") +
      geom_smooth(method="loess", se=FALSE) +
      labs(x="Fitted Values", y="Residuals") +
      theme(legend.position = "none",
            axis.title = element_text(size=14, face="bold"))
    p <- p + geom_label_repel(aes(label=date), data=myfor2 %>% filter(ranks %in% 1:3))
    print(p)
  }
  
  return(df.iso)
}

## This is another slope change point algorithm called segmented that is not used in the paper 
break_point_dection=function(D,k){
  tmax <- nrow(D)
  df.mds <- doMDS(D,doplot = F)
  mds <- df.mds$mds
  df.iso <- doIso(mds, mdsd=k)
  x=1:tmax/tmax
  y1=df.iso$iso
  os1<-segmented(lm(y1~x),psi=c(0.2))
  result=as.data.frame(matrix(0,1,3))
  result[1,1:4]=c((os1$psi[1,2]-1.95*os1$psi[1,3])*tmax,(os1$psi[1,2])*tmax,(os1$psi[1,2]+1.95*os1$psi[1,3])*tmax,os1$psi[1,3]*tmax)
  return(result)
  ## this result returns you the confidence interval of point estimation and the standard deviation
}


## Implementation of the 3rd step in Algorithm 2 by recasting it as a linear programming problem.
## This function will return the objective function value Sk in the paper for a given change point t 
linf_cp=function(t,y,cp){
  n=length(t)
  nl=sum(t<cp)+1
  XL=matrix(1,nrow = nl,ncol=4)
  XL[,4]=0
  XL[,3]=t[1:nl]-cp
  
  #XL; y[1:nl]
  
  XL2=XL
  XL2[,2:3]=-XL[,2:3]
  
  #rbind(XL,XL2); c(y[1:nl],-y[1:nl])
  
  XR=matrix(1,nrow = n-nl,ncol = 4)
  XR[,3]=0
  XR[,4]=t[(nl+1):n]-cp
  
  XR2=XR
  XR2[,c(2,4)]=-XR[,c(2,4)]
  
  
  X=rbind(XL,XR,XL2,XR2)
  Y=c(y,-y)
  
  library(lpSolveAPI)
  lprec <- make.lp(0,4)
  set.objfn(lprec,c(1,0,0,0))
  for (i in 1:(nrow(X)) ) {
    add.constraint(lprec, X[i,], ">=", Y[i])
  }
  
  set.bounds(lprec, lower = c(0,-Inf,-Inf,-Inf), columns = c(1,2,3,4))
  ColNames <- c('Z', "alpha", "bl","br")
  dimnames(lprec)[[2]] <- ColNames
  solve(lprec)
  return(get.variables(lprec))
}


