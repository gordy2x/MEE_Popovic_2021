###########################################################
## To compare speed
## simulate from gllvm 
## then fit cord, gllvm, nMDS and DCA, and time
###########################################################


# library(gllvm)
library(ecoCopula)
library(reshape2)
library(tictoc)
library(vegan)
library(MASS)

#data
data(spider)
abund <- mvabund(spider$abund)
X=spider$x

#### fit gllvm
dat_lv=gllvm::gllvm(abund, family = "negative.binomial")

# simulation 
pars=expand.grid(K=2^(0:4),bigP=c(TRUE,FALSE),sims=1:20)

res=data.frame(N=NULL,P=NULL,elapsed_gllvm=NULL,elapsed_cord=NULL,
               elapsed_metaMDS=NULL,elapsed_decorana=NULL)

for(isim in 1:nrow(pars)){
  K=pars$K[isim]
  bigP=pars$bigP[isim]
  
  #if we want more species, make wide
  if(bigP){
    Y=matrix()
    for(i in 1:K){
      Y <- cbind(Y,simulate(dat_lv))
    }
    Y=Y[,-1]
  }else{ #long
    Y <- simulate(dat_lv,nsim = K)
  }
  N=nrow(Y)
  P=ncol(Y)
  
  # Fit models and time
  
  #gllvm
  tic("gllvm")
  fitl0 <- try(gllvm::gllvm(Y, family = "negative.binomial"))
  times_gl=toc(quiet = TRUE)
  elapsed_gllvm=times_gl$toc-times_gl$tic
  tic.clearlog()
  
  tic("cord")
  spider_mod=stackedsdm(Y,formula_X=~1,data=X)
  spid_lv=try(cord(spider_mod))
  times_co=toc(quiet = TRUE)
  elapsed_cord=times_co$toc-times_co$tic
  tic.clearlog()
  
  tic("metaMDS")
  bray <- vegdist(Y^(1/4))
  fit.mds <- try(vegan::metaMDS(bray,k=2,trace = FALSE)) ## Fit MDS
  times_MDS=toc(quiet = TRUE)
  elapsed_metaMDS=times_MDS$toc-times_MDS$tic
  tic.clearlog()
  
  tic("decorana")
  fit.dca <- try(decorana(Y)) ## Fit DCA
  times_dec=toc(quiet = TRUE)
  elapsed_decorana=times_dec$toc-times_dec$tic
  tic.clearlog()

  res=rbind(res,data.frame(N=rep(N,4),P=rep(P,4),
                           elapsed=round(c(elapsed_gllvm,elapsed_cord,elapsed_metaMDS,elapsed_decorana),2)+0.01,
                           method=c("gllvm","cord","metaMDS","decorana")))
  
  print(c(N,P))
  }

library(dplyr)
library(ggplot2)
res %>%  
  filter(N==28) %>% 
  ggplot(aes(x=P,y=elapsed+1e-14,fill=method,color=method, shape=method))+
  geom_point()+
  scale_y_continuous(trans="log10")+
  scale_x_continuous(breaks=sort(unique(res$P)))

res %>%  
  filter(P==12) %>% 
  ggplot(aes(x=N,y=elapsed+1e-14,color=method,fill=method, shape=method))+
  geom_point()+
  scale_y_continuous(trans="log10")+
  scale_x_continuous(breaks=sort(unique(res$N)))
