###########################################################
## To compare ordination recovery
## simulate from gllvm (K=1), cord (K=2), or compass (K=3),
## then fit cord, gllvm, nMDS and DCA,
## and compare true ordinations to fitted ordinations
## using procrustes distance
###########################################################


library(vegan)
library(CommEcol)
library(mvabund)
library(gllvm)
library(MASS)
library(ecoCopula)
library(mvtnorm)

# baseline data for simulation
data(spider)
y <- mvabund(spider$abun)

# simulaton methods
sim_methods <- c("gllvm","cord","compas")
B <- 50  #number of simulations
K=3 # which simulation method to use, need to cycle though 1,2,3

sim_method<-sim_methods[K]

#calculate 'true' lvs by fitting spider data to model

if(sim_method=="gllvm"){
  true.mod <- gllvm(y, family = "negative.binomial")
  true.ords <- true.mod$lvs
}

if(sim_method=="cord"){
  mvabund_sim=manyglm(y~1)
  true.mod=cord(mvabund_sim)
  true.ords <- true.mod$scores
}

if(sim_method=="compas"){
  bray <- vegdist(y^(1/4)); 
  true.ords <- metaMDS(bray,k=2,trace = FALSE)$points
}

pro.errors <- matrix(NA,B,4); 
colnames(pro.errors) <- c("gllvm","nMDS","DCA","cord")

for(b in 1:B) {
	
	#simulate
	if(sim_method=="gllvm"){
	  cw.y <- y
  	for(j in 1:ncol(y)) {
  	  for(i in 1:nrow(y)) {
  		  etas <- sum(c(1,true.mod$params$theta[j,])*c(true.mod$params$beta0[j],  true.mod$lv[i,])) #+ true.mod$site[i]  for sute effect, i don't have that
  		  cw.y[i,j] <- rnbinom(1, mu=exp(etas), size = true.mod$params$inv.phi[j])
  	  }
  	}
	}

	if(sim_method=="cord"){
    sig=true.mod$sigma[[1]]
    eta=mvabund_sim$coefficients[rep(1,nrow(y)),]
    phi.inv=1/t(as.matrix(mvabund_sim$phi))[rep(1,nrow(y)),]
    true.load=as.matrix(true.mod$loadings)
    Psi = diag(diag(sig - true.load %*% t(true.load)))
    cx.z=scale(true.ords%*%t(true.load) + rmvnorm(nrow(y),rep(0,ncol(y)),Psi))
  	cw.y <-qnbinom(pnorm(cx.z),mu=exp(eta),size=phi.inv)
	}

	if(sim_method=="compas"){
  	cw.y <- 1
  	class(cw.y)[1]="try-error"
  	while(class(cw.y)[1]=="try-error"){
  	  # values taken from first 2d example in compass
  	  cw.y <- try(compas(S=200, dims=2, am=c(2,2), 
  	                     beta.R=c(1,1), coords=true.ords, 
  	                     n.quanti=5, n.quali=0.1, add1=0.1))
  	  cw.y=try(cw.y[,1:12]) #can't control exactly how many species
  	}
	}

  # fit models
  
  # distance based methods
  fit.mds <- try(metaMDS(cw.y^(1/4),distance = "bray",trace = FALSE))
	fit.dca <- try(decorana(cw.y)) ## Fit DCA

	# copula
	mvabund_mod=try(manyglm(cw.y~1))
	fit.cord=try(cord(mvabund_mod))
	
	# hierarchical 
	fitlvm.nbsite <- gllvm(cw.y, family = "negative.binomial")

	# calculate procrustes dist
	p_lvm <- try(procrustes(true.ords, fitlvm.nbsite$lvs)$ss)
	if(class(p_lvm)[1]=="try-error"){p_lvm=NA}
	
	p_mds <- try(procrustes(true.ords, fit.mds$points)$ss)
	if(class(p_mds)[1]=="try-error"){p_mds=NA}
	
	p_dca <- try(procrustes(true.ords,fit.dca$rproj[,1:2])$ss)
	if(class(p_dca)[1]=="try-error"){p_dca=NA}
	
	p_cord <- try(procrustes(true.ords,fit.cord$scores)$ss)
	if(class(p_cord)[1]=="try-error"){p_cord=NA}
	
	pro.errors[b,] <- c(p_lvm,p_mds,p_dca,p_cord)

}
		
boxplot(pro.errors)
