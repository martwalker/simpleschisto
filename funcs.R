#########################################
## utility functions
########################################

funcs <- list(
# mating probability (monogamous, May 1977)
matingf =  function(W,k){
  
  W <- max(0, W)
  
  alpha= W/(k+W)
  
  integrand <- function(theta)
  {(1-cos(theta))/(((1+alpha*cos(theta))^(1+k)))}
  intsol <- integrate(integrand, lower=0, upper=2*pi)
  
  int <- intsol$value
  I= (((1-alpha)^(1+k))/(2*pi))
  I1= 1-I*int
  
  return(I1)
},

ddred = function(n, b){
    n^(b-1) 
},

prob = function(n, W, k, rho) {
  dnbinom(x= n, size = k, mu= W*rho)
}, 

# density-dependent reduction in worm fecundity
densdep = function(W, k, b, rho){
  
  max <- max(qnbinom(0.999, size=k, mu=I(W*rho)),1)
  
  probs <- dnbinom(seq(1,max), size = k, mu= W*rho)
    
  sumsol <- sum( funcs$ddred(1:max, b=b)*probs ) / (1-(1+(W*rho)/k)^(-k))

  return(as.numeric(sumsol))
},

netred = function(W, k, b) {
  funcs$densdep2(W=W,k=k,b=b)*funcs$matingf(W=W,k=k)
},

## this is the analytical solution for a different DD functional form
densdep2 = function(W, k, b) {
  A=W/(W+k)
  B=exp(-b)
  p0 = (1+W/k)^(-k)
  (1-A)^k/(B*(1-p0))*((1-A*B)^(-k)-1)
},

# treatment event function
eventfun0 = function(t, y, par) {
  with(as.list(c(y,par)), {
    ymat <- matrix(y, ncol=1, nrow=nw)
    ymat[,1] <- ymat[,1]*(1-coverage*efficacy)
    
    ymat[,1] <- max(ymat[,1], par["tol"])
    
    return(c(ymat))
  })
},

eventfun = function(t, y, par) {
  with(as.list(c(y,par)), {
    
    t1 <- seq(start.tx, start.tx+(n.tx-1)*freq.tx, freq.tx)
    
    if (any(abs(t-t1) < dt)) {
      funcs$eventfun0(t, y, par)
    } else {
      y
    }
    
  })
  
},

# main function for running the simple model
runmod = function(par, W0=2) {
  
  ## initialize dynamic variable environment
  e <<- new.env()
  e$kdyn <- as.vector(par["k"])
  e$kinf <- as.vector(par["k"])
  e$kpost <- as.vector(par["k"])
  e$Winf <- W0
  e$Wpost <- W0
  
  with(as.list(c(par)), {
    
    
    if(is.na(W0)) {
      y0 <- ((R0/rho)*mu2-mu2)/(R0^R0_weight*mu1*N1N2) 
    } else {
      y0 <- W0
    }
    
    
    if (dotx==0) {
      t1 <- seq(0,stop.t,by=dt)
      out<-lsoda(y=y0, times=t1, func=funcs$mod, parms=par)
    } else {
      stop.tx <- start.tx + (n.tx-1)*freq.tx
      t1 <- seq(0, stop.t, by=dt)
      out <- lsoda(y=y0, times=t1, func=funcs$mod, parms = par,
                   events=list(func=funcs$eventfun, 
                               time=seq(start.tx, stop.tx, freq.tx), 
                               par=par, root=F))
    }

    out <- as.data.frame(out)
    
    out
  })
}, 

derivs = function(W, par) {
  funcs$mod(t=0, y=W, par=par)[[1]]
},

findbreak = function(par)
{
  W<-cntrlpar$accbreak
  while(W<cntrlpar$maxbreak)
  {
    W0 <- W
    W1 <- W0+cntrlpar$accbreak
    dW0 <- funcs$derivs(W0, par=parameters)
    dW1 <- funcs$derivs(W1, par=parameters)
    i <- diff(sign(c(dW0, dW1)))
    if(i!=0) {
      return(W0+(W1-W0)/2)
      break
    }
    W <- W + cntrlpar$accbreak
  }
  return(NA)
},

findstable = function(Wbreak, par) {
  if(is.na(Wbreak)) {
    return(NA)
  } else {
    W <- Wbreak+cntrlpar$accbreak/2
    while(W<cntrlpar$maxendem) {
      W0 <- W
      W1 <- W0+cntrlpar$accendem
      dW0 <- funcs$derivs(W0, par=parameters)
      dW1 <- funcs$derivs(W1, par=parameters)
      i <- diff(sign(c(dW0, dW1)))
      if(i!=0) {
        return(W0+(W1-W0)/2)
        break
      }
      W <- W + cntrlpar$accendem
    }
    return(NA)
  }
},

findstable2 = function(par) {
  maxW <- max(as.numeric((parameters["R0"]/parameters["rho"])*parameters["mu2"]-parameters["mu2"]/
    (parameters["R0hs"]*parameters["mu1"]*parameters["N1N2"])),1)
  {
    W <- maxW-cntrlpar$accendem/2
    while(I(W - cntrlpar$accendem)>0) {
      W0 <- W
      W1 <- W0-cntrlpar$accendem
      dW0 <- funcs$derivs(W0, par=parameters)
      dW1 <- funcs$derivs(W1, par=parameters)
      i <- diff(sign(c(dW0, dW1)))
      if(i!=0) {
        return(W0-(W1-W0)/2)
        break
      }
      W <- W - cntrlpar$accendem
      #print(W)
      
    }
    return(NA)
  }
},

findstable3 =function(par) {
  tmp <- funcs$runmod(par = par, W0=NA)
  W <- tmp[nrow(tmp), "W"]
  RE <- tmp[nrow(tmp), "RE"]
  if (RE<cntrlpar$REthresh) {
    return(NA)
  } else {
    return(W)
  }
},

findbreak2 = function(par, Wstar)
{
  if(is.na(Wstar)) {
    return(NA)
  } else {
    W<-1.0E-12
    while(W<Wstar)
    {
      W0 <- W
      W1 <- W0+cntrlpar$accbreak
      dW0 <- funcs$derivs(W0, par=par)
      dW1 <- funcs$derivs(W1, par=par)
      i <- diff(sign(c(dW0, dW1)))
      if(i!=0) {
        return(W0+(W1-W0)/2)
        break
      }
    W <- W + cntrlpar$accbreak
    }
  }
  return(NA)
},

PRCC=function(outcome,covariates){
  
  rank_outcome<-rank(outcome)
  rank_covariates<-as.data.frame(apply(covariates, 2, rank))
  
  PRCC_out<-c()
  for (par in 1:ncol(covariates)){
    xx<-rank_covariates[,par]		  
    xy<-rank_covariates[,-par]	
    
    xj<-xx-predict(lm(xx~.,data=xy))
    yy<-rank_outcome-predict(lm(rank_outcome~.,data=xy))
    
    PRCC_out[par]<-cov(xj,yy)/(sqrt(var(xj)*var(yy)))	
  }
  names(PRCC_out)<-colnames(covariates)
  PRCC_out
}, 

mod = function(t,y,par){
  with(as.list(c(y,par)),{
    
    # initialise matrices & vectors
    ymat<- matrix(y, ncol=1, nrow=nw)
    dymat <- matrix(0, ncol=1, nrow=nw)
    
    ## prevents numerical errors when worm numbers very low
    ymat[,1]<-max(ymat[,1], par["tol"])
    
    ## mean worm burden
    W = sum(ymat[,1]) 
    
    # dynamic overdispersion parameter for the distribution of worms (NBD)
    kdyn <- as.vector(e$kdyn)
    kpost <- as.vector(e$kpost)
    kinf <- as.vector(e$kinf)
    Wpost <- as.vector(e$Wpost)
    Winf <- as.vector(e$Winf)
    
    ## density-dependent functions
    mmat <- funcs$matingf(W,k)
    dd <- funcs$densdep(W, k, b, rho)
    
    ## eggoutput
    epgout= rho*W*dd*mmat*a
    
    # prevalence of female infection
    prev_f = (1-(1+(W*rho)/kdyn)^-kdyn)
    
    # prevalence of mated female infection
    prev_mf = prev_f*mmat
    
    # prevalence of female or male infection
    prev = (1-(1+(W)/kdyn)^-kdyn)
    
    # prevalence of (detectable) heavy infection
    prev_h = (1-pnbinom((z/a)^(1/b), size = kdyn, mu = W*rho))*mmat
    
    # assymtery weighting
    R0hs <- R0^R0_weight
    
    prop_inf_snails <- mu1*R0hs*N1N2*W*mmat*dd/(mu1*R0hs*N1N2*W*mmat*dd+mu2)
    
    # Re is anverge number of worms (either sex) produced my single mated female 
    RE <- (R0/rho)*mmat*dd*(1-prop_inf_snails)
    
    
    # model ODEs    
    for (i in 1:nw) { 
      if (i==1) { 
        dymat[i,1] = ((R0/rho)*mu2*W*mmat*dd) / 
          (R0hs*N1N2*W*mmat*dd + mu2/(nw*mu1)) - (nw*mu1)*ymat[i,1]
      } else if (i>1) {
        dymat[i,1] = nw*mu1*ymat[i-1,1] - (nw*mu1)*ymat[i,1]
      }
    }
    
    if (dotx==1) {
      ## treatment time?
      tmp <- funcs$eventfun(t,W,par=par)
      if (tmp!=W) {
        kdyn <- kinf*W/( (1+kinf)*Winf - kinf*W)
        Wpost <- W
        kpost <- kdyn
      }   else if (t>start.tx) {
        kdyn <- W^2*(Winf-Wpost)^2/
          ((Winf^2/kinf)*(W-Wpost)^2 + (Wpost^2/kpost)*(W-Winf)^2)
      } else {
        kdyn<-kdyn
      }
      
    }
    
    ## store updates in environment
    e$kdyn <<- kdyn
    e$kinf <<- kinf
    e$kpost <<- kpost
    e$Wpost <<- Wpost
    e$Winf <<- Winf 
    
    
    
    
    return(list(rbind(dymat),
                W=W, E=epgout, prev=prev, prev_f=prev_f, prev_mf=prev_mf, prevh=prev_h, kdyn=kdyn, RE=RE, y=prop_inf_snails))
  })
}, 

model_out = function(par, W) {
  with(as.list(c(par)),{
  
  mmat <- funcs$matingf(W,k)
  dd <- funcs$densdep(W, k, b, rho)
  
  ## eggoutput
  epgout= rho*W*dd*mmat*a
  
  # prevalence of female infection
  prev_f = (1-(1+(W*rho)/k)^-k)
  
  # prevalence of mated female infection
  prev_mf = prev_f*mmat
  
  # prevalence of female or male infection
  prev = (1-(1+(W)/k)^-k)
  
  # prevalence of (detectable) heavy infection
  prev_h = (1-pnbinom((z/a)^(1/b), size = k, mu = W*rho))*mmat


  return(list(epg=epgout, prev_f=prev_f, prev_mf=prev_mf, prev=prev, prev_h=prev_h))
  })
}

)


