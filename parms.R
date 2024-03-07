#########################################
## default parameters
########################################
parameters= c(R0=5,                   # basic reproduction number
              R0_weight = 0.5,         # asymmetry weighting - when 
              N1N2 = 1/10,            # host-to-snail population density
              mu1= 1/5,               # adult worm mortality rate (1/mu1=average lifespan)
              mu2= 1/(1/12),          # snail mortality rate (1/mu2=average snail lifespan)
              k= 0.5,                 # overdispersion adult worms among hosts
              coverage= 0.50,         # population treatment coverage
              efficacy= 0.94,         # efficacy of PZQ
              start.tx=1,             # treatment start time
              n.tx=25,                # number of treatments
              freq.tx=1,              # frequency of treatments
              nw=1,                   # adult worm structure (not used)
              dt=1/12,                # integration time step (fixed)
              stop.t=100,              # simulation stop time
              a= exp(3.24),           # maximum egg output from 1 Sh female (Neves et al 2021) - assumes Nmax = 300 (based on max observation from autopsy)
             # a= exp(1.69)           # maximum egg output for 1 Sm female (Neves et al 2021)
              b=0.25,                  # density-dependent constraint on fecundity for Sh (Neves et al 2021) - assumes Nmax = 300 (based on max observation from autopsy)
             # b = 0                  # density-dependent constraint on fecundity for Sm (Neves et al 2022) 
              z=50,                   # threshold egg density for heavy infection Sh
             # z = 400                # threshold eddg density for heavy infection Sm
              dotx=0,                 # toggle treatments on/off
              rho=0.5, 
              tol=1.0E-6)                 

cntrlpar = list(accbreak=0.0125, accendem=0.5, block=10, nsamp=5000, REthresh=0.99, sh=1)

if (cntrlpar$sh!=1) {
  parameters["a"] = exp(1.69)
  parameters["b"] = 0.26 # perhaps this should be 0 for no density dependence 
  parameters["z"] = 400
}
