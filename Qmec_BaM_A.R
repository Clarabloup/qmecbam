cat("\014") # Clear console
rm(list=ls())# Clean workspace

library(ggplot2) #graphiques
library(RBaM) 
library(dplyr) #manipuler un jeu de données
library(psych) #analyses stats, corrélation, tests...

set.seed(2024)
# directory BaM! results, générer toujours les mêmes résultats
dir <- 'C:/Users/famendezrios/Documents/Felipe_MENDEZ/Stagiaires/2024_maree/1_Clara_Chabrillangeas'
setwd(file.path(dir,'Qmec_BaM'))
workspace=file.path(dir,'Qmec_BaM','BaM_workspace')

# Read simulation from original Qmec model 
Qmec_original=read.table('C:/Users/famendezrios/Documents/Felipe_MENDEZ/GitHub/Qmec/results/F2/simulation.txt',header = T)

# Do calibration ?
run_option_calibration = T
# Do prediction ? 
run_option_prediction = T
# dt model
dt_model=60     # in seconds

# test posterior values Qmec (distributions = 'FIX')
test.original.values=T

if(test.original.values==F){
  Q0_only_parameters = F   # Only Q0 as parameter ?
}

# test Qmec prior information  ?
test.prior.qmec.info = F

# Read calibration data (h1,h2, Q)
dat=read.table('calibrationData.txt',header = T, sep = '\t')

D=dataset(X=dat[c('h1','h2')],Y=dat['Q'],Yu=dat['u_Q'],data.dir=workspace)

# Prior information of parameters
if(test.prior.qmec.info==T){
  prior.par.dist=data.frame(Be=c(2000,'LogNormal',7.6,0.1),
                            he=c(15,'LogNormal',2.7,0.15),
                            dzeta=c(-1,'Gaussian',-1,0.15),
                            ne=c(0.05,'LogNormal',-3,0.1),
                            Q0=c(0,'Gaussian',0,15000))
}else{
  prior.par.dist=data.frame(Be=c(2000,'LogNormal',7.6,0.3),
                            he=c(10,'LogNormal',2.5,0.3),
                            dzeta=c(-1,'Gaussian',-1,0.3),
                            ne=c(0.04,'LogNormal',-3,0.15),
                            Q0=c(0,'Gaussian',0,15000))
}


if(test.original.values==F){
  
  if(Q0_only_parameters==T){
    Be=parameter(name='Be',init=as.numeric(prior.par.dist$Be[1]),prior.dist='FIX')
    he=parameter(name='he',init=as.numeric(prior.par.dist$he[1]),prior.dist='FIX')
    dzeta=parameter(name='dzeta',init=as.numeric(prior.par.dist$dzeta[1]),prior.dist='FIX')
    ne=parameter(name='ne',init=as.numeric(prior.par.dist$ne[1]),prior.dist='FIX')
    
    prior.par.dist=data.frame(Q0=prior.par.dist$Q0)
    
    
    Q0=parameter(name='Q0',init=as.numeric(prior.par.dist$Q0[1]),prior.dist=prior.par.dist$Q0[2],prior.par=as.numeric(prior.par.dist$Q0[c(3,4)]))
  }else{
    
    Be=parameter(name='Be',init=as.numeric(prior.par.dist$Be[1]),prior.dist=prior.par.dist$Be[2],prior.par=as.numeric(prior.par.dist$Be[c(3,4)]))
    he=parameter(name='he',init=as.numeric(prior.par.dist$he[1]),prior.dist=prior.par.dist$he[2],prior.par=as.numeric(prior.par.dist$he[c(3,4)]))
    dzeta=parameter(name='dzeta',init=as.numeric(prior.par.dist$dzeta[1]),prior.dist=prior.par.dist$dzeta[2],prior.par=as.numeric(prior.par.dist$dzeta[c(3,4)]))
    ne=parameter(name='ne',init=as.numeric(prior.par.dist$ne[1]),prior.dist=prior.par.dist$ne[2],prior.par=as.numeric(prior.par.dist$ne[c(3,4)]))
    Q0=parameter(name='Q0',init=as.numeric(prior.par.dist$Q0[1]),prior.dist=prior.par.dist$Q0[2],prior.par=as.numeric(prior.par.dist$Q0[c(3,4)]))
  }
  
}else{ # test posterior values following Qmec original results
  
  # Article values
  Be=parameter(name='Be',init=1400,prior.dist='FIX')
  he=parameter(name='he',init=27.07,prior.dist='FIX')
  dzeta=parameter(name='dzeta',init=-0.59,prior.dist='FIX')
  ne=parameter(name='ne',init=0.0483,prior.dist='FIX')
  Q0=parameter(name='Q0',init=0,prior.dist='FIX')
  
  # Posterior values from codeOcean
  Be=parameter(name='Be',init=1406.7822,prior.dist='FIX')
  he=parameter(name='he',init=22.5903 ,prior.dist='FIX')
  dzeta=parameter(name='dzeta',init=-1.1744,prior.dist='FIX')
  ne=parameter(name='ne',init=0.048125,prior.dist='FIX')
  Q0=parameter(name='Q0',init=0,prior.dist='FIX')
  
}

d1=parameter(name='d1',init=-1.379,prior.dist='FIX')
d2=parameter(name='d2',init=-1.958,prior.dist='FIX')
c=parameter(name='c',init=4/3,prior.dist='FIX')
g=parameter(name='g',init=9.81,prior.dist='FIX')
dx=parameter(name='dx',init=38000,prior.dist='FIX')
dt=parameter(name='dt',init=dt_model,prior.dist='FIX')

# Model
M=model(ID='SFDTidal_Qmec',
        nX=2,nY=1, # number of input/output variables
        par=list(Be,he,dzeta,ne,d1,d2,c,g,Q0,dx,dt)) # list of model parameters
# Remnant error
remnant=remnantErrorModel(funk='Constant',
                          par=list(parameter(name='gamma1',init=1000,prior.dist='Uniform',prior.par=c(0,10000))))
# Settings of parameters for mcmc cooking
nCycles = 200 # Number of cycles
burn=0.5      # Percentage of data burned
nSlim=10      # Slim factor: after burning, only one iteration every Nslim is kept.

# Set up configuration files : careful, same quantity of remnants configuration files as number of observations
mcmc_temp=mcmcOptions(nCycles=nCycles)
cook_temp=mcmcCooking(burn=burn,nSlim=nSlim)
mcmcSummary_temp = mcmcSummary(fname = "Config_Summary.txt",
                               result.fname = "Results_Summary.txt",
                               DIC.fname = "Results_DIC.txt",
                               xtendedMCMC.fname = "Results_MCMC_Xtended.txt")
# Run BaM executable
BaM(mod=M,data=D,workspace=workspace,
    doCalib=TRUE,
    doPred=FALSE,
    run=run_option_calibration,
    mcmc=mcmc_temp,
    cook = cook_temp,
    remnant = list(remnant),
    summary = mcmcSummary_temp
)

# Analyse results
# Read 'cooked' MCMC file in the workspace
MCMC=readMCMC(file.path(workspace,'Results_Cooking.txt'))

png(file.path('Graphiques','MCMC_trace.png'),width=2500,height=2000,res=300)
# Trace plot for each parameter, useful to assess convergence.
plots=tracePlot(MCMC)
gridExtra::grid.arrange(grobs=plots,ncol=3)
dev.off()

# Have a look to prior and a posterior density for each parameter
if(test.original.values==F){
  
  if(Q0_only_parameters==T){
    prior.density <- data.frame(Q0=rnorm(nrow(MCMC),mean = as.numeric(prior.par.dist$Q0[3]), sd = as.numeric(prior.par.dist$Q0[4])))
  }else{
    prior.density <- data.frame(Be=rlnorm(nrow(MCMC), meanlog = as.numeric(prior.par.dist$Be[3]), sdlog = as.numeric(prior.par.dist$Be[4])),
                                he=rlnorm(nrow(MCMC),meanlog = as.numeric(prior.par.dist$he[3]), sdlog = as.numeric(prior.par.dist$he[4])),
                                dzeta=rnorm(nrow(MCMC),mean = as.numeric(prior.par.dist$dzeta[3]), sd = as.numeric(prior.par.dist$dzeta[4])),
                                ne=rlnorm(nrow(MCMC),meanlog = as.numeric(prior.par.dist$ne[3]), sdlog = as.numeric(prior.par.dist$ne[4])),
                                Q0=rnorm(nrow(MCMC),mean = as.numeric(prior.par.dist$Q0[3]), sd = as.numeric(prior.par.dist$Q0[4])))
  }
  
  DF <- c()
  for(i in 1:ncol(prior.density)){
    prior.par.DF <- data.frame(value=prior.density[,i],
                               Distributions=rep('Prior',nrow(prior.density)),
                               id=rep(colnames(prior.density)[i]))
    posterior.par.DF <- data.frame(value=MCMC[,i],
                                   Distributions=rep('Posterior',nrow(MCMC)),
                                   id=rep(colnames(MCMC)[i]))
    DF <- rbind(DF,prior.par.DF,posterior.par.DF)
  }
  
  g <- ggplot(DF,aes(x = value ,fill= Distributions))+
    geom_density(alpha=0.4,show.legend = TRUE)+
    facet_wrap(~ id,scales="free",ncol=3)+
    theme_bw()+
    theme(legend.position="top",legend.key.size = unit(0.7, 'cm'))+
    scale_y_continuous(name = "Probability density function")
  
  png(file.path('Graphiques','prior_posterior_density.png'),width = 2300,height = 1600,res=300)
  print(g)
  dev.off()
  
  # Scatterplot
  png(file.path('Graphiques','error_scatterplot_matrix.png'),width = 2200,height = 1600,res=300)
  scatterplot<- pairs.panels(MCMC[,1:ncol(prior.density)],
                             smooth=FALSE,
                             lm=T,
                             scale = T,
                             method = "pearson", # correlation method
                             hist.col = "#00AFBB",
                             cor=T,
                             ci=TRUE,
                             density = TRUE,  # show density plots
                             ellipses = TRUE # show correlation ellipses
  )
  dev.off()
}

# PREDICTIONS
# Define a 'prediction' object for total predictive uncertainty
totalU=prediction(X=dat[c('h1','h2')], # stage values
                  spagFiles='totalU.spag', # file where predictions are saved
                  data.dir=workspace, # a copy of data files will be saved here
                  doParametric=TRUE, # propagate parametric uncertainty, i.e. MCMC samples?
                  doStructural=TRUE) # propagate structural uncertainty ?

# Define a 'prediction' object for parametric uncertainty only - not the doStructural=FALSE
paramU=prediction(X=dat[c('h1','h2')],spagFiles='paramU.spag',data.dir=workspace,
                  doParametric=TRUE,doStructural=FALSE)

# Define a 'prediction' object with no uncertainty - this corresponds to the 'maxpost' discharge maximizing the posterior
maxpost=prediction(X=dat[c('h1','h2')],spagFiles='maxpost.spag',data.dir=workspace,
                   doParametric=FALSE,doStructural=FALSE)
# Define a 'prediction' object with no uncertainty - this corresponds to the discharge computed with prior information
priorU=prediction(X=dat[c('h1','h2')],spagFiles='priorU.spag',data.dir=workspace,
                  priorNsim = nrow(MCMC), doParametric=TRUE,doStructural=FALSE)
# Re-run BaM, but in prediction mode
BaM(mod=M,data=D,remnant=list(remnant), # model and data
    pred=list(totalU,paramU,maxpost,priorU), # list of predictions
    doCalib=FALSE,doPred=TRUE) # Do not re-calibrate but do predictions

# Transform -9999 to NA (missing values)
dat$Q[which(dat$Q==-9999)] <- NA
dat$u_Q[which(dat$u_Q==-9999)] <- NA

# Transform and add date in a different format
dat$date=as.POSIXct(paste(dat$year, dat$month, dat$day, dat$hour, dat$minute, dat$second), format="%Y %m %d %H %M %S")

# Plot spaghetti representing total uncertainty in red
env=read.table(file.path(workspace,'totalU.env'),header=TRUE)
g=ggplot()+geom_ribbon(data=cbind(env,date=dat$date),aes(x=date,ymin=q2.5,ymax=q97.5,fill='Total'),alpha=0.5)

env=read.table(file.path(workspace,'paramU.env'),header=TRUE)
g=g+geom_ribbon(data=cbind(env,date=dat$date),aes(x=date,ymin=q2.5,ymax=q97.5,fill='Parametric'),alpha=0.5)

env=read.table(file.path(workspace,'maxpost.spag'))
g=g+geom_line(data=cbind(env,date=dat$date),aes(x=date,y=V1,color='Simulation_BaM'), linewidth=1)

##### Load Qmec results
# Initial date to capture the tide
date_start <- dat[1,1:6]
# Final date to capture the tide
date_end <- dat[nrow(dat),1:6]

colnames(Qmec_original) <- c(colnames(date_end),'Q')

Qmec_original_gauging <- Qmec_original[dplyr::between(Qmec_original[,1:6],
                                                      date_start[1,],
                                                      date_end[1,]),]

# transform and add date in a different format
Qmec_original_gauging$date=as.POSIXct(paste(Qmec_original_gauging$year, 
                                            Qmec_original_gauging$month, 
                                            Qmec_original_gauging$day, 
                                            Qmec_original_gauging$hour, 
                                            Qmec_original_gauging$minute,
                                            Qmec_original_gauging$second), format="%Y %m %d %H %M %S")

g=g+geom_line(data=Qmec_original_gauging,aes(x=date,y=Q,color='Simulation_Qmec'), linewidth=1)

g=g+geom_ribbon(data=dat,aes(x = date,ymin=Q-u_Q,ymax=Q+u_Q,fill='Gaugings'),alpha=0.15)+
  geom_point(data=dat,aes(date,Q,color='Gaugings'),alpha=0.5)

export=g+labs(x='time (hours)',y='Discharge (m3/s)')+
  scale_color_manual(name=NULL,
                     values=c('Gaugings'='green','Simulation_BaM'='black','Simulation_Qmec'='blue'),
                     labels=c('Gaugings','Simulation BaM!','Simulation Qmec'))+
  scale_fill_manual(name='Uncertainty',
                    values=c('Total'='red','Parametric'='yellow','Gaugings'='green'),
                    labels=c('Gaugings','Parametric','Total'))+
  labs(title = paste0('Calibration between Neuville - Lauzon using Saint-Nicolas gaugings on ',
                      paste(date_start$day,date_start$month,date_start$year,sep='/')))+
  guides(color=guide_legend(override.aes = list(shape=c(16,NA,NA),
                                                linetype =c('blank','solid','solid'))))+
  theme_bw()+
  theme(plot.title = element_text(hjust=0.5,size=13))

png(file.path('Graphiques','simulations_vs_observations.png'),width = 2500,height = 1800,res=300)
print(export)
dev.off()

zoom.gauging=dat$date[!is.na(dat$Q)]

zoom_plot=export+
  coord_cartesian(xlim = c(zoom.gauging[1],rev(zoom.gauging)[1]))


png(file.path('Graphiques','simulations_vs_observations_zoom.png'),
    width = 2500,height = 1800,res=300)
print(zoom_plot)
dev.off()
