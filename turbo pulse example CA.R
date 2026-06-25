# turbo example for luminophore pulse experiments
# author: Clemens Abraham

################################################################################
### Housekeeping ###
################################################################################

# library
library(deSolve)
library(rootSolve)
library(coda)
library(FME)
library(tidyverse)
library(Matrix)
library(plotrix)

# load customized functions
load(file = "turbo pulse functions.RData")

################################################################################
### Example ###
################################################################################

# data
data_all <- read.table("luminophore.txt", header = TRUE)

# profiles that will be modeled
# profilenames <- c("a")               # select profiles manually OR...
profilenames <- unique(data_all$profile) # model every profile

# global parameters
parms <- list(
  db = 0.01,            # adjust if unlikely
  days = 14,
  dx = 0.5,
  cakethickness = 0.5,
  sl = 2,               # adjust if unlikely
  wt =20,               # adjust if unlikely
  k =0,
  flux = 0,
  fluxintroduction = 0
)

# initialize summary table
summary_final <- list()

### automated Modeling ####

# loop over all elements in "profilenames"
for (profilename in profilenames) {
  
  data <- subset(data_all, profile == profilename)
  
  parms$depth <- max(data$end)
  parms$slicenumber <- parms$depth/parms$dx
  parms$times <- c(0,parms$days)
  
  # data profile
  datalimits <- unique(c(data$start,data$end))
  datamidpoints <- midpoints(datalimits)
  conc_profile <- numtoconc(data$lum,datalimits) # if values is an absolute unit
  # conc_profile <- data$lum
  dataprofileframe <- data.frame(depth=datamidpoints, concentration=conc_profile)
  initprof <- initialprofile(conc_profile,datalimits,parms$slicenumber,
                             parms$dx,parms$cakethickness)
  initlimits <- seq(0,parms$slicenumber*parms$dx,by=parms$dx)
  initmidpoints <- midpoints(initlimits)
  initprofileframe <- data.frame(depth=initmidpoints, concentration=initprof)
  
  ### diffusive model ###
  
  # fit diffusive model (ODE)
  fit <- modFit(p=parms$db,f=diffobjective,lower=c(0))
  fitdb <- as.numeric(fit$par[1])
  
  diff_modelprofileframe <- modelprofileframe_diff(fitdb)
  
  ### non local model ###
  
  # fit non local model (CTRW)
  fitnl <- modFit(p = c(sl = parms$sl, wt = parms$wt), f = nlobjective_analytic,
                  lower = c(0,0))
  
  fitsl <- as.numeric(fitnl$par[1])                                               # fitted step length (sl)
  fitwt <- as.numeric(fitnl$par[2])                                               # fitted waiting time (wt)
  fitnldb <- fitsl^2/(2*fitwt)                                                    # calculates db (sl^2/2*wt)
  
  nl_modelprofileframe <- modelprofileframe_nl(sl = fitsl, wt = fitwt)
  
  ### results ###
  
  # plot the results
  par(lwd = 1, cex = 1, cex.lab = 1, cex.axis = 1, mar = c(8, 5, 2, 2))
  revaxis(dataprofileframe$concentration, dataprofileframe$depth,
          xlab = "concentration of luminophores", ylab = "depth in cm",
          pch = 4)
  
  title(main = paste0("profile ", profilename), line = 0)
  
  lines(initprofileframe$concentration, -initprofileframe$depth, lty = 1)
  lines(diff_modelprofileframe$concentration, -diff_modelprofileframe$depth,
        lty = 2)
  lines(nl_modelprofileframe$concentration, -nl_modelprofileframe$depth,
        lty = 3)
  
  legend("right",
        legend = c("initial profile", "observed data",
                  paste0("diffusive (Db = ", round(fitdb, 3), ")"),
                  paste0("non-local (Db = ", round(fitnldb, 3), ")")),
       lty = c(1, NA, 2, 3),
       pch = c(NA, 4, NA, NA),
       bty = "n")
 
  summary_final[[profilename]] <- list(
    fitdb   = round(fitdb, 4),
    fitsl   = round(fitsl, 4),
    fitwt   = round(fitwt, 4),
    fitnldb = round(fitnldb, 4)
  )
}

# summary table
results <- do.call(rbind, lapply(names(summary_final), function(name) {
  data.frame(
    profile = name,
    fitdb   = summary_final[[name]]$fitdb,
    fitsl   = summary_final[[name]]$fitsl,
    fitwt   = summary_final[[name]]$fitwt,
    fitnldb = summary_final[[name]]$fitnldb
  )
}))

print(results)
