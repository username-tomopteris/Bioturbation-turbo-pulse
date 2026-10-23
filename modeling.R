# modeling PS152 MUC LUM bioturbation data
# author: Clemens Abraham

################################################################################
### Housekeeping ###
################################################################################

# library
library(dplyr)
library(deSolve)
library(rootSolve)
library(coda)
library(FME)
library(tidyverse)
library(Matrix)
library(plotrix)
library(readxl)
library(outliers)
library(readODS)

# load customized functions
load(file = "turbo pulse functions.RData")

################################################################################
### data and settings ###
################################################################################

# metadata
metadata <- read_ods("metadata PS 152.ods")

# data
data_all <- read_ods("data PS152 MUC LUM.ods")
data_all$profile <- paste(data_all$station, data_all$`liner id`, data_all$treatment)
data_all$dataname <- paste(data_all$profile, data_all$depth)

# basic calculations
data_all$`RFU-blank`[data_all$`RFU-blank` < 0] <- 0
data_all$`layer RFU` <- data_all$`RFU-blank`*data_all$`thickness coefficient`
data_all$`mean layer RFU` <- ave(data_all$`layer RFU`, data_all$dataname, FUN = mean)

# model data
data_all %>%
  group_by(dataname) %>%
  summarise(profile = unique(profile), station = unique(station),
            treatment = unique(treatment),
            start = unique(start), end = unique(end), depth = unique(depth),
            thickness_coefficient = unique(`thickness coefficient`),
            RFU_mean = mean(`layer RFU`)) -> data_model

# station and treatment vector
data_all %>%
  group_by(profile) %>%
  summarise(profile = unique(profile), station = unique(station),
            treatment = unique(treatment)) -> core_list

# calculating tracer sum per core
data_model %>%
  group_by(profile) %>%
  summarise(sum_lum = sum(RFU_mean)) -> sum_lum

data_model <- left_join(data_model, sum_lum, by = "profile")

data_model$init.conc <- ifelse(data_model$start == 0.0, 100, 0) # initial profile (%)

data_model$end.conc <- data_model$RFU_mean*100/data_model$sum_lum # end concentration (%)

### settings ####

# profiles that will be modeled

# profilenames <- c("KN_10 AWI_21 control", "KN_10 AWI_17 natural",
#                   "KN_15 UHRO_111 control", "KN_15 AWI_13 natural")  # select profiles manually OR...
profilenames <- unique(data_model$profile) # model every profile

# global parameters
parms <- list(
  db = 0.005,           # adjust if unlikely
  days = 14,
  dx = 0.1,             # model slice thickness (more slices = more precision?)
  cakethickness = 0.1,  # >= dx
  sl = 0.5,             # adjust if unlikely
  wt =2,                # adjust if unlikely
  k =0,
  flux = 0,
  fluxintroduction = 0
)

# initialize summary table
summary_final <- list()

### automated modeling ####

# pdf("plots_of_hope.pdf", width = 6, height = 6)
# par(mfrow = c(1, 1))

# loop over all elements in "profilenames"
for (profilename in profilenames) {
  
  data <- subset(data_model, profile == profilename)
  
  parms$depth <- max(data$end)
  parms$slicenumber <- parms$depth/parms$dx
  parms$times <- c(0,parms$days)
  
  # data profile
  datalimits <- unique(c(data$start, data$end))
  datamidpoints <- midpoints(datalimits)
  conc_profile <- numtoconc(data$end.conc,datalimits, convert = FALSE) # no conversion
  dataprofileframe <- data.frame(depth = datamidpoints,
                                 concentration = conc_profile)
  initprof <- initialprofile(conc_profile, datalimits, parms$slicenumber,
                             parms$dx, parms$cakethickness)
  
  initlimits <- seq(0,parms$slicenumber*parms$dx, by=parms$dx)
  initmidpoints <- midpoints(initlimits)
  initprofileframe <- data.frame(depth=initmidpoints, init.conc=initprof)
  
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
  # par(lwd = 1, cex = 1, cex.lab = 1, cex.axis = 1, mar = c(8, 5, 2, 2))
  # 
  # revaxis(c(dataprofileframe$concentration, NA), c(dataprofileframe$depth, 15.1),
  #         yrev = TRUE, xlab = "% luminophores", ylab = "depth in cm",
  #         pch = 4, xlim = c(0, 10), ylim = c(-8, 0), bty = "o")
  # 
  # title(main = profilename, line = 0)
  # 
  # # lines(initprofileframe$init.conc, -initprofileframe$depth, lty = 1)
  # lines(diff_modelprofileframe[[2]]$profile, -diff_modelprofileframe[[2]]$depth,
  #       lty = 1)
  # lines(nl_modelprofileframe[[2]]$profile, -nl_modelprofileframe[[2]]$depth,
  #       lty = 2)
  # 
  # legend("right",
  #       legend = c("observed data",
  #                 paste0("diffusive (Db = ", round(fitdb, 4), ")"),
  #                 paste0("non-local (Db = ", round(fitnldb, 4), ")")),
  #      lty = c(NA, 1, 2),
  #      pch = c(4, NA, NA),
  #      bty = "n")
 
  summary_final[[profilename]] <- list(
    fitdb   = round(fitdb, 4),
    fitsl   = round(fitsl, 4),
    fitwt   = round(fitwt, 4),
    fitnldb = round(fitnldb, 4),
    station = unique(data$station),
    treatment = unique(data$treatment)
  )
}

# dev.off()

# summary table
results <- do.call(rbind, lapply(names(summary_final), function(name) {
  data.frame(
    profile = name,
    station = summary_final[[name]]$station,
    treatment = summary_final[[name]]$treatment,
    fitdb   = summary_final[[name]]$fitdb,
    fitsl   = summary_final[[name]]$fitsl,
    fitwt   = summary_final[[name]]$fitwt,
    fitnldb = summary_final[[name]]$fitnldb
  )
}))

results <- unique.data.frame(results)
results <- left_join(results, metadata)

print(results)
write_ods(results, "results.ods")

### current tasks ####

# box plots (Db~station...)
par(mar = c(5, 5, 2, 2), mfrow = c(1, 2))

boxplot(data = results, fitdb ~ depth, subset = treatment == "natural",
        xlab = "depth in m", ylab = "Db", 
        main = "natural", ylim = c(0, 0.01))
boxplot(data = results, fitnldb ~ depth, subset = treatment == "natural",
        xlab = "depth in m", ylab = "non-local Db",
        main = "natural", ylim = c(0, 0.4))

boxplot(data = results, fitdb ~ depth, subset = treatment == "control",
        xlab = "depth in m", ylab = "Db",
        main = "control", ylim = c(0, 0.01))
boxplot(data = results, fitnldb ~ depth, subset = treatment == "control",
        xlab = "depth in m", ylab = "non-local Db",
        main = "control", ylim = c(0, 0.1))

# alternative 1
par(mfrow = c(1, 2))

boxplot(data = results, fitdb ~ treatment:depth, col = c("grey", "lightblue"),
        xlab = "depth in m", ylab = "Db", main = "Biodiffusion",
        ylim = c(0, 0.01), xaxt = "n") +
  axis(1, at = seq(1.5, 7.5, 2), labels = c("287", "923", "1355", "2159"))

boxplot(data = results, fitnldb ~ treatment:depth, col = c("grey", "lightblue"),
        xlab = "depth in m", ylab = "non-local Db", main = "CTRW",
        ylim = c(0, 0.4), xaxt = "n")
  legend("bottomright", c("control", "natural"), fill = c("grey", "lightblue")) +
  axis(1, at = seq(1.5, 7.5, 2), labels = c("287", "923", "1355", "2159"))
  

# alternative 2
par(mfrow = c(1, 1))

plot(data = results, fitdb ~ depth, subset = treatment == "natural",
     xlab = "depth in m", ylab = "") +
  points(data = results, fitnldb ~ depth)



### results plot ####
core_list <- left_join(core_list, results, by = "profile")
core_list <- filter(core_list, treatment.x == "natural")
core_list <- data.frame(profile = core_list$profile,
                        station = core_list$station.x,
                        fitdb = core_list$fitdb,
                        fitsl = core_list$fitsl,
                        fitwt = core_list$fitwt,
                        fitnldb = core_list$fitnldb)
core_list <- unique.data.frame(core_list)

dixon.test(core_list$fitdb)
dixon.test(core_list$fitnldb)
dixon.test(core_list$fitsl)
dixon.test(core_list$fitwt)

# remove outliers
core_list$fitnldb[9] <- NA
core_list$fitsl[9] <- NA

par(mfrow = c(2,1))

barplot(core_list$fitdb, names.arg = core_list$station) +
  title(main = "Db (diffusion)")
barplot(core_list$fitnldb, names.arg = core_list$station) +
  title(main = "Db (CTRW)")
barplot(core_list$fitsl, names.arg = core_list$station) +
  title(main = "step length (CTRW)")
barplot(core_list$fitwt, names.arg = core_list$station) +
  title(main = "waiting time (CTRW)")

### scatter plot ####

# pdf("RFU profiles PS152.pdf", width = 8, height = 6)

for (profilename in profilenames) {
  sub_data <- subset(data_all, profile == profilename)
  par(lwd = 1, cex = 1, cex.lab = 1, cex.axis = 1, mar = c(8, 5, 2, 2))
  revaxis(sub_data$RFU, sub_data$depth,
          xlab = "RFU", ylab = "depth in cm",
          pch = 4, xlim = c(0, 1e+07))
  
  title(main = profilename, line = 0)
}

dev.off()



### box plot without top layer ####

for (profilename in profilenames) {
  
  sub_data <- subset(data_all, profile == profilename)
  sub_data <- sub_data[-c(1:3), ]
  
  par(lwd = 1, cex = 1, cex.lab = 1, cex.axis = 1, mar = c(5, 5, 4, 2))
  
  # Achsengrenzen bestimmen
  xlim <- range(0, max(sub_data$RFU))
  ylim <- rev(range(0, max(sub_data$depth)))        # umgekehrt für Tiefenachse
  
  # leeren Plot mit umgekehrter y-Achse initialisieren
  plot(NULL, xlim = xlim, ylim = ylim,
       xlab = "RFU", ylab = "depth in cm")
  
  # Boxplots je Tiefe einzeichnen
  depths <- sort(unique(sub_data$depth))
  for (d in depths) {
    vals <- sub_data$RFU[sub_data$depth == d]
    boxplot(vals, at = d, add = TRUE,
            horizontal = TRUE,                             # Boxen horizontal (x = RFU)
            boxwex = 0.3,                                  # Breite der Box in Tiefeneinheiten
            axes = FALSE)
  }
  
  title(main = paste0("profile ", profilename), line = 1)
}

