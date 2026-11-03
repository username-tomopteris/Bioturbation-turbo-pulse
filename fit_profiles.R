################################################################################
### Luminophoren-Pulse: diffusives und non-lokales Modell an Profile fitten
### Voraussetzung: FME sowie die turbo-Funktionen midpoints(), initialprofile(),
### roughprofile(), diffusive(), nonlocal() und revaxis() sind geladen.
################################################################################

load(file = "turbo.RData")
library(FME)

# ---- Einstellungen -----------------------------------------------------------

datafile      <- "data_PS152_lum.txt"
slicenumber   <- 150          # Modellschichten
dx            <- 0.1          # Schichtdicke (cm)  -> Modelltiefe = slicenumber * dx = 15 cm
cakethickness <- 0.1          # Dicke der Anfangsschicht (cm)
days          <- 14           # Inkubationszeit (Tage)

# Fitgrenzen (anpassen!). Db in cm^2/Jahr, sl in cm.
db_min_year <- 0.01
db_max_year <- 100
sl_min      <- 0.01
sl_max      <- 3            # Obergrenze der Sprungweite; deutlich unter der Profiltiefe wählen
sl_starts   <- c(0.3, 1, 3) # mehrere Startwerte für sl, bester Fit wird behalten

# Plot
xlim_plot <- c(0, 20)         # % luminophores
ytick_step <- 2               # Abstand der y-Achsenbeschriftung (cm)
pdf_file  <- "plots_of_hope_condensed.pdf"
csv_file  <- "results.csv"

model_depth <- slicenumber * dx
db_min <- db_min_year / 365   # intern cm^2/Tag
db_max <- db_max_year / 365

# ---- Hilfsfunktionen (außerhalb der Schleife, alle Daten als Argument) -------

# Kostenfunktion diffusives Modell; x = log(Db)
obj_diff <- function(x, ctx) {
  final <- diffusive(c(0, ctx$days), ctx$initprofile,
                     list(db = exp(x[1]), dx = ctx$dx))
  model <- final[2, 2:(ctx$slicenumber + 1)]
  rough <- roughprofile(model, ctx$datalimits, ctx$dx)
  modCost(obs   = ctx$obs,
          model = data.frame(depth = ctx$datamidpoints, concentration = rough),
          x = "depth")
}

# Kostenfunktion non-lokales Modell; x = (log(Db), log(sl)), wt = sl^2 / (2 Db)
obj_nl <- function(x, ctx) {
  db <- exp(x[1]); sl <- exp(x[2]); wt <- sl^2 / (2 * db)
  final <- nonlocal(c(0, ctx$days), ctx$initprofile,
                    list(waitingtime = wt, steplength = sl, dx = ctx$dx))
  model <- final[2, 2:(ctx$slicenumber + 1)]
  rough <- roughprofile(model, ctx$datalimits, ctx$dx)
  modCost(obs   = ctx$obs,
          model = data.frame(depth = ctx$datamidpoints, concentration = rough),
          x = "depth")
}

# Fit des diffusiven Modells
fit_diffusive <- function(ctx) {
  modFit(p = log(0.01), f = obj_diff, ctx = ctx,
         lower = log(db_min), upper = log(db_max))
}

# Fit des non-lokalen Modells mit mehreren Startwerten
fit_nonlocal <- function(ctx, db_start) {
  best <- NULL
  for (sl0 in sl_starts) {
    f <- try(modFit(p = log(c(db_start, min(sl0, sl_max))), f = obj_nl, ctx = ctx,
                    lower = log(c(db_min, sl_min)), upper = log(c(db_max, sl_max))),
             silent = TRUE)
    if (inherits(f, "try-error")) next
    if (is.null(best) || f$ssr < best$ssr) best <- f
  }
  if (is.null(best)) stop("non-lokaler Fit für alle Startwerte fehlgeschlagen")
  best
}

aic <- function(ssr, n, k) n * log(ssr / n) + 2 * k

# ---- Daten einlesen ----------------------------------------------------------

rawdata      <- read.csv(datafile, header = TRUE)
profilenames <- as.vector(unique(rawdata$profile))

results <- vector("list", length(profilenames))
names(results) <- profilenames

pdf(pdf_file, width = 12, height = 12)
par(mfrow = c(2, 2))

# ---- Schleife über Profile ---------------------------------------------------

for (profilename in profilenames) {

  cat(profilename, "\n")

  res <- tryCatch({

    # alle Tiefen des Profils verwenden (kein Abschneiden bei 10 cm)
    dataframe <- subset(rawdata, profile == profilename)
    dataframe <- dataframe[order(dataframe$depth), ]

    datalimits    <- sort(unique(c(dataframe$start, dataframe$end)))
    datamidpoints <- midpoints(datalimits)
    dataprofile   <- dataframe$end.conc

    # Konsistenzprüfungen
    if (length(datalimits) != nrow(dataframe) + 1)
      warning("Lücken oder Überlappungen zwischen den Messschichten in ", profilename)
    if (max(datalimits) > model_depth)
      stop("Messprofil (", max(datalimits), " cm) tiefer als Modell (", model_depth, " cm)")

    dataprofileframe <- data.frame(depth = datamidpoints, concentration = dataprofile)

    initprofile <- initialprofile(dataprofile, datalimits, slicenumber, dx, cakethickness)
    initlimits    <- seq(0, model_depth, by = dx)
    initmidpoints <- midpoints(initlimits)

    ctx <- list(initprofile = initprofile, datalimits = datalimits,
                datamidpoints = datamidpoints, obs = dataprofileframe,
                slicenumber = slicenumber, dx = dx, days = days)

    # Fits
    fit_d  <- fit_diffusive(ctx)
    fitdb  <- exp(fit_d$par[1])                    # cm^2/Tag
    DB     <- signif(fitdb * 365, 3)               # cm^2/Jahr

    fit_n   <- fit_nonlocal(ctx, db_start = fitdb)
    fitnldb <- exp(fit_n$par[1])                   # cm^2/Tag
    DBNL    <- signif(fitnldb * 365, 3)            # cm^2/Jahr
    fitsl   <- exp(fit_n$par[2])                   # cm
    fitwt   <- fitsl^2 / (2 * fitnldb)             # Tage

    final_diff <- diffusive(c(0, days), initprofile, list(db = fitdb, dx = dx))
    final_nl   <- nonlocal(c(0, days), initprofile,
                           list(waitingtime = fitwt, steplength = fitsl, dx = dx))

    # Gütemaße
    n <- nrow(dataprofileframe)
    at_bound <- function(v, lo, hi) isTRUE(v < lo * 1.01 || v > hi / 1.01)

    list(
      profile      = profilename,
      db_cm2_yr    = fitdb * 365,
      nldb_cm2_yr  = fitnldb * 365,
      sl_cm        = fitsl,
      wt_days      = fitwt,
      ssr_diff     = fit_d$ssr,
      ssr_nl       = fit_n$ssr,
      aic_diff     = aic(fit_d$ssr, n, 1),
      aic_nl       = aic(fit_n$ssr, n, 2),
      db_at_bound  = at_bound(fitdb, db_min, db_max),
      nl_at_bound  = at_bound(fitnldb, db_min, db_max) || at_bound(fitsl, sl_min, sl_max),
      plotdata     = list(dataprofile = dataprofile, datamidpoints = datamidpoints,
                          diff = final_diff[2, 2:(slicenumber + 1)],
                          nl   = final_nl[2, 2:(slicenumber + 1)],
                          depth = initmidpoints,
                          fitdb = fitdb, fitnldb = fitnldb)
    )

  }, error = function(e) {
    message("  Fehler bei ", profilename, ": ", conditionMessage(e))
    NULL
  })

  results[[profilename]] <- res

  # ---- Plot (Modelllinien bis 15 cm) -----------------------------------------
  if (!is.null(res)) {
    pd <- res$plotdata
    par(lwd = 1, cex = 1, cex.lab = 1, cex.axis = 1, mar = c(8, 5, 4, 2))

    # Tiefe nach unten (ylim umgekehrt), y-Achse manuell in Zweierschritten
    plot(pd$dataprofile, pd$datamidpoints,
         xlab = "% luminophores", ylab = "depth in cm",
         pch = 4, xlim = xlim_plot, ylim = c(model_depth, 0),
         yaxt = "n", bty = "o")
    axis(2, at = seq(0, model_depth, by = ytick_step), las = 1)
    title(main = profilename, line = 2)

    lines(pd$diff, pd$depth, lty = 1)
    lines(pd$nl,   pd$depth, lty = 2)

    legend("bottomright",
           legend = as.expression(c(
             "observed data",
             bquote("diffusive (" * D[b] * " = " * .(DB) * " cm"^2 * "/year)"),
             bquote("CTRW ("      * D[b] * " = " * .(DBNL) * " cm"^2 * "/year)")
           )),
           lty = c(NA, 1, 2), pch = c(4, NA, NA), bty = "n")
  }
}

dev.off()

# ---- Ergebnistabelle ---------------------------------------------------------

rows <- lapply(results, function(r) {
  if (is.null(r)) return(NULL)
  as.data.frame(r[setdiff(names(r), "plotdata")], stringsAsFactors = FALSE)
})
allresults <- do.call(rbind, rows)

print(allresults)
write.csv(allresults, file = csv_file, row.names = FALSE)

# Hinweis: Profile mit nl_at_bound = TRUE haben einen Wert an der Fitgrenze.
# Db ist dann nur als Unter-/Obergrenze zu interpretieren.
