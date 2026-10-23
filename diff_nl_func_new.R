# diff
modelprofileframe_diff <- function(db) {
  dx <- parms$dx
  slicenumber <- parms$slicenumber
  final             <- diffusive(times = parms$times, initprof = initprof,
                                 parameters = list(db = db, dx = dx))
  modelprofile <- final[2,2:(slicenumber+1)]
  roughmodelprofile <- roughprofile(modelprofile,datalimits,dx)
  modelprofileframe <- list(
    rough <- data.frame(
               depth=datamidpoints,
               concentration=roughmodelprofile),
    fine <-   data.frame(
               profile=modelprofile,
               depth=(seq_len(parms$slicenumber) - 0.5) * parms$dx
               )
    )
  
  return(modelprofileframe)
}

# nl
modelprofileframe_nl <- function(sl, wt) {
  
  final             <- nonlocal_analytic(parameters = list(sl = sl, wt = wt))
  modelprofileframe <- list(
    rough <- data.frame(
              depth = datamidpoints,
              concentration = roughprofile(final, datalimits, parms$dx)),
    fine <- data.frame(
             profile=as.numeric(final),
             depth=(seq_len(parms$slicenumber) - 0.5) * parms$dx
             )
  )
    
  return(modelprofileframe)
}
