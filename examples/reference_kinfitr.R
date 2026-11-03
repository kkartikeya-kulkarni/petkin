args <- commandArgs(trailingOnly=TRUE)
# args: R library directory, fetched R source directory, PBR28 RData, output directory
.libPaths(c(args[1],.libPaths()))
dir.create(args[4],recursive=TRUE,showWarnings=FALSE)
for (name in c('kinfitr_miscfuncs.R','kinfitr_bloodfuncs.R','kinfitr_weights.R','kinfitr_1tcm.R','kinfitr_2tcm.R','kinfitr_ma1.R','kinfitr_loganplot.R')) source(file.path(args[2],name))
load(args[3])
results <- data.frame()
for (i in seq_len(nrow(pbr28))) {
  tac <- as.data.frame(pbr28$tacs[[i]])
  actual <- tac[tac$Duration>0,]
  start <- actual$StartTime[1]/60
  times <- c(0, if(start>0) start/2, actual$Times/60)
  duration <- c(0, if(start>0) start, actual$Duration/60)
  activity <- c(0, if(start>0) 0, actual$FC)
  inp <- pbr28$input[[i]]
  m <- ma1(times,activity,inp,tstar=40,tstar_type='time',vB=.05,dur=duration)
  l <- Loganplot(times,activity,inp,tstar=40,tstar_type='time',vB=.05,dur=duration)
  result <- data.frame(scan=as.character(pbr28$PET[i]),ma1_VT=m$par$Vt,logan_VT=l$par$Vt)
  results <- rbind(results,result)
  print(result)
}
write.csv(results,file.path(args[4],'kinfitr_graphical.csv'),row.names=FALSE)
# Unmodified kinfitr 2TCM on three illustrative cases; midpoint observation
# differs from frame averaging, so these are diagnostic comparisons only.
comp <- data.frame()
for (i in c(1,3,11)) {
  tac <- as.data.frame(pbr28$tacs[[i]])
  result <- tryCatch(twotcm(tac$Times/60,tac$FC,pbr28$input[[i]],weights=tac$Weights,
                           vB=.05,inpshift=0,K1.lower=1e-6,k2.lower=1e-6,k3.lower=1e-6,k4.lower=1e-6,
                           K1.upper=3,k2.upper=3,k3.upper=3,k4.upper=3),error=function(e)e)
  if(inherits(result,'error')) row <- data.frame(scan=as.character(pbr28$PET[i]),K1=NA,k2=NA,k3=NA,k4=NA,VT=NA,status=conditionMessage(result)) else
    row <- data.frame(scan=as.character(pbr28$PET[i]),K1=result$par$K1,k2=result$par$k2,k3=result$par$k3,k4=result$par$k4,VT=result$par$Vt,status='completed')
  comp <- rbind(comp,row);print(row)
}
write.csv(comp,file.path(args[4],'kinfitr_compartment.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(args[4],'kinfitr_R_session.txt'))
