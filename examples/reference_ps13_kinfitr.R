options(warn=1)
args <- commandArgs(trailingOnly=TRUE)
.libPaths(c(args[1],.libPaths()))
for (name in c('kinfitr_miscfuncs.R','kinfitr_bloodfuncs.R','kinfitr_weights.R','kinfitr_ma1.R','kinfitr_loganplot.R')) source(file.path(args[2],name))
results <- data.frame()
for (path in list.files(args[3],pattern='_tac.csv$',full.names=TRUE)) {
 scan <- sub('_tac.csv','',basename(path));tac <- read.csv(path);inp0 <- read.csv(file.path(args[3],paste0(scan,'_input.csv')))
 inp <- data.frame(Time=inp0$time,Blood=inp0$blood,Plasma=inp0$plasma,ParentFraction=1,AIF=inp0$plasma);class(inp)<-c('interpblood',class(inp))
 times <- c(0,(tac$start+tac$end)/2);activity <- c(0,tac$activity);duration<-c(0,tac$end-tac$start)
 for (ts in c(20,30,40)) {
  result<-tryCatch({m<-ma1(times,activity,inp,tstar=ts,tstar_type='time',vB=.05,dur=duration);l<-Loganplot(times,activity,inp,tstar=ts,tstar_type='time',vB=.05,dur=duration);data.frame(scan=scan,tstar=ts,ma1_VT=m$par$Vt,logan_VT=l$par$Vt,status='completed')},error=function(e)data.frame(scan=scan,tstar=ts,ma1_VT=NA,logan_VT=NA,status=conditionMessage(e)))
  results<-rbind(results,result);print(result)
 }
}
write.csv(results,file.path(args[3],'kinfitr_graphical.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(args[3],'kinfitr_R_session.txt'))
