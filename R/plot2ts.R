#' @title Plot two time series with different Y scales
#' 
#' @description Plot two time series with different Y scales
#' @param ts1 time series 1 (ts1); object of class "ts"
#' @param ts2 time series 2 (ts2); object of class "ts"
#' @param ylim_ts1 the y limits (x1, x2) for ts1. The default value is NULL, indicating that the ylim are set to the min and max value of ts1
#' @param ylim_ts2 the y limits (x1, x2) for ts2. The default value is NULL, indicating that the ylim are set to the min and max value of ts2
#' @param lab_ts1 a title for ts1
#' @param lab_ts2 a title for ts2
#' @param col_ts1 colour for ts1. The default value is "black"
#' @param col_ts2 colour for ts2. The default value is "blue" 
#' @param points add points at time series observations. The default value is TRUE.
#' @param title plot title
#' 
#' @author Johannes Reiche
#' 
#' @export 

plot2ts <- function(ts1,ts2,ylim_ts1=NULL,ylim_ts2=NULL,lab_ts1="",lab_ts2="",col_ts1="black",col_ts2="blue", points=TRUE, title="") {
  
  #omit NAs from irregular TS
  zts1 <- na.omit(as.zoo(ts1))
  zts2 <- na.omit(as.zoo(ts2))

  #set par for the plot
  par(mar=c(5, 4, 4, 6))
  
  #set y data range for ts1
  #a - take it from data range (in case of an undefined range)
  if (is.null(ylim_ts1)){
    ylim_ts1 = c(range(zts1))
  #b - use defined range  
  } else {
    ylim_ts1 <- ylim_ts1
  }
  #set y data range for ts2
  #a - take it from data range (in case of an undefined range)
  if (is.null(ylim_ts2)){
    ylim_ts2 = c(range(zts2))
  #b - use defined range  
  } else {
    ylim_ts2 <- ylim_ts2
  }
  
  #plot ts1
  p <- plot(zts1,xlim=range(time(ts1)),axes=F, ylim=ylim_ts1, xlab="", ylab="",type="l",col=col_ts1,pch=".")
  #add points
  if(points==TRUE) {p <- p + points(zts1,xlim=range(time(ts1)),ylim=ylim_ts1, pch = 19,cex=0.5)}
  p <- axis(2, ylim=ylim_ts1,col=col_ts1,lwd=1)
  #add ylab for ts1
  mtext(2,text=lab_ts1,line=2)
  #add time lables to xaxis
  axis(1,pretty(range(time(ts1)),10))
  #add title
  mtext(3,text=paste(title,sep=""),line=1)
  
  #plot ts2
  par(new=T)
  p <- plot(zts2,xlim=range(time(ts1)),ylim=ylim_ts2, axes=F, xlab="", ylab="",type="l",col=col_ts2,pch=".")
  #add points for ts2
  if(points==TRUE) {p <- p + points(zts2,xlim=range(time(ts2)),ylim=ylim_ts2, pch = 19,cex=0.5,col="blue")}
  p <- axis(4, ylim=ylim_ts2,col=col_ts2,col.axis=col_ts2,lwd=1)
  #add ylab for ts2
  mtext(4,text=lab_ts2,line=2,col=col_ts2)
}


