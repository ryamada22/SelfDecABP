% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/manual R codes.r
\name{color_scale}
\alias{color_scale}
\title{Rescale the Colors Based on Raw Dataset.}
\usage{
color_scale(obs, neg.cols, pos.cols)
}
\arguments{
\item{obs}{a verctor of raw data}

\item{neg.cols}{the colors to be assigned for negatives.}

\item{pos.cols}{the colors to be assigned for positives.}
}
\value{
the re-scaled colors corresponding the raw data points.
}
\description{
Color-scale based on raw data points with your visualized colors.
}
\details{
rgb() is used for rescaling the colors and "gdata" package is necessary to be installed.
}
\examples{
load("mean.success.per.n.exp.100.201.RData")
load("mean.success.per.n.wpbeta.100.201.RData")  
N <- 100
ps<-combn(seq(0,1,by=0.01),2)
w<-seq(from=-1,to=1,by=0.01)
Difs<-array(0,c(length(w),length(ps[1,]),N))
   for(i in 1:N){
     Difs[,,i]<-t(t(mean.success.per.n.wpbeta[,,i+1])-mean.success.per.n.exp[i+1,])
   }
   w.id<-c(21,61,141,181)
   #the positition of w values in the vector w, such as w=0.8,-0.4,0.4,or 0.8
   N.id<-c(10,30,100)
   Difs.<-Difs[w.id,,N.id]
   cols<-color_scale(obs=c(Difs.),neg.cols=c("blue","white"),pos.cols=c("white","red"))
   par(mfrow=c(1,1),mar=c(3,3,3,3))
   image(matrix(seq(0,1,length.out=length(c(Difs.))),ncol=1),col=cols[order(c(Difs.))], xlab=" ", ylab=" ",xaxt="n",yaxt="n")
   axis(3,at=c(0.3),labels=c("Negatives"),cex.axis=1.4,font=4,las=1,lwd=4,tick=FALSE,line=-0.2)
   axis(3,at=c(0.85),labels=c("Positives"),cex.axis=1.4,font=4,las=1,lwd=4,tick=FALSE,line=-0.2)
   se<-range(c(Difs.))
   axis(3,at=c(0,0.5,1),labels=round(c(min(Difs.),median(Difs.),max(Difs.)),3),cex.axis=1.1,font=4,las=1,line=-0.2,lwd=4,col.ticks="green",col="green")
   colbar<-seq(0,1,length.out=length(cols))
   axis(3,at=colbar[length(which(Difs.<0))],labels=c("0"),cex.axis=1.4,font=4,las=1,lwd=4,line=-0.2,col.ticks="green")
}
\keyword{NA}

