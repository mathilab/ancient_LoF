## Created: 06/16/21 by Nina Tansey
## Description: density and number line plots for oe for LoF variants 

#Preamble
#check your working directory and set if if appropriate
getwd()
setwd()
#setwd('/Users/ntansey/Desktop/PENNworks/Research')
library("xlsx")
library("plotrix")
#

allLoFs<-read.csv("outputs/sgdp_allLoFs.csv")
hetLoFs<-read.csv("outputs/sgdp_hetLoFs.csv")
homLoFs<-read.csv("outputs/sgdp_homLoFs.csv")
transitionLoFs<-read.csv("outputs/sgdp_transitionLoFs.csv")
transversionLoFs<-read.csv("outputs/sgdp_transversionLoFs.csv")

#number line plot (hom, het) 
pdf("figs/oe_sgdp_homhet.pdf")
xmin = min(homLoFs$oe_lof,hetLoFs$oe_lof)
xmax = max(homLoFs$oe_lof,hetLoFs$oe_lof)
plot(0:1, axes=FALSE,type="n",xlab="oe_lof",ylab="", main="SGDP oe_lof Visualization", xlim=c(xmin,xmax))
axis(1,pos=0) 
points(hetLoFs$oe_lof, (hetLoFs$oe_lof-hetLoFs$oe_lof), col="gold")
points(homLoFs$oe_lof, (homLoFs$oe_lof-homLoFs$oe_lof), col="firebrick")
legend("topleft",
       c("homLoFs","hetLoFs"),
       fill=c("firebrick", "gold"))
dev.off()


#density plot (all, hom, het)
pdf("figs/doe_sgdp_homhet.pdf")
alld<-density(allLoFs$oe_lof)
hetd<-density(hetLoFs$oe_lof)
homd<-density(homLoFs$oe_lof)
xmin = min(hetd$x,homd$x,alld$x)
xmax = max(hetd$x,homd$x,alld$x)
plot(hetd, col="gold", main="SGDP Density Plot", xlim=c(xmin,xmax))
lines(homd$x, homd$y, col="firebrick")
lines(alld$x, alld$y, col="cornflowerblue")
legend("topleft",
       c("allLoFs","homLoFs","hetLoFs"),
       fill=c("cornflowerblue","firebrick", "gold"))
dev.off()

#number line plot (transition, transversion) 
pdf("figs/oe_sgdp_tt.pdf")
xmin = min(transitionLoFs$oe_lof,transversionLoFs$oe_lof)
xmax = max(transitionLoFs$oe_lof,transversionLoFs$oe_lof)
plot(0:1, axes=FALSE,type="n",xlab="oe_lof",ylab="", main="SGDP oe_lof Visualization", xlim=c(xmin,xmax))
axis(1,pos=0) 
points(transitionLoFs$oe_lof, (transitionLoFs$oe_lof-transitionLoFs$oe_lof), col="darkorchid")
points(transversionLoFs$oe_lof, (transversionLoFs$oe_lof-transversionLoFs$oe_lof), col="forestgreen")
legend("topleft",
       c("transitionLoFs","transversionLoFs"),
       fill=c("darkorchid", "forestgreen"))
dev.off()


#density plot (all, transition, transversion)
pdf("figs/doe_sgdp_tt.pdf")
alld<-density(allLoFs$oe_lof)
transitiond<-density(transitionLoFs$oe_lof)
transversiond<-density(transversionLoFs$oe_lof)
xmin = min(transitiond$x,transversiond$x,alld$x)
xmax = max(transitiond$x,transversiond$x,alld$x)
ymax = max(transitiond$y,transversiond$y,alld$y)
plot(transitiond, col="darkorchid", main="SGDP Density Plot", xlim=c(xmin,xmax), ylim=c(0,ymax))
lines(transversiond$x, transversiond$y, col="forestgreen")
lines(alld$x, alld$y, col="cornflowerblue")
legend("topleft",
       c("allLoFs","transitionLoFs","transversionLoFs"),
       fill=c("cornflowerblue","darkorchid", "forestgreen"))
dev.off()


#density plot (all, hom, het, transition, transversion)
pdf("figs/doe_sgdp_all.pdf")
alld<-density(allLoFs$oe_lof)
hetd<-density(hetLoFs$oe_lof)
homd<-density(homLoFs$oe_lof)
transitiond<-density(transitionLoFs$oe_lof)
transversiond<-density(transversionLoFs$oe_lof)
xmin = min(hetd$x,homd$x,alld$x,transitiond$x,transversiond$x)
xmax = max(hetd$x,homd$x,alld$x,transitiond$x,transversiond$x)
ymax = max(hetd$y,homd$y,alld$y,transitiond$y,transversiond$y)
plot(hetd, col="gold", main="SGDP Density Plot", xlim=c(xmin,xmax), ylim=c(0,ymax))
lines(homd$x, homd$y, col="firebrick")
lines(alld$x, alld$y, col="cornflowerblue")
lines(transitiond$x, transitiond$y, col="darkorchid")
lines(transversiond$x, transversiond$y, col="forestgreen")
legend("topleft",
       c("allLoFs","homLoFs","hetLoFs", "transitionLoFs", "transversionLoFs"),
       fill=c("cornflowerblue","firebrick", "gold", "darkorchid", "forestgreen"))
dev.off()
