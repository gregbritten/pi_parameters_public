library(viridis)
library(randomForest)
source('r/main_process_nc.R')

FITS <- readRDS("results/FITS.rds")
RSQ  <- readRDS("results/RSQ.rds")

day_Ek    <- 25
day_PBmax <- 30

#############################################
## RSQ ######################################
#############################################
pdf('plots/RSQ.pdf',height=4,width=6.5)
par(mfrow=c(1,2),mar=c(2,1,2,1.5),oma=c(3,3,3,3),cex.axis=0.8)
p1 <- unlist(RSQ[["reg"]][["PBmax"]])
p2 <- unlist(RSQ[["lm"]][["PBmax"]])
  plot(p1,type='l',ylim=c(0.3,0.8),bty='n')
  lines(p2,lty=2)
  abline(v=day_PBmax,lty=1)
  abline(v=which(p2==max(p2)),lty=2)
  mtext(side=2,expression(italic('R'^2)),line=2.5,cex=1.2)
  mtext(expression(italic('a) P'['max']^'B')~'[mg C (mg chla)'^{-1}~'h'^{-1}*']'),adj=0)

p1 <- unlist(RSQ[["reg"]][["Ek"]])
p2 <- unlist(RSQ[["lm"]][["Ek"]])
  plot(p1,type='l',ylim=c(0.3,0.8),yaxt='n',cex.axis=0.8,bty='n'); axis(side=2,labels=NA)
  lines(p2,lty=2)
  abline(v=day_Ek,lty=1)
  abline(v=which(p2==max(p2)),lty=2)
  mtext(side=1,outer=TRUE,'Environmental Averaging Timescale [Days]',line=0.5)
  mtext(expression(italic('b) E'['k'])~'['*mu*'mol quanta m'^{-2}~'s'^{-1}*']'),adj=0)
dev.off()

#############################################
## EFFECTS ##################################
#############################################
#day_Ek    <- 10
#day_PBmax <- 10


cols   <- viridis::turbo(5)
factor <- 1E6/(86400)

d_PBmax <- D_nc[[day_PBmax]] %>% filter(is.finite(par) & is.finite(sst) & is.finite(pico) & is.finite(depth) & depth < mld)
d_Ek    <- D_nc[[day_Ek]]    %>% filter(is.finite(par) & is.finite(sst) & is.finite(pico) & is.finite(depth) & depth < mld)

#pdf('plots/partial_effects_10day.pdf',height=4,width=8)
pdf('plots/partial_effects.pdf',height=4,width=8)
par(mfrow=c(2,4),mar=c(0.5,0.5,0.5,0.5),oma=c(4,5,3,3),cex.lab=1,cex.axis=1)

##--PBmax--###########################
##par
plot(-999,xlim=c(0,80)*factor,ylim=c(2,5),xlab='',ylab='',xaxt='n'); axis(side=1,labels=FALSE); axis(side=2,labels=FALSE)
coefs <- c()
for(i in 1:length(regions)){
  xx     <- partialPlot(FITS[["reg"]][['PBmax']][[day_PBmax]],x.var='par',pred.data= d_PBmax %>% filter(region==regions[i]) ,plot=FALSE)
  lines(xx$x*factor,xx$y,col=cols[i])
  fit <- lm(xx$y ~ I(scale(xx$x*factor)))
  coefs <- c(coefs,summary(fit)$coefficients[2,1])
}
legend('bottomright',cex=0.7,legend=round(coefs,4),bty='n',col=cols,lty=1)
mtext(side=2,expression(italic('P'['max']^'B')~'[mg C (mg chla)'^{-1}~'h'^{-1}*']'),line=2.5,cex=0.75)
mtext('a)',adj=0.05,line=-1.5)
legend('topright',legend=region_long,lty=1,col=cols,bty='n',cex=0.8)

##sst
plot(-999,xlim=c(-2,30),ylim=c(2,5),xlab='',ylab='',yaxt='n',xaxt='n'); axis(side=2,labels=FALSE); axis(side=1,labels=FALSE)
coefs <- c()
for(i in 1:length(regions)){
  xx <- partialPlot(FITS[["reg"]][["PBmax"]][[day_PBmax]],x.var='sst',pred.data=d_PBmax %>% filter(region==regions[i]),plot=FALSE)
  lines(xx$x,xx$y,col=cols[i])
  fit <- lm(xx$y ~ I(scale(xx$x*factor)))
  coefs <- c(coefs,summary(fit)$coefficients[2,1])
}
legend('bottomright',cex=0.7,legend=round(coefs,4),bty='n',col=cols,lty=1)
mtext('b)',adj=0.05,line=-1.5)

##pico
plot(-999,xlim=c(0,0.8),ylim=c(2,5),xlab='',ylab='',yaxt='n',xaxt='n'); axis(side=2,labels=FALSE); axis(side=1,labels=FALSE)
coefs <- c()
for(i in 1:length(regions)){
  xx <- partialPlot(FITS[["reg"]][["PBmax"]][[day_PBmax]],x.var='pico',pred.data=d_PBmax %>% filter(region==regions[i]),plot=FALSE)
  lines(xx$x,xx$y,col=cols[i])
  fit <- lm(xx$y ~ I(scale(xx$x*factor)))
  coefs <- c(coefs,summary(fit)$coefficients[2,1])
}
legend('bottomright',cex=0.7,legend=round(coefs,4),bty='n',col=cols,lty=1)
mtext('c)',adj=0.05,line=-1.5)

##depth
plot(-999,xlim=c(0,50),ylim=c(2,5),xlab='',ylab='',yaxt='n',xaxt='n'); axis(side=2,labels=FALSE); axis(side=1,labels=FALSE)
coefs <- c()
for(i in 1:length(regions)){
  xx <- partialPlot(FITS[["reg"]][["PBmax"]][[day_PBmax]],x.var='depth',pred.data=d_PBmax %>% filter(region==regions[i]),plot=FALSE)
  lines(xx$x,xx$y,col=cols[i])
  fit <- lm(xx$y ~ I(scale(xx$x*factor)))
  coefs <- c(coefs,summary(fit)$coefficients[2,1])
}
legend('topright',cex=0.7,legend=round(coefs,4),bty='n',col=cols,lty=1)
mtext('d)',adj=0.05,line=-1.5)

##--Ek--######################
##par
plot(-999,xlim=c(0,80)*factor,ylim=c(0,200),xlab='',ylab='',xaxt='n'); axis(side=1)
coefs <- c()
for(i in 1:length(regions)){
  xx <- partialPlot(FITS[["reg"]][["Ek"]][[day_Ek]],x.var='par',pred.data=d_Ek %>% filter(region==regions[i]),plot=FALSE)
  lines(xx$x*factor,xx$y,col=cols[i])
  fit <- lm(xx$y ~ I(scale(xx$x*factor)))
  coefs <- c(coefs,summary(fit)$coefficients[2,1])
}
legend('bottomright',cex=0.7,legend=round(coefs,3),bty='n',col=cols,lty=1)
mtext(side=2,expression(italic('E'['k'])~'['*mu*'E m'^{-2}~'s'^{-1}*']'),line=2.5,cex=0.75)
mtext(side=1,expression(italic('PAR')~'['*mu*'E m'^{-2}~'s'^{-1}*']'),line=2.75)
mtext('e)',adj=0.05,line=-1.5)

##sst
plot(-999,xlim=c(-2,30),ylim=c(0,200),xlab='',ylab='',xaxt='n',yaxt='n'); axis(side=1); axis(side=2,labels=FALSE)
coefs <- c()
for(i in 1:length(regions)){
  xx <- partialPlot(FITS[["reg"]][["Ek"]][[day_Ek]],x.var='sst',pred.data=d_Ek %>% filter(region==regions[i]),plot=FALSE)
  lines(xx$x,xx$y,col=cols[i])
  fit <- lm(xx$y ~ I(scale(xx$x*factor)))
  coefs <- c(coefs,summary(fit)$coefficients[2,1])
}
legend('bottomright',cex=0.7,legend=round(coefs,3),bty='n',col=cols,lty=1)
mtext(side=1,expression(italic('SST')~'['*degree*'C]'),line=2.75)
mtext('f)',adj=0.05,line=-1.5)

##pico
plot(-999,xlim=c(0,0.8),ylim=c(0,200),xlab='',ylab='',xaxt='n',yaxt='n'); axis(side=1); axis(side=2,labels=FALSE)
coefs <- c()
for(i in 1:length(regions)){
  xx <- partialPlot(FITS[["reg"]][["Ek"]][[day_Ek]],x.var='pico',pred.data=d_Ek %>% filter(region==regions[i]),plot=FALSE)
  lines(xx$x,xx$y,col=cols[i])
  fit <- lm(xx$y ~ I(scale(xx$x*factor)))
  coefs <- c(coefs,summary(fit)$coefficients[2,1])
}
legend('bottomright',cex=0.7,legend=round(coefs,3),bty='n',col=cols,lty=1)
mtext(side=1,expression(italic('c'['pico']*'/c'['total'])),line=2.75)
mtext('g)',adj=0.05,line=-1.5)

##depth
plot(-999,xlim=c(0,50),ylim=c(0,200),xlab='',ylab='',xaxt='n',yaxt='n'); axis(side=1); axis(side=2,labels=FALSE)
coefs <- c()
for(i in 1:length(regions)){
  xx <- partialPlot(FITS[["reg"]][["Ek"]][[day_Ek]],x.var='depth',pred.data=d_Ek %>% filter(region==regions[i]),plot=FALSE)
  lines(xx$x,xx$y,col=cols[i])
  fit <- lm(xx$y ~ I(scale(xx$x*factor)))
  coefs <- c(coefs,summary(fit)$coefficients[2,1])
}
legend('right',cex=0.7,legend=round(coefs,3),bty='n',col=cols,lty=1)
mtext(side=1,expression(italic('Depth')~'[m]'),line=2.75)
mtext('h)',adj=0.05,line=-1.5)
dev.off()


####################################################
## MIXED EFFECTS REGRESSION EFFECTS ################
####################################################
plotlm <- function(fit,add=FALSE,ylims,labs){
  coef <- summary(fit)$coefficients[c(1,3,2,5,4),1]  #reordering variables to be consistent across plots
  ses  <- summary(fit)$coefficients[c(1,3,2,5,4),2]
  n    <- length(coef)
  if(add==FALSE){
    #plot(1:n,coef,ylim=c(min(coef-2*ses),max(coef+2*ses)),pch=19,xaxt='n')
    plot(1:n,coef,pch=19,xaxt='n',ylim=ylims,xlim=c(1,n+0.5),xaxt='n',bty='n',yaxt='n',cex=0.8)
    segments(x0=1:n,x1=1:n,y0=coef-2*ses,y1=coef+2*ses)
  }
  if(add==TRUE){
    points(1:n+0.25,coef,ylim=c(min(coef-2*ses),max(coef+2*ses)),pch=19,col='red')
    segments(x0=1:n+0.25,x1=1:n+0.25,y0=coef-2*ses,y1=coef+2*ses,col='red')
  }
  abline(h=0,lty=2)
  #axis(side=1,at=1:n,labels=row.names(summary(fit)$coefficients),las=2)
  axis(side=1,at=1:n,labels=labs,las=2)
}

pdf('plots/lm_effects.pdf',height=3.75,width=7)
par(mfrow=c(1,2),mar=c(2,2,2,0),oma=c(4,4,2,2),cex.axis=0.8)
plotlm(FITS[["lm"]][["PBmax"]][[34]],ylim=c(-0.8,0.8),labs=c(expression(italic("Intercept")),
                                                             expression(italic("PAR")),
                                                             expression(italic("SST")),
                                                             expression(italic('c'['pico']*'/c'['total'])),
                                                             expression(italic("Depth"))))
axis(side=2,at=seq(-0.8,0.8,0.2))
mtext(expression(italic('a) P'['max']^'B')~'[mg C (mg chla)'^{-1}~'h'^{-1}*']'),adj=0)

mtext(side=2,expression('Standardized Slope ['*Delta*'sdY/'*Delta*'sdX]'),line=2.5)
plotlm(FITS[["lm"]][["Ek"]][[24]],ylim=c(-0.8,0.8),labs=c(expression(italic("Intercept")),
                                                   expression(italic("PAR")),
                                                   expression(italic("SST")),
                                                   expression(italic('c'['pico']*'/c'['total'])),
                                                   expression(italic("Depth"))))
mtext(expression(italic('b) E'['k'])~'['*mu*'mol quanta m'^{-2}~'s'^{-1}*']'),adj=0)
dev.off()


######################################################
## RAW PBmax vs. Ek ##################################
######################################################
pdf('plots/PBmax_vs_Ek.pdf',height=5,width=6)
rsqs <- numeric(5)
par(mfrow=c(1,1))
plot(-999,xlim=c(-0,15),ylim=c(0,500),xlab='',ylab='')
for(i in 1:length(regions)){
  dd <- d %>% filter(region==regions[i])
  points(dd$PBmax,dd$Ek,col=adjustcolor(cols[i],alpha.f=0.4),pch=19,lwd=0)
  fit <- lm(Ek ~ PBmax,data=dd)
  if(regions[i]%in%c('spac','scot')){
    PBmax_in <- seq(min(dd$PBmax),max(dd$PBmax)-2,length.out=1000)
  }else{PBmax_in <- seq(min(dd$PBmax),max(dd$PBmax+3),length.out=1000)}
  preds <- predict(fit,newdata=list(PBmax=PBmax_in))
  lines(PBmax_in,preds,col=cols[i],lwd=3)
  rsqs[i] <- round(summary(fit)$r.squared,3)
}
mtext(expression(italic('P'['max']^'B')~'[mg C (mg chla)'^{-1}~'h'^{-1}*']'),side=1,line=2.5)
mtext(expression(italic('E'['k'])~'['*mu*'mol quanta m'^{-2}~'s'^{-1}*']'),side=2,line=2.5)

legend('topright',legend=c(bquote(.(region_long[1])~'R'^2~'='~.(rsqs[1])),
                           bquote(.(region_long[2])~'R'^2~'='~.(rsqs[2])),
                           bquote(.(region_long[3])~'R'^2~'='~.(rsqs[3])),
                           bquote(.(region_long[4])~'R'^2~'='~.(rsqs[4])),
                           bquote(.(region_long[5])~'R'^2~'='~.(rsqs[5]))),cex=0.7,bty='n',col=adjustcolor(cols,alpha.f=0.5),pch=19)
dev.off()

