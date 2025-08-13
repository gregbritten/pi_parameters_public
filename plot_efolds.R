source('r/f_plot_efolds.R')

lets <- c('a)','b)','c)','d)','e)','f)')
nms <- numeric(6)
for(i in 1:6) nms[i] <- paste(lets[i],region_long[i])
cols <- c(turbo(4)[c(3,4)],'dark green')

pdf('plots/efolds.pdf',height=3,width=11)
par(mfrow=c(1,5),mar=c(0,2,0,0),oma=c(4,5,3,3),cex.axis=0.8)
f_plot_efolds(D_nc,breaks=seq(0,60,4),ylim=c(0,0.2),xlim=c(0,60))
dev.off()

