##--function to compute e-folding timescale---###
efold <- function(alpha){
  t    <- seq(0,60,0.001)
  cors <- alpha^t
  e    <- (cors - (1/exp(1)))^2
  return(t[e==min(e)][1])
}

f_plot_efolds <- function(D,xlim,breaks,ylim){
  for(j in 1:length(regions)){
    d     <- D[[1]] %>% filter(region==regions[j])#,complete.cases(chl,kd_490,sst,par,lat,lon,depth,month,daylength))
    E     <- array(NA,dim=c(dim(d)[1],dim(d)[2],30))
    for(i in 1:30) E[,,i] <- data.matrix(D[[i]] %>% filter(region==regions[j]))
    colnames(E) <- colnames(d)
    
    dd <- t(E[,colnames(d)==vars[1],30:1]) %>% .[,apply(.,2,function(x) sum(is.finite(x))>5)]
    cors <- apply(dd,2,function(x){
      n <- length(x)
      return(cor(x[2:n],x[1:(n-1)],use='complete.obs'))
    })
    
    es <- unlist(lapply(cors,function(x) efold(x)))
    
    if(j==1) hist(es,main='',breaks=breaks,xlim=xlim,freq=FALSE,ylim=ylim,col=cols[1])
    if(j!=1) hist(es,main='',breaks=breaks,xlim=xlim,freq=FALSE,ylim=ylim,col=cols[1],yaxt='n')
    mtext(nms[j],adj=0)
    abline(v=mean(es,na.rm=TRUE),lty=2,col=cols[1])
    if(j==1) legend('top',bty='n',legend=c('Chl','SST','PAR'),col=cols[c(3,2,1)],pch=15,cex=1.5)
    
    for(p in 3:2){
      dd <- t(E[,colnames(d)==vars[p],30:1]) %>% .[,apply(.,2,function(x) sum(is.finite(x))>5)]
      cors <- apply(dd,2,function(x){
        n <- length(x)
        return(cor(x[2:n],x[1:(n-1)],use='complete.obs'))
      })
      es <- unlist(lapply(cors,function(x) efold(x)))
      
      hist(es,breaks=breaks,add=TRUE,col=adjustcolor(cols[p],alpha.f=0.7),freq=FALSE)
      if(vars[p]!="sst") abline(v=mean(es,na.rm=TRUE),col=cols[p],lty=2)
    }
  }
  mtext(outer=TRUE,side=1,expression(italic(e)*'-Folding Timescale [days]'),line=2.5)
  mtext(outer=TRUE,side=2,'Frequency Density',line=0.5)
}
