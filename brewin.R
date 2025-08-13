
##--parameters from Brewin paper--############
G <- c(-1.51, -1.25, 14.95, 0.25)
H <- c( 0.29,  3.05, 16.24, 0.56)
J <- c( 3.70,  1.13, 14.89, 0.569)
K <- c( 0.503, 1.33, 17.31, 0.258)

##--intermediate functions--########
Cm_pn <- function(G, SST) return(1 - (G[1]/(1+exp(-G[2]*(SST-G[3]))) + G[4]))
Cm_p  <- function(H, SST) return(1 - (H[1]/(1+exp(-H[2]*(SST-H[3]))) + H[4]))
D_pn  <- function(J, SST) return(J[1]/(1+exp(-J[2]*(SST-J[3]))) + J[4])
D_p   <- function(K, SST) return(K[1]/(1+exp(-K[2]*(SST-K[3]))) + K[4])

##--functions to compute fractions and ratios--############
brewin <- function(G, H, J, K, C, SST){
  C_pn <- Cm_pn(G,SST)*(1-exp(-(D_pn(J,SST)/Cm_pn(G,SST))*C))
  C_p  <- Cm_p(H,SST)*(1-exp(-(D_p(K,SST)/Cm_p(H,SST))*C))
  C_n  <- C_pn - C_p
  C_m  <- C - C_pn
  return(c(C_p,C_n,C_m))
}

pico <- function(G,H,J,K,C,SST){
  bb <- brewin(G,H,J,K,C,SST)
  return(max(bb[1],0,na.rm=TRUE)/C)
}

nano_pico <- function(G,H,J,K,C,SST){
  bb <- brewin(G,H,J,K,C,SST)
  bb[2]/bb[1]
}

micro_nano <- function(G,H,J,K,C,SST){
  bb <- brewin(G,H,J,K,C,SST)
  bb[3]/bb[2]
}

micro_pico <- function(G,H,J,K,C,SST){
  bb <- brewin(G,H,J,K,C,SST)
  max(bb[3],0,na.rm=TRUE)/max(bb[1],0,na.rm=TRUE)
}








