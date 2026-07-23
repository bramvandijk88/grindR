# Load the package (installed). During development use instead:
#   devtools::load_all("~/grindR")
library(grindR)

model <- function(t, state, parms) {
  with(as.list(c(state,parms)), {
    dR <- b*R*(1 - R/K) - d*R - a*R*N/(R+h)
    dN <- c*a*R*N/(R+h) - delta*N
    return(list(c(dR, dN)))  
  }) 
}  
p <- c(b=1,d=0.1,K=2,a=1,c=1,delta=0.3,c=1,h=1)  #parameters
s <- c(R=1.1,N=0.1) #begin-dichtheden van R en N  
par(mfrow=c(1,2)) # optioneel, maakt 2 plotjes naast elkaar

plane(xmax=1.8,ymax=1.5,vector=1,grid=3)
run(tstep=0.5, legend=F, draw=points)
