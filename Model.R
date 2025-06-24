rm(list=ls())

sir <- function(t,y,parms){

  with(c(as.list(y),parms),{
    dSdT <- -Lambda*S*I/N -Tau*S +Psi*V+B*N- Mu1 * S
    dEdT <- Lambda*S*I/N - Gamma*E -Mu1 * E
    dIdT <- Gamma*E-Mu1*I-Mu2*I
    dVdT <- Tau*S - Psi * V - Mu1 * V
    
    return(list(c(dSdT,dEdT,dIdT,dVdT)))
  })
}




time <- 0



N0 <- 1000000


pop.SI <- c(S = (1-0.02)*N0,
            E = (0.01)*N0,
            I = (0.01)*N0, 
            V=0 )       



values <- c(
  B = 0.0005,
  Mu1=0.005,
  Mu2=0.01,
  Psi=0.05,
  Tau=0.1,
  Lambda=0.49, # be
  Gamma=0.5, 
  N=N0)



sir(t=time,y=pop.SI,parms=values)


delta.t <- 0.1                    

pop.next <- pop.SI + unlist(sir(t=time,y=pop.SI,parms=values)) * delta.t
pop.SI
pop.next



library(deSolve)               
time.out <- seq(0,365,0.1)  

lsoda(
  y = pop.SI,               
  times = time.out,             
  func = sir,                   
  parms = values                
)

ts.sir <- data.frame(lsoda(
  y = pop.SI,               
  times = time.out,         
  func = sir,               
  parms = values            
))

head(ts.sir)

subset(ts.sir,time==365)

plot(ts.sir$time,               
     ts.sir$I,                  
     xlab = "Time in days",     
     ylab = "Number infected",  
     main = "Rabies",           
     xlim = c(0,400),           
     type = "l",                
     bty = "n")                 

