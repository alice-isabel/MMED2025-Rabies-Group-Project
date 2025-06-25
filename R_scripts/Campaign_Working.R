library(tidyverse)
library(ggplot2)
library(tidyr)
library(dplyr)

event_sir <- function(time, S, E, I ,V, params, t_end) { #JDC: count.inf removed from args
  
  with(as.list(params), {
    N <- S+I+E+V
    
    print(time)
    print(floor(time))
    print(floor(time)%%365))
  
    if ((floor(time)+365)&&365 < 31){
    rates <- c(
      infect = Lambda*S*I/N,
      infectious = Gamma*E,
      vaccinated = Tau*S/31,
      waned = Psi*V,
      birth = B*N,
      diseasedDeath = Mu2*I,
      backDeathS = Mu1 * S,
      backDeathV = Mu1 * V,
      backDeathE = Mu1 * E,
      backDeathI = Mu1 * I
      
    )} else {
      rates <- c(
        infect = Lambda*S*I/N,
        infectious = Gamma*E,
        vaccinated = 0,
        waned = Psi*V,
        birth = B*N,
        diseasedDeath = Mu2*I,
        backDeathS = Mu1 * S,
        backDeathV = Mu1 * V,
        backDeathE = Mu1 * E,
        backDeathI = Mu1 * I
        )}
    # print(rates)
    
    total_rate <- sum(rates)
    # print(total_rate)
    print(total_rate)
    print(total_rate==0)
    if (total_rate < 10^(-6)) {
      
      count.I <- 0
      count.S <- 0
      count.E <- 0 
      count.V <-0
      event_time <- t_end
      
    } else {
      
      event_time <- time + rexp(1, total_rate) 
      
      event_type <- sample(c("infect",
                             "infectious",
                             "vaccinated" ,
                             "waned",
                             "birth",
                             "diseasedDeath",
                             "backDeathS",
                             "backDeathV",
                             "backDeathE",
                             "backDeathI"
      ), 1, prob = rates / total_rate)
      
      switch(event_type,
             "infect" = {
               S <- S - 1
               E <- E + 1
               
               count.I <- 0
               count.S <- 0
               count.E <- 1 
               count.V <-0
             },
             "infectious" = {
               E <- E - 1
               I <- I + 1
               count.I <- 1
               count.S <- 0
               count.E <- 0 
               count.V <-0
             },
             "waned" = {
               V <- V - 1
               S <- S + 1
               count.I <- 0
               count.S <- 1
               count.E <- 0 
               count.V <-0
             },
             "birth" = {
               S <- S + 1
               count.I <- 0
               count.S <- 1
               count.E <- 0 
               count.V <-0
             },
             "diseasedDeath" = {
               I <- I - 1
               count.I <- 0
               count.S <- 0
               count.E <- 0 
               count.V <-0
             },
             "backDeathS" = {
               S <- S - 1
               count.I <- 0
               count.S <- 0
               count.E <- 0 
               count.V <-0
             },
             "backDeathV" = {
               V <- V - 1
               count.I <- 0
               count.S <- 0
               count.E <- 0 
               count.V <-0
             },
             "backDeathE" = {
               E <- E - 1
               count.I <- 0
               count.S <- 0
               count.E <- 0 
               count.V <-0
             },
             "backDeathI" = {
               I <- I - 1
               count.I <- 0
               count.S <- 0
               count.E <- 0 
               count.V <-0
             },
             "vaccinated" = {
               S <- S - 1
               V <- V + 1
               count.I <- 0
               count.S <- 0
               count.E <- 0 
               count.V <-0
             }
      )
    }
    
    return(data.frame(time = event_time, S = S, I = I, E = E,V=V, 
                      count.I = count.I,
                      count.S = count.S,
                      count.E = count.E,
                      count.V = count.V))
  })
}

## Function to simulate states from time 0 to t_end:

simulate_sir <- function(t_end, y, params) {
  # with(as.list(y), {
  events <- vector("list", 100000000)
  idx<- 1
  count.I2 <- 0
  count.S2 <- 0
  count.E2 <- 0 
  count.V2 <-0
  
  # initial state
  S <- y["S"]
  E<- y["E"]
  I<- y["I"]
  V<- y["V"]
  
  time<-0
  # store inital row
  events[[idx]]<- c(time, S,E,I,V,0,0,0,0)
  idx<- idx+1
  
  # This will be the cumulative count
  
  #count.int initialisation changed. At present, it doesn't include the infection events that seed the outbreak...
  # ts <- data.frame(time = 0, S = S, I = I, E=E,V = V, count.I = 0,
  #                  count.S = 0,
  #                  count.E = 0,
  #                  count.V = 0) 
  # next_event <- ts
  # print(ts)
  while (time < t_end && idx<=length(events)) {
    next_event <- event_sir(time, S, E, I,V, params, t_end)#JDC: count.inf removed from args
    time <- next_event$time
    S<- next_event$S
    E<- next_event$E
    I<- next_event$I
    V<- next_event$V
    count.I2 <- count.I2 + next_event$count.I
    count.S2 <- count.S2 + next_event$count.S
    count.E2 <- count.E2 + next_event$count.E
    count.V2 <- count.V2 + next_event$count.V
    
    next_event$count.I <- count.I2 # replace with the cumulative version
    next_event$count.S <- count.S2
    next_event$count.E <- count.E2
    next_event$count.V <- count.V2
    
    events[[idx]]<- unlist(next_event)
    idx<- idx+1
    
    
    # print(next_event)
    print(next_event$time)
    # ts <- rbind(ts, next_event)
  }
  #print(idx)
  out<- do.call(rbind, events[1:(idx-1)])
  colnames(out)<- c("time", "S", "E", "I", "V", "count.I", "count.S", "count.E", "count.V")
  return(as.data.frame(out))
  
}
## Run the model for specified inputs:

pop <- 100      
# population size
params <- c(
  B = 0.003, # unknown
  Mu1=1/(25*30), 
  Mu2=1/(3.1),
  Psi=1/(2.5*365),
  Tau=70/100, #time independent vaccination rate
  Lambda=2, # unknown
  Gamma=1/22.3, 
  N=pop,
  Tau2=5/100 #off-campaign tau
)    	# parameter values
final_time <- 365*1                      	# end time
y0 <- c(
  S = (1-0.02)*pop,
  E = (0.01)*pop,
  I = (0.01)*pop, 
  V=0
)        	# initial state

ts1 <- simulate_sir(final_time, y0, params)

## And plot:

ggplot(ts1, aes(x = time, y = count.I)) +
  geom_step(color = "purple", linewidth = 1.2) +
  labs(title = "Cumulative New Infections Over Time",
       y = "Cumulative Infections",
       x = "Time (days)") +
  theme_minimal(base_size = 14)
view(ts1)

ggplot(ts1, aes(x = time, y = S)) +
  geom_step(color = "darkgreen", linewidth = 1.2) +
  labs(title = "Number of Susceptible Dogs Over Time",
       y = "Susceptible Count",
       x = "Time (days)") +
  theme_minimal(base_size = 14)


ggplot(ts1, aes(x = time, y = I)) +
  geom_step(color = "red", linewidth = 1.2) +
  labs(title = "Number of Infectious Dogs Over Time",
       y = "Infectious Count",
       x = "Time (days)") +
  theme_minimal(base_size = 14)

ggplot(ts1, aes(x = time, y = V)) +
  geom_step(color = "blue", linewidth = 1.2) +
  labs(title = "Number of Vaccinated Dogs Over Time",
       y = "Vaccinated Count",
       x = "Time (days)") +
  theme_minimal(base_size = 14)

ts1 <- ts1 %>%
  mutate(Total = S + E + I + V)

ggplot(ts1, aes(x = time, y = Total)) +
  geom_step(color = "black", linewidth = 1.2) +
  labs(title = "Total Population Over Time",
       y = "Population Size",
       x = "Time (days)") +
  theme_minimal(base_size = 14)
