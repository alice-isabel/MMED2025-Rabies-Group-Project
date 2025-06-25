library(tidyverse)
library(ggplot2)
library(tidyr)
library(dplyr)


?sample

event_sir <- function(time, S, E, I ,V, params, t_end) { #JDC: count.inf removed from args
  
  with(as.list(params), {
    N <- S+I+E+V
    
    
    rates <- c(
      infect = Lambda*S*I/N,
      infectious = Gamma*E,
      vaccinated = Tau*S,
      waned = Psi*V,
      birth = B*N,
      diseasedDeath = Mu2*I,
      backDeathS = Mu1 * S,
      backDeathV = Mu1 * V,
      backDeathE = Mu1 * E,
      backDeathI = Mu1 * I
      
    )
    # print(rates)
    
    total_rate <- sum(rates)
    # print(total_rate)
    if (total_rate == 0) {
      
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
    
    return(data.frame(time = event_time, S = S, I = I, E = E,V=V, count.I = count.I,
                      count.S = count.S,
                      count.E = count.E,
                      count.V = count.V))
  })
}

## Function to simulate states from time 0 to t_end:











simulate_sir <- function(t_end, y, params) {
  with(as.list(y), {
    count.I2 <- 0
    count.S2 <- 0
    count.E2 <- 0 
    count.V2 <-0
    
    # This will be the cumulative count
    
    #count.int initialisation changed. At present, it doesn't include the infection events that seed the outbreak...
    ts <- data.frame(time = 0, S = S, I = I, E=E,V = V, count.I = 0,
                     count.S = 0,
                     count.E = 0,
                     count.V = 0) 
    next_event <- ts
    print(ts)
    while (next_event$time < t_end) {
      next_event <- event_sir(next_event$time, S = next_event$S, I = next_event$I, V = next_event$V, E = next_event$E , params, t_end)#JDC: count.inf removed from args
      count.I2 <- count.I2 + next_event$count.I
      count.S2 <- count.S2 + next_event$count.S
      count.E2 <- count.E2 + next_event$count.E
      count.V2 <- count.V2 + next_event$count.V
      
      
      next_event$count.I <- count.I2 # replace with the cumulative version
      next_event$count.S <- count.S2
      next_event$count.E <- count.E2
      next_event$count.V <- count.V2
      # print(next_event)
      print(next_event$time)
      ts <- rbind(ts, next_event)
    }
    print(ts)
    
    return(ts)
  })
}
## Run the model for specified inputs:

pop <- 2300000      
# population size
params <- c(
   B = 0.003, # unknown
   Mu1=1/(25*30), 
   Mu2=1/(3.1),
   Psi=1/(2.5*365),
   Tau=70/(100*365),
   Lambda=0.49, # unknown
   Gamma=1/22.3, 
   N=pop
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

ggplot(ts1, aes(x = time, y = I)) +
  geom_step(color = "red", linewidth = 1.2) +
  labs(title = "Number of Infectious Individuals Over Time",
       y = "Infectious Count",
       x = "Time (days)") +
  theme_minimal(base_size = 14)

ggplot(ts1, aes(x = time, y = V)) +
  geom_step(color = "blue", linewidth = 1.2) +
  labs(title = "Number of Vaccinated Individuals Over Time",
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







# Generate the data and produce the plot multiple times, by running the relevant
# lines above. What do you see?

## PART 2: SIR with spillover
## --------------------------------------------------------------------

# We now introduce the possibility of infections occurring in the
# population due to other exposures - for example, due to infected animals
# in some other maintenance population briefly entering the territory. 

# TASK 1: Make a copy of all the code presented for part 1 and modify it to 
# introduce spillover events - ADDCODE
# make functions event_sirspill and simulate_sirspill.
# More specifically, assume that, in addition to the transmission already 
# occurring, there is an additional rate of infection of susceptibles of 
# lambda*S/N (total rate). Here lambda is the rate at which, for example,
# animals from outside of the population manage to make contact 
# with animals in our population of interest. 

# When you are ready, scroll down to the bottom of the tutorial to find an
# example solution, and then move onto the next parts of the tutorial using it.

# Here are some hints if needed: 
#   Update the list of '## Transitions'
#   In the function event_sirspill, include a third event type 'Spillover'
#   Also try to include a counter for spillovers (count.spillovers) - 
#     this may prove handy if you get time to explore. 
#   To run the model, think about whether any of the inputs needs to change.


## PART 3: Explore patterns
## --------------------------------------------------------------------

# Let's compare outputs for different parameter values. .

# To make this easier, let's write a function that plots I over time,
#   showing several different possible trajectories, for a set of input 
#   parameter values

plot_sirspill <- function(nTraj = 16, params, final_time, y0){
  
  ts_collect <- list()
  for (ii in 1:nTraj){
    ts_collect[[ii]] <- simulate_sirspill(final_time, y0, params)
  }
  
  ts_collect_df <- (bind_rows(ts_collect, .id = "Simulation")
                    |> mutate(Simulation = factor(Simulation, levels = 1:nTraj))
  )
  
  gg <- ggplot(ts_collect_df, aes(x = time, y = I)) +
    geom_step(linewidth = 1.2) +
    facet_wrap(~Simulation) +
    labs(title = "I over time", y = "Count", x = "Time") +
    scale_x_continuous(breaks = NULL) +
    coord_cartesian(xlim=c(0, final_time)) +
    theme_minimal(base_size = 14)

  print(gg)
}

# We can now easily now easily plot trajectories - here we plot 
# 16 trajectories for a set of inputs

plot_sirspill(nTraj=16
	, params=c(beta = 0.3, gamma = 0.1,lambda = 0.01)
	, final_time=400
	, y0 = c(S = pop - 1, I = 1, R = 0)
)

# TASK 2: 
# Use the function plot_sirspill to compare outputs for different parameter
# values (specified below) and think about the patterns you see. Also think 
# about how the  patterns you see relate to R0 = beta/gamma, and, later in
# time, Reff = R0*S/N
# Pay attention to the scale of y-axis.

# Don't forget to erase your graphs every now and then using graphics.off().

# (1) 1 infected animal to begin, no spillover, with different beta and gamma values
# c(beta = 0.3, gamma = 0.1, lambda = 0.00), 400, c(S = pop - 1, I = 1, R = 0))
# versus 
# c(beta = 0.2, gamma = 0.2, lambda = 0.00), 400, c(S = pop - 1, I = 1, R = 0))

# (2) 0 infected animals to begin, with spillover, with different beta and gamma values 
# c(beta = 0.3, gamma = 0.1, lambda = 0.01), 400, c(S = pop, I = 0, R = 0))
# versus 
# c(beta = 0.2, gamma = 0.2, lambda = 0.01), 400, c(S = pop, I = 0, R = 0))

# (3) as for (2) but with a higher spillover rate
# c(beta = 0.3, gamma = 0.1, lambda = 0.04), 400, c(S = pop, I = 0, R = 0))
# versus 
# c(beta = 0.2, gamma = 0.2, lambda = 0.04), 400, c(S = pop, I = 0, R = 0))

# Lastly, would we have been able to study these patterns using 
# a deterministic model? If you have time, refer back to previous labs, 
# and try fit and plot a corresponding deterministic model. 

## PART 2: Solution 
## --------------------------------------------------------------------

## Our version of the functions for part 2 can be found in
## ICI3D_Ex1_StochasticSpillover_functions.R in the tutorials repo: 
## https://raw.githubusercontent.com/ICI3D/RTutorials/refs/heads/master/ICI3D_Ex1_StochasticSpillover_functions.R

# Compartments:
# (S,I,R) = (susceptible, infectious, removed)

# Transitions:
# Event                           Change        		 Rate
# Spillover (S)			  (S,I,R)->(S-1,I+1,R)		 lambda*S/N
# Infection (S)                   (S,I,R)->(S-1,I+1,R)           beta*I*S/N
# Recovery/Removal (I)            (S,I,R)->(S,I-1,R+1)           gamma*I


## Run the model for specified inputs:

params <- c(beta = 0.3, gamma = 0.1
            , lambda = 0.02)            	# parameter values
final_time <- 400                       	# end time
y0 <- c(S = pop - 1, I = 1, R = 0)        	# initial state

ts1 <- simulate_sirspill(final_time, y0, params)

## And plot:

ts1_long <- (ts1 
             |> pivot_longer(cols = c(S, I, R), names_to = "Compartment", values_to = "count")
             |> mutate(Compartment = factor(Compartment, levels = c('S','I','R')))
)

ggplot(ts1_long, aes(x = time, y = count, color = Compartment)) +
  geom_step(linewidth = 1.2) +
  labs(title = "SIR dynamics with spillover", y = "Count", x = "Time") +
  theme_minimal(base_size = 14)

