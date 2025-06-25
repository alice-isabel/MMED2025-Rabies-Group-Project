library(dplyr)
library(ggplot2)
library(tidyr)

rm(list=ls())
getwd()
source("./R_scripts/Stochastic_Model.R")

run_multiple_sims <- function(n_sim, t_end, y, params, time_step = 1) {

  
  
  all_sims <- vector("list", n_sim)
  
  for (i in 1:n_sim) {
    sim_result <- simulate_sir(t_end, y, params)
    
  
    sim_result$time <- floor(sim_result$time / time_step) * time_step
    sim_result$sim <- i
    
    all_sims[[i]] <- sim_result
  }
  
  combined <- bind_rows(all_sims)
  
  
  summary_stats <- combined %>%
    group_by(time) %>%
    summarise(
      S_mean = mean(S),
      E_mean = mean(E),
      I_mean = mean(I),
      V_mean = mean(V)
      
    )
  
  
  return(summary_stats)
}

n_sim <- 15
pop <- 100      
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
)
final_time <- 365*1                      	
y0 <- c(
  S = (1-0.02)*pop,
  E = (0.01)*pop,
  I = (0.01)*pop, 
  V=0
)        

average_results <- run_multiple_sims(
  n_sim = n_sim,
  t_end = final_time,
  y = y0,
  params = params,
  time_step = 1
)

average_results

average_results <- average_results %>%
  pivot_longer(cols = c(S_mean, E_mean, I_mean, V_mean),
               names_to = "Compartment", values_to = "count")
ggplot(average_results, aes(x = time, y = count, color = Compartment)) +
  geom_line(linewidth = 1) +
  theme_minimal() +
  labs(title = paste("Average of", n_sim, "Simulations"),
       x = "Time", y = "Average Count")



tau_values <- c(15/(100*365), 40/(100*365), 70/(100*365), 100/(100*365))


all_tau_results <- list()

for (tau in tau_values) {
  

  params_current <- params
  params_current["Tau"] <- tau
  

  avg_res <- run_multiple_sims(
    n_sim = n_sim,
    t_end = final_time,
    y = y0,
    params = params_current,
    time_step = 1
  )
  
  avg_res <- avg_res %>%
    mutate(Tau = tau) %>%
    pivot_longer(cols = c(S_mean, E_mean, I_mean, V_mean),
                 names_to = "Compartment", values_to = "count")
  
  all_tau_results[[as.character(tau)]] <- avg_res
}


combined_tau_results <- bind_rows(all_tau_results)
combined_tau_results

ggplot(combined_tau_results, aes(x = time, y = count, color = Compartment)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ Tau, scales = "free_y") +
  theme_minimal() +
  labs(title = paste("Average of", n_sim, "Simulations for Different Tau Values"),
       x = "Time", y = "Average Count")

















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
  
  #print(gg)
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

