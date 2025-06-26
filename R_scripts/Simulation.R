library(dplyr)
library(ggplot2)
library(tidyr)

rm(list=ls())
getwd()
source("./R_scripts/Stochastic_Model.R")

run_multiple_sims <- function(n_sim, t_end, y, params, time_step = 1,local_sim=simulate_sir) {
  
  
  
  all_sims <- vector("list", n_sim)
  
  for (i in 1:n_sim) {
    sim_result <- local_sim(t_end, y, params)
    
    
    sim_result$time <- floor(sim_result$time / time_step) * time_step
    sim_result$sim <- i
    
    all_sims[[i]] <- sim_result
  }
  
  combined <- bind_rows(all_sims)
  print("ddddddddddddddddddddddddddddddddddddd")
  
  summary_stats <- combined %>%
    group_by(time) %>%
    summarise(
      S_mean = mean(S),
      E_mean = mean(E),
      I_mean = mean(I),
      V_mean = mean(V),
      Total_mean = mean(S)+mean(E)+mean(I)+mean(V)
      
    )
  
  
  return(summary_stats)
}

