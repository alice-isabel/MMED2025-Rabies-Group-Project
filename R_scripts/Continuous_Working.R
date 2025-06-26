library(dplyr)
library(ggplot2)
library(tidyr)

rm(list=ls())
getwd()
source("./R_scripts/Simulation.R")

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
  pivot_longer(cols = c(S_mean, E_mean, I_mean, V_mean,Total_mean),
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
    pivot_longer(cols = c(S_mean, E_mean, I_mean, V_mean,Total_mean),
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

