



prefix_in_path 	<- ""
prefix_out_path <- ""

# Load functions for subsequent scripts
source(paste0(prefix_in_path, 'functions_simulations.R'))


# model parameters to define for generations simulation

# pop size (estimated effective size of the whole pop)
N = 1000
# sample size
sample_size <- 1000
# number of generations
gen = 10^3
# burnin
burnin = gen/2
# rate of htz loss
gamma <- 1/3
# clonality rate
C = 0.99
# run
run <- 1
# mutation rate
mu <- 10^-4
# generations to keep for illustrating stochasticity
stocha_range <- seq((burnin+1), (gen), by=100)




# (re)initialize summary files
summary <- data.frame()


# Run simulations
source(paste0(prefix_in_path, 'simulated_generations.R'))



summary <- data.frame(Run = run,
                      ngen = gen,
                      eff_gen = length(kept_gen),
					  Clonality_rate = C,
                      N = N,
                      Htz_loss_rate = gamma,
                      Mutation_rate = mu,
                      Hosim = Hosim,
                      Hosim_sd = Hosim_sd,
                      Hesim = Hesim,
                      Hesim_sd = Hesim_sd,
                      Fsim = Fsim,
                      Fsim_sd = Fsim_sd)




# txt outputs
write.table(summary, paste0(prefix_out_path, 'summary.txt'), append = T, col.names = F, row.names = F, sep = "\t")

write.table(marker_stats[stocha_range,], paste0(prefix_out_path, 'stats_C', C, '_Y', gamma, '_mu',mu, '_run',run,
													'.txt'), sep = "\t")












