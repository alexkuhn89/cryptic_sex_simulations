


prefix_in_path 	<- ""
prefix_out_path <- ""

library('minqa')

# Load functions for subsequent scripts
source(paste0(prefix_in_path,'functions_simulations.R'))

# Read file with Loiselle's kinships (output Spagedi) to compute marker's recombination and mutation rates
empirical_data_file <- "empirical_data.csv"
empirical_data_table <- read.csv(paste0(prefix_in_path, empirical_data_file), header = T)

# Define color for rij distributions
color1 <- rgb(1,0,0,alpha = 0.4) # red 40% opacity
color2 <- rgb(0,0,1,alpha = 0.4) # blue 40% opacity





# Model parameters to define for generations simulation

# Population size (estimated effective size of the whole population)
N = 250
# number of generations
gen = 5000
# burnin
burnin = gen/2
# migration rate
m = 0
m = m/100
# microsatellite motif size in bp
microsat_motif_size = 2
# generations to compute rij
rij_seq <- seq(2510,5000,by=10)


# lineages names
lineage_list <- c("R1", "R2")
# markers names
marker_list <- c("CH12", "CH08", "CH01", "CC54", "CH11", "CC11", "CH22")




# 'Clonality rate' (rate of parthenogenetic production of new queens)
C = 0.95


run = 1
print(paste0("clonality_rate = ", C, " run",run))


# (Re)initialize summary files
summary_file <- data.frame()
Fis_sim_tot <- as.data.frame(matrix(nrow = gen, ncol = 0))
Hs_sim_tot <- as.data.frame(matrix(nrow = gen, ncol = 0))




for (lineage in lineage_list) {
  
  
  # general 'mother choice' for all markers
  mother_reprod <- matrix(NA, nrow = gen, ncol = N)
  mother_reprod <- apply(mother_reprod,1, function (x) ifelse(runif(N) > C, ceiling(runif(N)*N), 
                                                              -ceiling(runif(N)*N))) 
  mother_reprod <- t(mother_reprod)
  
  
  for (marker in marker_list) {
    # read empirical allelic data
    path_all1 <- paste0(prefix_in_path, marker, lineage, "_all1.txt")
    path_all2 <- paste0(prefix_in_path, marker, lineage, "_all2.txt")
    marker.all1_emp <- read.table(path_all1, header=T)
    marker.all2_emp <- read.table(path_all2, header=T)
    
    marker_Fis_empirical <- empirical_data_table[empirical_data_table$marker == paste0(marker,lineage), "Fis"]
    marker_Hs_empirical <- empirical_data_table[empirical_data_table$marker == paste0(marker,lineage), "Hs"]
    
    # extract allelic frequencies of 'unassigned males' for migration rate (see main text)
    migration_path <- paste0(prefix_in_path, marker, lineage, "_mig_freq.txt")
    migration_freq <- read.table(migration_path, header = T, check.names = F)
    migration_cumfreq <- cumsum(c(migration_freq[1,]))
    
    
    # Generate empirical data
    source(paste0(prefix_in_path,'empirical_data.R'))
    
    
    # Run simulations
    source(paste0(prefix_in_path,'Simulation_generations.R'))
    
    
    # Data analyses
    source(paste0(prefix_in_path,'Data_analyses.R'))
    # main results in 'marker_stats', 'pval_F' and 'pval_He'
    
    
    # save F evolution across generations
    Fis_sim_current <- data.frame(marker_stats[,'Fis'])
    colnames(Fis_sim_current) <- c(paste0(marker,lineage))
    Fis_sim_tot <- cbind(Fis_sim_tot, Fis_sim_current)
    
    
    # save He evolution across generations
    Hs_sim_current <- data.frame(marker_stats[,'Hs'])
    colnames(Hs_sim_current) <- c(paste0(marker,lineage))
    Hs_sim_tot <- cbind(Hs_sim_tot, Hs_sim_current)
    
    
    
    # Rij computation (coeff. after Li et al. 1993)
    source(paste0(prefix_in_path,'Rij.R'))
    
    summary_out_current <- data.frame(Run = run,
                                      ngen = gen,
                                      Clonality_rate = C,
                                      N = N,
                                      Lineage = lineage,
                                      Marker = marker,
                                      Recomb_rate = gamma,
                                      Mutation_rate = mu,
                                      Migration_rate = m,
                                      Hs_emp = marker_stats_emp[,'Hs'],
                                      Hs_sim = Hs_sim_mean,
                                      Fis_emp = marker_stats_emp[,'Fis'],
                                      Fis_sim = Fis_sim_mean)
    
    summary_file <- rbind(summary_file, summary_out_current)
    
    
    
  } # out of markers loop
  
  # compute rij for empirical data by lineage (over all loci)
  rij_emp <- sum_rijl_emp/sum_wijl_emp
  name_rij_emp <- paste0('rij_emp_',lineage)
  assign(name_rij_emp,rij_emp)
  
  
  
  # compute rij for simulated data for by generation and by lineage (over all loci)
  for (generation in 1:length(rij_seq)) {
    rij_gen[generation,] <- sum_rijl_gen[generation,]/sum_wijl_gen[generation,]
  }
  
  name_rij_gen <- paste0('rij_gen_',lineage)
  assign(name_rij_gen,rij_gen)
  
  remove(rij_gen, sum_rijl_gen, sum_wijl_gen, sum_rijl_emp, sum_wijl_emp)
  
  
  source(paste0(prefix_in_path, 'rij.R'))
  
  
  name_binned_rij_emp <- paste0('binned_rij_emp_',lineage)
  assign(name_binned_rij_emp,binned_rij_emp)
  
  name_binned_rij_sim <- paste0('binned_rij_sim_',lineage)
  assign(name_binned_rij_sim,binned_rij_sim)
  
  
  
} # out of lineage loop


## OUTPUTS

# Fis and Hs over generations
write.table(Fis_sim_tot, paste0(prefix_out_path,"Fis_sim_tot_N",N,"_C",C,"_m",m,"_run", run, ".txt"),sep="\t")
write.table(Hs_sim_tot, paste0(prefix_out_path,"Hs_sim_tot_N",N,"_C",C,"_m",m,"_run", run, ".txt"),sep="\t")

# Binned rij distributions
for (lineage in lineage_list) {
  binned_rij_sim_lineage <- get(paste0("binned_rij_sim_",lineage))
  binned_rij_emp_lineage <- get(paste0("binned_rij_emp_",lineage))
  bins <- round(seq(1-(length(binned_rij_sim_lineage[1,])-1)/10,1,0.1),digits=1)
  colnames(binned_rij_sim_lineage) <- bins
  write.table(binned_rij_sim_lineage, paste0(prefix_out_path,"rij_",lineage,"sim_N",N,"_C", C,"_m",m,"_run",run,".txt"), sep="\t")
  write.table(binned_rij_emp_lineage, paste0(prefix_out_path,"rij_",lineage,"emp_N",N,"_C", C,"_m",m,"_run",run,".txt"), sep="\t")
}

# summary of the simulations' parameters
write.table(summary_file, paste0(prefix_out_path,"summary.txt"), append = T, col.names = F, row.names = F, sep = "\t")










  
  
  







