

marker_reprod_mode <- matrix(NA, nrow = gen, ncol = N)
marker.all1_gen <- matrix(NA, nrow = gen, ncol = N)
marker.all2_gen <- matrix(NA, nrow = gen, ncol = N)


# Starting parameter values for mu and gamma, respectively, for the 'bobyqa' function
start_par <- c(0.001,0.1)

# Iterative search of optimal values of mu and gamma according to the empirical values of Hi, Hs and Fis 
mu_gamma <- bobyqa(start_par, mu_gamma_func,lower=c(1/100000,0),upper = c(0.005,1/3),control = list(rhobeg=0.0000001))
mu <- mu_gamma$par[1]
gamma <- mu_gamma$par[2]


# Alleles randomly assigned to first generation individuals according to the empirical allelic frequencies
# (each row corresponds to a generation, each column to an individual)
allele_index <- findInterval(runif(N),marker_prop$cumfreq)+1
marker.all1_gen[1,] <- marker_prop$alleles[allele_index]
allele_index <- findInterval(runif(N),marker_prop$cumfreq)+1
marker.all2_gen[1,] <- marker_prop$alleles[allele_index]


# Loop to simulate new generations

for (i in 2:gen) {
  
  # Mother reproduction for according to the marker's rate of heterozygosity loss
  marker_reprod_mode[i,] <- sapply(mother_reprod[i,], function(x) marker_recomb_func(x, 
                                                                                     htz_loss_rate = gamma))
  
  marker.all1_2.nextgen <- sapply(marker_reprod_mode[i,], 
                                function(x) both_alleles_func(x,
                                                              migration_cumfreq = migration_cumfreq,
                                                              left_alleles = marker.all1_gen[i-1,],
                                                              right_alleles = marker.all2_gen[i-1,]))
  
   # results of 'both_alleles_func' put the first alleles on the first line and the second alleles on the second line
  
   marker.all1_gen[1:i,] <- rbind(marker.all1_gen[1:i-1,], marker.all1_2.nextgen[1,])
   marker.all2_gen[1:i,] <- rbind(marker.all2_gen[1:i-1,], marker.all1_2.nextgen[2,])
  
}






