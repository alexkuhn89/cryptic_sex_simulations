
### Analysis of simulated data


marker_htz <- matrix(NA, nrow = gen, ncol = sample_size)


# Loop for heterozygotes counting 
for (i in 1:gen) {
  for (j in 1:sample_size) {
    if(marker.all1_gen[i,j] == marker.all2_gen[i,j]) {
      marker_htz[i,j] = 0
    } else marker_htz[i,j] = 1
  }
}


# Compute observed heterozygosity for each generation
Hi <- t(apply(marker_htz, MARGIN = 1, FUN = mean))
Hi <- t(Hi)

# Compute allelic frequencies for each generation (correction for sample size using Nei 1983 formula)
pooled_alleles_gen <- cbind(marker.all1_gen[,1:sample_size], marker.all2_gen[,1:sample_size])
frequency_table <- apply(pooled_alleles_gen, 1, allele_frequency_func)

# Compute expected heterozygosity for each generation
Hs_uncorrected <- lapply(frequency_table, function(x) {
  result <- 1-sum(x^2)
  return(result)
})
Hs_uncorrected <- unlist(Hs_uncorrected)

Hs <- (sample_size/(sample_size - 1)) * (Hs_uncorrected - (Ho/(2*sample_size)))




Fis <- 1-(Hi/Hs)

marker_stats <- cbind(Hs,Fis)
colnames(marker_stats) <- c('Hs','Fis')

# Mean and SD of Fis and Hs per run (for informative purpose)
kept_gen <- which(!is.na(marker_stats[,'Fis']))
Fis_sim_mean <- mean(marker_stats[kept_gen[kept_gen > burnin],'Fis'])
Fis_sim_sd <- sd(marker_stats[kept_gen[kept_gen > burnin],'Fis'])
Hs_sim_mean <- mean(marker_stats[kept_gen[kept_gen > burnin],'Hs'])
Hs_sim_sd <- sd(marker_stats[kept_gen[kept_gen > burnin],'Hs'])