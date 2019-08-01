
## Initial generation

# random generation of alleles
marker.all1_gen <- ceiling(runif(N)*N)
marker.all2_gen <- ceiling(runif(N)*N)



marker_htz <- matrix(NA, nrow = 1, ncol = sample_size)

# Heterozygosity calculation
for (j in 1:sample_size) {
  if(marker.all1_gen[j] == marker.all2_gen[j]) {
    marker_htz[j] = 0
  } else marker_htz[j] = 1
}

# Observed heterozygosity estimation
Ho <- mean(marker_htz)

# Table of allelic frequencies
pooled_alleles_gen <- cbind(marker.all1_gen[1:sample_size], marker.all2_gen[1:sample_size])
frequency_table <- allele_frequency_func(pooled_alleles_gen)

# Expected heterozygosity estimation
He <- 1-sum(frequency_table^2)
He <- (sample_size/(sample_size - 1)) * (He - (Ho/(2*sample_size)))

# Individual inbreeding coefficient estimation
Fcoef <- 1-(Ho/He)

marker_stats <- cbind(Ho,Fcoef,He) 
colnames(marker_stats) <- c('Ho','F','He')

## Loop for generating new generations

for (i in 2:gen) {
  
  # Mother reproduction (sexual, or parthenogenetic with or without heterozygosity loss according to gamma) 
  mother_reprod <- ifelse(runif(N) > C, ceiling(runif(N)*N), -ceiling(runif(N)*N))
  marker_reprod_mode <- sapply(mother_reprod, function(x) marker_recomb_func(x, 
                                                                             htz_loss_rate = gamma))
  
  
  # Alleles of the next generation according to the 'mother reproduction' and to potential mutations under an IAM model
  marker.all1_2.nextgen <- sapply(marker_reprod_mode, 
                                  function(x) both_alleles_IAM_func(x,
                                                                    left_alleles = marker.all1_gen,
                                                                    right_alleles = marker.all2_gen))
  
  
  
  marker.all1_nextgen <- marker.all1_2.nextgen[1,]
  marker.all2_nextgen <- marker.all1_2.nextgen[2,]
  
  marker_htz <- matrix(NA, nrow = 1, ncol = sample_size)
  
  
  for (j in 1:sample_size) {
    if(marker.all1_nextgen[j] == marker.all2_nextgen[j]) {
      marker_htz[j] = 0
    } else marker_htz[j] = 1
  }
  
  Ho <- mean(marker_htz)
  
  pooled_alleles_gen <- cbind(marker.all1_nextgen[1:sample_size], marker.all2_nextgen[1:sample_size])
  frequency_table <- allele_frequency_func(pooled_alleles_gen)
    
  He <- 1-sum(frequency_table^2)
  He <- (sample_size/(sample_size - 1)) * (He - (Ho/(2*sample_size)))
  
  Fcoef <- 1-(Ho/He)
    
  marker_stats <- rbind(marker_stats,c(Ho,Fcoef,He))
    
  marker.all1_gen <- marker.all1_nextgen
  marker.all2_gen <- marker.all2_nextgen
}



marker_stats <- as.data.frame(marker_stats)


# keep only generations without 'NA' values (and after burnin generations) for 
# the estimation of the mean of each statistic 
kept_gen <- which(!is.na(marker_stats[,'F']))

# mean Fsim 
Fsim <- mean(marker_stats[kept_gen[kept_gen > burnin],'F'])
Fsim_sd <- sd(marker_stats[kept_gen[kept_gen > burnin],'F'])

# mean He
Hesim <- mean(marker_stats[kept_gen[kept_gen > burnin],'He'])
Hesim_sd <- sd(marker_stats[kept_gen[kept_gen > burnin],'He'])

# mean Ho
Hosim <- mean(marker_stats[kept_gen[kept_gen > burnin],'Ho'])
Hosim_sd <- sd(marker_stats[kept_gen[kept_gen > burnin],'Ho'])






