
sample_size <- length(marker.all1_emp)

# Observed heterozygosity
marker_htz_obs <- matrix(NA, ncol = sample_size)
marker_htz_obs <- ifelse(marker.all1_emp == marker.all2_emp, 0, 1)
Hi_emp <- mean(marker_htz_obs[1,])


# convert to matrix
marker.all1_emp <- data.matrix(marker.all1_emp)
marker.all2_emp <- data.matrix(marker.all2_emp)

# Compute allelic frequencies
frequency_table_emp <- allele_frequency_func(c(marker.all1_emp,marker.all2_emp))

# Expected heterozygosity estimation (correction for sample size using Nei 1983 formula)
Hs_emp_uncorrected <- 1 - sum(frequency_table_emp^2)
Hs_emp <- (sample_size/(sample_size - 1)) * (Hs_emp_uncorrected - (Hi_emp/(2*sample_size)))

# Individual inbreeding coefficient
Fis_emp <- 1-(Hi_emp/Hs_emp)

marker_stats_emp <- cbind(Hi_emp,Fis_emp,Hs_emp) 
colnames(marker_stats_emp) <- c('Hi','Fis','Hs')

# Marker properties
marker_cumfreq <- cumsum(frequency_table_emp)
marker_prop <- data.frame(alleles = as.numeric(names(frequency_table_emp)), freq = c(frequency_table_emp), 
                     cumfreq = marker_cumfreq, row.names = NULL)

