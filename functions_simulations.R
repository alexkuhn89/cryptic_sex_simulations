#FUNCTIONS sex-frequency simulations
#########################

# function to read microsatellites data files from the prefix (eg. 'CH22R1')
method_extract_obs_data <- function(prefix_allele_file, queen_lineage,prefix) {
  stopifnot(is.character(prefix_allele_file))
  # paste0 to couple string objects
  pattern_file <- paste0(prefix, prefix_allele_file, queen_lineage, "_all1.txt")
  # create lists with the file names and stop if more than one file or if it doesn't exist
  file_all1 <- list.files(pattern = pattern_file) ; stopifnot(length(file_all1) == 1)
  # gsub to replace characters
  file_all2 <- gsub("all1", "all2", file_all1)    ; stopifnot(file.exists(file_all2))
  
  result <- list()
  result$all1 <- read.table(file_all1,header = T)
  result$all2 <- read.table(file_all2,header = T)
  return(result)
}




# Generate table of allelic frequencies
allele_frequency_func <- function(alleles)  {
  alleles_name <- sort(unique(alleles[alleles != 0]))
  frequency <- table(alleles[alleles != 0]) / length(alleles[alleles != 0])
}




# function to generate the 'mother choice' (sex or parthenogenesis - with or without recombination)
# for each marker according to its recombination rate
# sex if ]0, N] ; parth. without recombi if [-N, 0] ; parth with recombi if [-2N, -N[ 
marker_recomb_func <- function(mother_reprod_index,
                               htz_loss_rate,
                               n = N) {
  if(mother_reprod_index > 0) {
    marker_reprod <-  mother_reprod_index
  } else if(runif(1) < htz_loss_rate) {
    marker_reprod <-  -n -abs(mother_reprod_index)
  } else {
    marker_reprod <- -abs(mother_reprod_index)
  }
  
  return(marker_reprod)
}


# Input function for 'bobyqa' to optimize 'mu' and 'gamma' for each marker
mu_gamma_func <- function(param, n=N, Hs_emp= marker_Hs_empirical, Fis_emp=marker_Fis_empirical, c=C) {
  # Hs estimation using eq. 4 from main text
  Hs_IAM <- (param[1]*(2-param[1])*(2*n+c*(1-2*param[1])*(2*n*(1-param[2])-1)))/(1+2*n*param[1]*(2-param[1])+c*(1-2*param[1])*(1-2*param[2]+2*n*param[1]*(2-param[1])*(1-param[2]))) 
  # Conversion of the IAM mutation model to a SMM mutation model using eq. 10 from main text
  Hs_SMM <- 1 - 1/((1+Hs_IAM)/(1-Hs_IAM))^(1/2)
  # Hi estimation using eq. 2 from main text
  Hi <- (2*param[1] + Hs_IAM*(1-c)*(1-1/n)*(1-2*param[1]))/(1-(1-2*param[1])*(c*(1-param[2])+(1-c)/(2*n)))
  
  Fcalc <- 1-Hi/Hs_IAM
  result <- sum((Hs_SMM - Hs_emp)^2 + (Fcalc - Fis_emp)^2)
}


# function to generate stepwise mutation in microsatellite markers
mutation_func <- function(allele, mutation_rate = mu, motif_size = microsat_motif_size) {
  if(runif(1) < mutation_rate) {
    new_allele <- c(allele - motif_size, allele + motif_size)
    result <- new_allele[sample(1:2,1)]
  } else {
    result <- allele
  }
  return(result)
}



# Function to generate the first ('left') allele for each individual of a new generation based on the 
# 'mother choice' index of a given marker. Alleles can mutate under a stepwise mutation model.
left_allele_mu_func <- function(marker_reprod_index,
                                left_alleles,
                                right_alleles,
                                n = N,
                                mutation_rate) {
  if(marker_reprod_index > 0) { # sex
    if(runif(1)<0.5) {
      result <- mutation_func(left_alleles[marker_reprod_index])
    }
    else {
      result <- mutation_func(right_alleles[marker_reprod_index])
    }
  } else if(marker_reprod_index >= -n) { # partheno without recombi
    result <- mutation_func(left_alleles[abs(marker_reprod_index)])
  } else if(runif(1)<0.5) {  # partheno with recombi
    result <- mutation_func(left_alleles[abs(marker_reprod_index)-n])
  } else {
    result <- mutation_func(right_alleles[abs(marker_reprod_index)-n])
  }
  return(result)
}





# Function to generate the second ('right') allele for each individual of a new generation based on the 
# 'mother choice' index of a given marker. Alleles can mutate under a stepwise mutation model. 
# When reproduction is sexual, alleles may come from migrant males (freq of unassigned males; see main text)
# with a probability m (the migration rate).
right_allele_mu_func <- function(marker_reprod_index,
                                 migration_cumfreq,
                                 left_alleles,
                                 right_alleles,
                                 n = N,
                                 mutation_rate,
                                 migration_rate = m) {
  if(marker_reprod_index > 0) {
    # probability m of migration event 
    if(runif(1)< migration_rate) {
      index <- findInterval(runif(1),migration_cumfreq)+1
      result <- as.numeric(names(migration_cumfreq[index]))
    } else {
      # if sexual reproduction without migration, second allele randomly chosen in the precedent generation
      male <- ceiling(runif(1)*n)
      # resample to avoid selfing
      if(male == marker_reprod_index) male <- ceiling(runif(1)*n)
      if(runif(1)<0.5) {
        
        result <- mutation_func(left_alleles[male])
      }
      else {
        result <- mutation_func(right_alleles[male])
      }
    } 
  } else if(marker_reprod_index >= -n) {
    result <- mutation_func(right_alleles[abs(marker_reprod_index)])
  } else {
    result <- c('duplicate')
  }
  return(result)
}



# Function to combine functions for left and right alleles because if recombination occurs, the right
# allele depends on the value of the left allele
both_alleles_func <- function(marker_reprod_index,
                              migration_cumfreq,
                              left_alleles,
                              right_alleles,
                              n=N) {
  la <- left_allele_mu_func(marker_reprod_index,
                            left_alleles,
                            right_alleles,
                            n = N)
  ra <- right_allele_mu_func(marker_reprod_index, 
                             migration_cumfreq,
                             left_alleles,
                             right_alleles,
                             n = N)
  if( ra == "duplicate") {
    ra <- la
  }
  results <- c(la,ra)
  return(results)
}


Hs_corrected_func <- function(allele_freq, sample_size, Hi) {
  Hs_uncorrected <- 1 - sum(allele_freq^2)
  Hs <- (sample_size/(sample_size - 1)) * (Hs_uncorrected - (Hi/(2*sample_size)))
  return(Hs)
}


# function for similarity coefficient sij
sij_func <- function(alleles,
                     sample_comparisons) {
  index_pw_comparison <- t(combn(c(1:sample_comparisons), 2))
  sapply(1:nrow(index_pw_comparison), function(x) {
    index_i <- index_pw_comparison[x, 1]
    index_j <- index_pw_comparison[x, 2]
    i <- c(alleles[index_i],alleles[index_i+sample_comparisons])
    j <- c(alleles[index_j],alleles[index_j+sample_comparisons])
    res <- sum(i %in% j, j %in% i)/4
    return(res)
  })
}
