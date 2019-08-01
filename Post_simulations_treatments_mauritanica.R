

library(EmpiricalBrownsMethod)
library(tidyr)
library(dplyr)



# Generating pvalues by comparing empirical to simulated data over loci and over runs

# Define the number of runs simulated per parameter set 
nb_runs <- c(1:2)



# Read the simulations' parameters
summary_all <- read.table("summary.txt", header=F)
colnames(summary_all) <- c("Run",
                           "ngen",
                           "Clonality_rate",
                           "N",
                           "Lineage",
                           "Marker",
                           "Recomb_rate",
                           "Mutation_rate",
                           "Migration_rate",
                           "Hs_emp",
                           "Hs_sim",
                           "Fis_emp",
                           "Fis_sim")
 
detailed_pval <- data.frame()
global_pval <- data.frame()
global_pval_rij <- data.frame()

# Check the simulated parameter sets 
N_range <- sort(unique(summary_all$N))
clonality_range <- sort(unique(summary_all$Clonality_rate))
m_range <- sort(unique(summary_all$Migration_rate))
lineage_list <- paste0(unique(summary_all$Lineage))
marker_list <- paste0(unique(summary_all$Marker))


for(N in N_range) {
  for(C in clonality_range) {
    for(m in m_range) {
      
      if(sum(file.access(paste0("Fis_sim_tot_N",N,"_C",C,"_m",m,"_run",nb_runs,".txt"))) != -length(nb_runs)) {
        
        for(lineage in lineage_list){
          
          ebm_matrix_Fis_global <- NULL
          ebm_matrix_Hs_global <- NULL
          
          for(run in nb_runs){
            
            path_Fis <- paste0("Fis_sim_tot_N",N,"_C",C,"_m",m,"_run",run,".txt")
            path_Hs <- paste0("Hs_sim_tot_N",N,"_C",C,"_m",m,"_run",run,".txt")
            
            if(file.access(path_Fis) != -1 & file.access(path_Hs) != -1){
              Fis_sim <- read.table(path_Fis)
              Hs_sim <- read.table(path_Hs)
              
              exclude <- NULL
              
              for(marker in marker_list){
                # exclude generations with NA values for Fis (due to fixed alleles)
                exclude <- c(exclude, which(is.na(Fis_sim[,paste0(marker, lineage)])))
                
              } # out of marker list
              
              exclude <- unique(exclude)
              exclude <- exclude[exclude > burnin]
              
              ebm_matrix_Fis <- t(Fis_sim[,which(substr(colnames(Fis_sim),5,6) == lineage)])
              ebm_matrix_Hs <- t(Hs_sim[,which(substr(colnames(Fis_sim),5,6) == lineage)])
              
              if(length(exclude) > 0) {
                ebm_matrix_Fis <- ebm_matrix_Fis[,-exclude]
                ebm_matrix_Hs <- ebm_matrix_Hs[,-exclude]
              }
              
              ebm_matrix_Fis <- ebm_matrix_Fis[,(burnin+1):length(ebm_matrix_Fis[1,])]
              ebm_matrix_Hs <- ebm_matrix_Hs[,(burnin+1):length(ebm_matrix_Hs[1,])]
              
              ebm_matrix_Fis_global <- cbind(ebm_matrix_Fis_global, ebm_matrix_Fis)
              ebm_matrix_Hs_global <- cbind(ebm_matrix_Hs_global, ebm_matrix_Hs)
              
              
            }
          } # out of run loop
            
            
          
          
          for(marker in marker_list){
            
            # mean and sd for simulated Fis values over kept generations (over runs)
            mean_Fis_sim <- mean(ebm_matrix_Fis_global[paste0(marker,lineage),])
            sd_Fis_sim <- sd(ebm_matrix_Fis_global[paste0(marker,lineage),])
            
            # pval for Fis
            Fis_emp <- mean(summary_all[summary_all$Lineage == lineage & summary_all$Marker == marker, 'Fis_emp'])
            
            inf_Fis <- sum(ebm_matrix_Fis_global[paste0(marker,lineage),] >= Fis_emp)
            sup_Fis <- sum(ebm_matrix_Fis_global[paste0(marker,lineage),] <= Fis_emp)
            pval_Fis <- 2*min(c(inf_Fis, sup_Fis))/length(ebm_matrix_Fis_global[1,])
            
            # mean and sd for simulated Hs values over kept generations (over runs)
            mean_Hs_sim <- mean(ebm_matrix_Hs_global[paste0(marker,lineage),])
            sd_Hs_sim <- sd(ebm_matrix_Hs_global[paste0(marker,lineage),])
            
            # pval for Hs
            Hs_emp <- mean(summary_all[summary_all$Lineage == lineage & summary_all$Marker == marker, 'Hs_emp'])
            
            inf_Hs <- sum(ebm_matrix_Hs_global[paste0(marker,lineage),] >= Hs_emp)
            sup_Hs <- sum(ebm_matrix_Hs_global[paste0(marker,lineage),] <= Hs_emp)
            pval_Hs <- 2*min(c(inf_Hs, sup_Hs))/length(ebm_matrix_Hs_global[1,])
            
            recomb_rate <- mean(summary_all[summary_all$Lineage == lineage & summary_all$Marker == marker, 'Recomb_rate'])
            mutation_rate <- mean(summary_all[summary_all$Lineage == lineage & summary_all$Marker == marker, 'Mutation_rate'])
            
            detailed_pval_temp <- data.frame(N =N,
                                         Clonality_rate = C,
                                         Migration_rate = m,
                                         Lineage = lineage,
                                         Marker = marker,
                                         Recomb_rate = recomb_rate,
                                         Mutation_rate = mutation_rate,
                                         sampled_gen = length(ebm_matrix_Fis_global[paste0(marker,lineage),]),
                                         Hs_emp = Hs_emp,
                                         mean_Hs_sim = mean_Hs_sim,
                                         sd_Hs_sim = sd_Hs_sim,
                                         pval_Hs = pval_Hs,
                                         Fis_emp = Fis_emp,
                                         Mean_Fis_sim = mean_Fis_sim,
                                         sd_Fis_sim = sd_Fis_sim,
                                         pval_Fis = pval_Fis,
                                         stringsAsFactors = F)
            
            detailed_pval <- rbind(detailed_pval,detailed_pval_temp)
            
            
            
            
            
          } # out of marker loop (2nd)
          
          #prepare data for global pval
          ebm_pval_Fis <- as.numeric(detailed_pval[detailed_pval$Clonality_rate == C &
                                               detailed_pval$N == N &
                                               detailed_pval$Migration_rate == m &
                                               detailed_pval$Lineage == lineage 
                                             , 'pval_Fis'])
          
          ebm_pval_Hs <- as.numeric(detailed_pval[detailed_pval$Clonality_rate == C &
                                                detailed_pval$N == N &
                                                detailed_pval$Migration_rate == m &
                                                detailed_pval$Lineage == lineage 
                                              , 'pval_Hs'])
          
          
          # combine pval over markers (same size, dependent pval)
          global_pval_Fis <- empiricalBrownsMethod(data_matrix = ebm_matrix_Fis_global, p_values = ebm_pval_Fis, extra_info = F)
          global_pval_Hs <- empiricalBrownsMethod(data_matrix = ebm_matrix_Hs_global, p_values = ebm_pval_Hs, extra_info = F)
          
          global_pval_temp <- data.frame(N = N,
                                        Clonality_rate = C,
                                        Migration_rate = m,
                                        Lineage = lineage,
                                        global_pval_Hs = global_pval_Hs,
                                        global_pval_Fis = global_pval_Fis, 
                                        stringsAsFactors = F)
          
          global_pval <- rbind(global_pval, global_pval_temp)
          
          
          ### Compute pval for rij (through a pseudo-chisquare test; see main text) over runs
          source("chisquare-like_rij.R")
          
          
          
          
        } # out of lineage loop
      }
      
      
      
    }
  }
}


write.table(detailed_pval, "detailed_pval.txt", sep = "\t", row.names = F)

if(length(global_pval[,1]) == length(global_pval_rij[,1])) {
  global_pval <- cbind(global_pval, pval_Rij = global_pval_rij$pval_rij)
  write.table(global_pval, "global_pval.txt", sep = "\t", row.names = F)
} else print("incongruences between rij and Fis/Hs pvalues")
       




