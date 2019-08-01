


if(sum(file.access(paste0("rij_",lineage,"sim_N",N,"_C",C,"_m",m,"_run",
                          nb_runs,".txt"))) != -length(nb_runs)) {
  
  
  
  tot_rij_sim_wide <- data.frame()
  tot_rij_sim_long <- data.frame()
  rij_obs <- data.frame(ncol=1)
  
  for(run in nb_runs){
    path_rij_sim <- paste0("rij_",lineage,"sim_N",N,"_C",C,"_m",m,"_run",run,".txt")
    if(file.access(path_rij_sim) != -1){
      
      rij_sim <- read.table(path_rij_sim, header = T,check.names = F)
      
      rij_sim_long <- gather(rij_sim, breaks, rij, names(rij_sim), factor_key=T)
      rij_sim_long <- data.frame(rij_sim_long, stringsAsFactors = F)
      tot_rij_sim_long <- rbind(tot_rij_sim_long, rij_sim_long)
      
      path_rij_emp <- paste0("rij_",lineage,"emp_N",N,"_C",C,"_m",m,"_run",run,".txt")
      rij_emp_tmp <- read.table(path_rij_emp, header = T,check.names = F)
      # keep rij_obs with the higher nb of bins (to be consistent with next chisq computation)
      if(length(rij_emp[,1]) < length(rij_emp_tmp[,1])) rij_emp <- rij_emp_tmp
    }
    
  } # out of run loop
  
  
  tot_rij_sim_long <- tot_rij_sim_long %>% 
    group_by(breaks) %>% 
    mutate(grouped_id = row_number())
  
  
  tot_rij_sim_wide <- tot_rij_sim_long %>% 
    spread(breaks, rij) %>% 
    select(-grouped_id)
  
  tot_rij_sim_wide[is.na(tot_rij_sim_wide)] <- 0
  
  
  
  
  # mean rij by bin over sim generations
  mean_sim <- apply(tot_rij_sim_wide, 2, mean)
  
  # if too few categories in rij_sim compare to rij_emp
  # pool the last cat of rij_emp if freq < 0.05
  if(mean_sim[1] == 0){
    rij_emp[2,] <- rij_emp[1,] + rij_emp[2,]
    rij_emp <- rij_emp[-1,]
    tot_rij_sim_wide <- tot_rij_sim_wide[,-1]
  }
  mean_sim <- apply(tot_rij_sim_wide, 2, mean)
  
  
  
  # variance over generations by bin
  var_sim <- apply(tot_rij_sim_wide, 2, var)
  
  # first step for pseudo chi2 calculation
  pre_chi2like_sim <- apply(tot_rij_sim_wide, 1, function(x) {
    (x - mean_sim)^2/var_sim
  })
  
  pre_chi2like_sim <- t(pre_chi2like_sim)
  chi2like_emp <- sum((rij_emp - mean_sim)^2/var_sim)
  chi2like_sim <- apply(pre_chi2like_sim, 1, sum)
  
  # pval calculation via chi2sim distribution
  inf_rij <- sum(chi2like_sim <= chi2like_emp)
  sup_rij <- sum(chi2like_sim >= chi2like_emp)
  pval_rij <- 2*min(inf_rij,sup_rij)/length(chi2like_sim)
  if(pval_rij == 'NA') pval_rij <- 0
  
  global_pval_rij_temp <- data.frame(N = N,
                                  Clonality_rate = C,
                                  Migration_rate = m,
                                  Lineage = lineage,
                                  pval_rij = pval_rij,
                                  stringsAsFactors = F)
  global_pval_rij <- rbind(global_pval_rij, global_pval_rij_temp)
}