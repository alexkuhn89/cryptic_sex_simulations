

current_rij_sim <- get(paste0('rij_gen_',lineage))
current_rij_sim <- na.omit(current_rij_sim)
# round rij values < -1.5 to -1.5 to avoid subsequent issues due to bin with very low representation
current_rij_sim <- replace(current_rij_sim, current_rij_sim < -1.5, -1.5)

current_rij_emp <- get(paste0('rij_emp_',lineage))

# define the minimal bin (over both obs and sim rij)
min_break <- min(floor(min(current_rij_sim)*10)/10, floor(min(current_rij_emp)*10)/10)
# define the bins
breaks <- seq(min_break, 1, by = 0.1)

# mean rij by bin over sim generations
mean_sim <- hist(current_rij_sim,plot=F, breaks=breaks)$counts/sum(hist(current_rij_sim,plot=F)$counts)
# obs rij divided with the new bins
binned_rij_emp <- hist(current_rij_emp,plot=F, breaks = breaks)$counts/sum(hist(current_rij_emp,
                                                                                plot=F, breaks = breaks)$counts)

# sim rij divided with the new bins for each generation
binned_rij_sim <- apply(current_rij_sim, 1, function(x) {
  hist(x,plot=F, breaks = breaks)$counts/sum(hist(x,plot=F, breaks = breaks)$counts)
})
binned_rij_sim <- t(binned_rij_sim)

