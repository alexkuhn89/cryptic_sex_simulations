# cryptic_sex_simulations
R code for simulations of "Detection of cryptic sex in automictic populations: theoretical expectations and a case study in *Cataglyphis* desert ants".
---------------------------------------------------
## Files : ` CXXXRX_allX.txt `

Allelic data for microsatellite marker CcXX or ChXX of *Cataglyphis mauritanica* queens from the study of Kuhn et al. (Kuhn, A., Bauman, D., Darras, H., & Aron, S. (2017). Sex-biased dispersal creates spatial genetic structure in a parthenogenetic ant with a dependent-lineage reproductive system. Heredity, 119(4), 207-213.)

`R1` and `R2` corresponds to the lineage of the queens. `all1` and `all2` corresponds to the allele (first or second; queens are diploids in ants).

## Files : ` CXXXRX__mig_freq.txt `

Allelic frequencies of the putative migrant males from the study of Kuhn et al. 2017.

## File : ` Data_analyses_mauritanica.R `

R script to compute observed and expected heterozygosity and individual inbreeding coefficient on the *C. mauritanica* simulated data.

## File : ` Post_simulations_treatments_mauritanica.R `

Generate p-values by comparing empirical to simulated data over loci and over runs for the *C. mauritanica* based dataset.

## File : ` Simulation_generations_mauritanica.R `

R script to simulate new generation of C. mauritanica queens according to their reproductive mode (rate of automixis and rate of sexual reproduction) and based on empirical data from Kuhn et al. 2017.

## File : ` chisquare-like_rij_mauritanica.R `

R script to compute the chisquare-like statistics based on *r*<sub>ij</sub> to compare empirical and simulated datasets.

## File : ` empirical_data.csv `

Loiselle's kinships (Fis) and expected heterozygosity (Hs) from empirical queen genotypes computed with SPAGeDi.

## File : ` empirical_data_mauritanica.R `

R script to compute observed and expected heterozygosity and individual inbreeding coefficient on the *C. mauritanica* empirical data.

## File : ` functions_simulation.R `

R script with the functions used for simulating generations (in the general model or for the simulated *C. mauritanica* dataset).

## File : ` param_script.R `

R script to introduce simulation parameters and to run generations simulation (through the script `simulated_generations.R`) for the general model.

## File : ` param_script_mauritanica.R `

R script to introduce simulation parameters and to run generations simulation (through the script ` Simulation_generations_mauritanica.R `) for the *C. mauritanica* simulated dataset.

## File : ` rij_mauritanica.R `

## File : ` simulated_generations.R `

R script to run generations simulation for the general model.
