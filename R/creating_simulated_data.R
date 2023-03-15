# generate fake data
library(tidyselect)
library(tidyverse)
library(data.table)
library(devtools)
#compile
library(rstan)
perc_chrom = rbeta(100,2.5,2.5)# generates 100 protein samples with some % of pure chromatin (enrichment)
cell_types = rep(c(1:50), times = 2) #two replicates of each cell line
final_df = data.frame(perc_chrom = perc_chrom) 

N_prot= 20 #twenty proteins in dummy test set
Total_prot = rnorm(N_prot,20,1) #total abundance of each protein (chrom+non-chrom)
Prot_chrom = c(
    runif(N_prot, min = 5, max = 16) # general protein abundance in pure (hypotherical) chromatin samples

)
#making a chromatin and cytoplasm (non-chrom) matrix 
#col = proteins, row = cell lines, 
#if a protein is more chromatin only in a specific cell line, it will have a positive value in this matrix

chrom_prot_cellline = matrix(rep(0,length(unique(cell_types))* length(Prot_chrom)), #making a 0 matrix with rows T of tissue and Col P of PRoteins
                           nrow = length(unique(cell_types)), ncol = length(Prot_chrom))
cyto_prot_cellline = matrix(rep(0,length(unique(cell_types))* length(Prot_chrom)), #making a 0 matrix with rows T of tissue and Col P of PRoteins
                           nrow = length(unique(cell_types)), ncol = length(Prot_chrom))

#making some proteins cell line specific
chrom_prot_cellline[1,1] = 4
cyto_prot_cellline[1,1] = -3
chrom_prot_cellline[1,20] = -2
cyto_prot_cellline[1,20] = 3


#the cytoplasmic/non-chrom fraction of the protein is total-chrom
Prot_cyto = Total_prot - Prot_chrom
for(i in 1:length(Prot_chrom)){
    Prot_id = glue::glue("Prot_{i}")
    Pr_cyto = glue::glue("Pr_cyto_{i}")
    Pr_nucl = glue::glue("Pr_nucl{i}")
    test_dat = data.frame(perc_chrom = perc_chrom) 
    test_dat_tmp = test_dat|> 
        mutate(Tissue = tissue_type) |> 
        mutate(
          
            !!Pr_nucl  := rnorm( n() , mean= Prot_chrom[i] +chrom_prot_tissue[Tissue,i]  , sd=0.2 )
            ,
            !!Pr_cyto := rnorm( n(), mean=Prot_cyto[i] +cyto_prot_tissue[Tissue,i]  , sd=0.2 )
            ,
            #the OBSERVED abundance of the protein in the sample, depends on how much that protein is on pure chromatin and what perce of that sample is pure chromatin
            !!Prot_id  := .data[[Pr_cyto]]*(1-perc_chrom)+ .data[[Pr_nucl]]*perc_chrom
        ) |> as.data.table()
    final_df  = cbind(final_df,test_dat_tmp)
}


test_dat = dplyr::select(final_df, matches("Prot")) |> as.matrix()
library(data.table)
test_dat <- as.data.table(test_dat)


#check stan is install3d

model_orig <-  stan_model('multiple_prots_linear_logistic_enrichment_data_hierarchical.stan')


#pass data to stan and run model
parallel::detectCores()
options(mc.cores = 12)
# init_fun <- function(...) list(b=0.1)
input_enrichment = scale((final_df$perc_chrom*2.8) +3) #in my actual experiment i don't have the true percentages, which I need to estimate
#so I'm taking the true percentages and transforming the so that the model needs to find the original percentages of the simulation
input_enrichment = input_enrichment[,1]
matrix_input = test_dat
fit3 <- sampling(model_orig, list(N= nrow(matrix_input),
                             P = ncol(matrix_input),
                             input_enrichment =input_enrichment,
                             prot_intensity = matrix_input,
                             cell_line= cell_types,
                             cell_lines= 50), 
                 iter = 1000, chains = 4
                 # , init = init_fun
)
#diagnose
shinystan::launch_shinystan(fit3)
print(fit3)
params <- rstan::extract(fit3)