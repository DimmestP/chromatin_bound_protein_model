data {
  int<lower=0> cell_lines; //number of different cell_lines
  int<lower=0> proteins; //number of proteins 
  int<lower=0> replicates; //number of replicates
  int<lower=0> tissues; //number of replicates
  matrix[cell_lines,proteins] prot_intensity[replicates,tissues]; // array of matricies of protein intensity from MS data per cell line (across replicates)
  }
  
parameters {
  vector<lower=0,upper=1>[cell_lines] enrichment[replicates,tissues];
  
  real<lower=0> sigma;
  real<lower=0> sigma_enrich;
  real<lower=0> mu_enrich;
  //vector<lower=0>[proteins] mu_chrom[tissues]; 
  //vector<lower=0>[proteins] mu_cyto[tissues]; 
  real<lower=0> mu_cell_line[tissues,cell_lines];
  real<lower=0> cyto_extra[tissues,cell_lines];
    
  // swapped dimensions to allow for easy vector addition (protein number is always column number of matrix)
  matrix<lower=0>[cell_lines,proteins] chrom_C[tissues]; 
  matrix<lower=0>[cell_lines,proteins] cyto_C[tissues];
  }

transformed parameters{
  matrix[cell_lines,proteins] prot_intensity_est[replicates,tissues];
   
  for(r in 1:replicates){
    for(c in 1:cell_lines){
      for(t in 1:tissues){
        prot_intensity_est[r,t,c,] = enrichment[r,t,c]*chrom_C[t,c,] + 
        (1-enrichment[r,t,c])*cyto_C[t,c,];
      }
    }
  }
}


model {
  
  // chromatin_fraction ~ beta(2.5, 2.5); removed for same issue with enrichment
  sigma ~ exponential(1); // generally should define priors for all variables
  sigma_enrich ~ exponential(1); // generally should define priors for all variables
  mu_enrich ~ beta(3,2);
  
  for(t in 1:tissues){
    mu_cell_line[t,] ~ normal(20,1);
    cyto_extra[t,] ~ normal(15,1);
    
    for(c in 1:cell_lines){ 
      
      //mu_cyto[t,c] ~ normal(mu_cell_line[t,c],1); 
      //mu_chrom[t,c] ~ normal(mu_cell_line[t,c] ,1); 
      chrom_C[t,c,] ~ normal(mu_cell_line[t,c], 1);
      cyto_C[t,c,] ~ normal(cyto_extra[t,c], 1); 
    }
  }

  for(r in 1:replicates){
    for(t in 1:tissues){
      enrichment[r,t] ~ normal(mu_enrich, sigma_enrich);

      for(p in 1:proteins){ 
        // generally discouraged to have complex maths in call to sample from distribution
        prot_intensity[r,t,,p] ~ normal(prot_intensity_est[r,t,,p], sigma);
      }
    }
  }
} 

generated quantities {
  real gen_prot_intensity[replicates, tissues, cell_lines, proteins];
  
  for(r in 1:replicates){
    for(p in 1:proteins){ 
      for(t in 1:tissues){
        // generally discouraged to have complex maths in call to sample from distribution
        gen_prot_intensity[r,t,,p] = normal_rng(prot_intensity_est[r,t,,p], sigma);
      }
    }
  }
}
