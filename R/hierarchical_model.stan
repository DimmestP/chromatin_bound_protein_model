data {
  int<lower=0> cell_lines; //number of different cell_lines
  int<lower=0> proteins; //number of proteins 
  int<lower=0> replicates; //number of replicates
  int<lower=0> tissues; //number of replicates
  matrix<lower=0>[cell_lines,proteins] prot_intensity[replicates,tissues]; // array of matricies of protein intensity from MS data per cell line (across replicates)
  }
  
parameters {
  vector<lower=0,upper=1>[cell_lines] enrichment[replicates,tissues];
  
  real<lower=0> sigma;
  real<lower=0> sigma_enrich;
  real<lower=0,upper=1> mu_enrich;
  
  vector<lower=0>[proteins] mu_chrom[tissues];
  vector<lower=0>[proteins] mu_cyto[tissues];
    
  // swapped dimensions to allow for easy vector addition (protein number is always column number of matrix)
  matrix<lower=0>[cell_lines,proteins] chrom_intensity[tissues]; 
  matrix<lower=0>[cell_lines,proteins] cyto_intensity[tissues];
  }

transformed parameters{
  matrix[cell_lines,proteins] prot_intensity_est[replicates,tissues];
   
  for(r in 1:replicates){
    for(c in 1:cell_lines){
      for(t in 1:tissues){
        prot_intensity_est[r,t,c,] = enrichment[r,t,c]*chrom_intensity[t,c,] + 
        (1-enrichment[r,t,c])*cyto_intensity[t,c,];
      }
    }
  }
}


model {
  
  // chromatin_fraction ~ beta(2.5, 2.5); removed for same issue with enrichment
  sigma ~ exponential(1); // generally should define priors for all variables
  sigma_enrich ~ exponential(1); // generally should define priors for all variables
  mu_enrich ~ beta(4,2);
  
  for(t in 1:tissues){
  mu_chrom[t] ~ normal(20,2); 
  mu_cyto[t] ~ normal(18,1); 
    
    for(c in 1:cell_lines){
      
      chrom_intensity[t,c,] ~ normal(mu_chrom[t], 1);
      cyto_intensity[t,c,] ~ normal(mu_cyto[t], 1); 
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
