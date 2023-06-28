data {
  int<lower=0> cell_lines; //number of different cell_lines
  int<lower=0> proteins; //number of proteins 
  int<lower=0> replicates; //number of replicates
  matrix[cell_lines,proteins] prot_intensity[replicates]; // array of matricies of protein intensity from MS data per cell line (across replicates)
  }
  
parameters {
  vector<lower=0,upper=1>[cell_lines] enrichment[replicates];
  
  
  
  real<lower=0> sigma;
  real<lower=0> sigma_enrich;
  real<lower=0> mu_enrich;
  vector<lower=0>[proteins] mu_chrom; 
  vector<lower=0>[proteins] mu_cyto; 
  
    
  // swapped dimensions to allow for easy vector addition (protein number is always column number of matrix)
  matrix[cell_lines,proteins] chrom_C; 
  matrix[cell_lines,proteins] cyto_C;
  
 
 
  }

transformed parameters{
  matrix[cell_lines,proteins] prot_intensity_est[replicates];
   
  for(r in 1:replicates){
    for(c in 1:cell_lines){
      prot_intensity_est[r,c,] = enrichment[r,c]*chrom_C[c,]+
        (1-enrichment[r,c])*cyto_C[c,];
    }
  }
}


model {
  
  // chromatin_fraction ~ beta(2.5, 2.5); removed for same issue with enrichment
  sigma ~ exponential(1); // generally should define priors for all variables
  
  mu_cyto ~ normal(20,1) ; 
  mu_chrom ~ normal(15,1); 
    
  sigma_enrich ~ exponential(1); // generally should define priors for all variables
  mu_enrich ~ beta(3,2);
  
  for(C in 1:cell_lines){ 
    // remove nested for loop by calls cols of matrix
   
    chrom_C[C,] ~ normal(mu_chrom,1);
    cyto_C[C,] ~ normal(mu_cyto, 1); 
  }

  for(r in 1:replicates){
    enrichment[r] ~ normal(mu_enrich, sigma_enrich);

    for(p in 1:proteins){ 
      // generally discouraged to have complex maths in call to sample from distribution
      prot_intensity[r,,p] ~ normal(prot_intensity_est[r,,p], sigma);
    }
  }
} 

generated quantities {
  real gen_prot_intensity[replicates, cell_lines, proteins];
  
  for(r in 1:replicates){
    for(p in 1:proteins){ 
      // generally discouraged to have complex maths in call to sample from distribution
      gen_prot_intensity[r,,p] = normal_rng(prot_intensity_est[r,,p], sigma);
    }
  }
}
