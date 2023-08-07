data {
  int<lower=0> cell_lines; //number of different cell_lines
  int<lower=0> proteins; //number of proteins 
  int<lower=0> replicates; //number of replicates
  int<lower=0> tissues; //number of replicates
  matrix<lower=0>[cell_lines,proteins] prot_intensity[replicates,tissues]; // array of matricies of protein intensity from MS data per cell line (across replicates)
  }
  
parameters { 
  vector<lower=0>[cell_lines] enrichment[replicates,tissues];
  
  real<lower=0> sigma;
  
  vector[proteins] mu_chrom[tissues];
  
  vector<lower=0>[tissues] sigma_chrom;
  
  vector<lower=0>[proteins] chrom_intensity_protein;
    
  matrix[cell_lines,proteins] chrom_intensity_cell_line[tissues]; 
  }

transformed parameters{
  matrix[cell_lines,proteins] prot_intensity_est[replicates,tissues];
   
  for(r in 1:replicates){
    for(c in 1:cell_lines){
      for(t in 1:tissues){
        for(p in 1:proteins){
        prot_intensity_est[r,t,c,p] = enrichment[r,t,c] * (chrom_intensity_cell_line[t,c,p] + chrom_intensity_protein[p]);
        }
      }
    }
  }
}


model {
  
  sigma ~ gamma(2,2);
  sigma_chrom ~ gamma(2,2);
  
  chrom_intensity_protein ~ normal(15,2);
  
  for(t in 1:tissues){
  mu_chrom[t] ~ normal(0,4); 
    
    for(p in 1:proteins){
      
      chrom_intensity_cell_line[t,,p] ~ normal(mu_chrom[t,p], sigma_chrom[t]);
    }
  }

  for(r in 1:replicates){
    for(t in 1:tissues){
      enrichment[r,t] ~ normal(1, 0.3);

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
