data {
  int<lower=1> num_samples; // total of different samples (cell lines * rep)
  int<lower=1> num_cell_lines; // tnumber of cell lines
  int<lower=1> num_proteins; //number of proteins 
  int<lower=1> cell_line[num_samples]; //  which cell line the sample comes from
  int<lower=1> num_tissue; // number of  tissues
  int<lower=1> tissue[num_cell_lines]; // which tissue the sample comes from
  vector<lower=0>[num_proteins] prot_intensity[num_samples]; // array of matricies of protein intensity from MS data per cell line (across replicates)
  }
  
parameters {
  vector<lower=0, upper=1>[num_samples] enrichment;
  
  vector<lower=0>[num_proteins] prot_intensity_est_cyto;
  
  vector<lower=0>[num_proteins] prot_intensity_est_chrom[num_samples];
  
  real<lower=0> sigma;
  
  vector<lower=0>[num_proteins] mu_chrom[num_tissue];
  
  vector<lower=0>[num_tissue] sigma_chrom;
  
  real<lower=0> sigma_chrom_est;
  
  vector<lower=0>[num_proteins] chrom_intensity_protein;
    
  vector[num_proteins] chrom_intensity_cell_line[num_cell_lines]; 
  }

transformed parameters{
  vector[num_proteins] prot_intensity_est[num_samples];
  vector[num_proteins] prot_intensity_sum_chrom[num_samples];
   
  for(s in 1:num_samples){
      prot_intensity_sum_chrom[s] =  chrom_intensity_cell_line[cell_line[s]] + chrom_intensity_protein;
      
      prot_intensity_est[s] = enrichment[s] * prot_intensity_est_chrom[s] + (1 - enrichment[s]) * prot_intensity_est_cyto;
  }
}


model {
  
  sigma ~ gamma(2,2);
  sigma_chrom ~ gamma(2,2);
  sigma_chrom_est ~ gamma(2,2);
  
  chrom_intensity_protein ~ normal(15,2);
  prot_intensity_est_cyto ~ lognormal(13,2);
  enrichment ~ beta(5, 2);
  
  for(t in 1:num_tissue){
  mu_chrom[t] ~ normal(0,0.3); 
  }
  for(c in 1:num_cell_lines){
      
      chrom_intensity_cell_line[c] ~ normal(mu_chrom[tissue[c]], sigma_chrom[tissue[c]]);
  }

  for(s in 1:num_samples){
    prot_intensity_est_chrom[s] ~ lognormal(prot_intensity_sum_chrom[s], sigma_chrom_est);
        // generally discouraged to have complex maths in call to sample from distribution
    prot_intensity[s] ~ normal(prot_intensity_est[s], sigma);
  }
} 

generated quantities {
  real gen_prot_intensity[num_samples, num_proteins];
  
  for(s in 1:num_samples){
    gen_prot_intensity[s] = normal_rng(prot_intensity_est[s], sigma);
  }
}
