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
  vector[num_samples - 1] enrichment_n_min_1;
  
  real<lower=0> sigma;
  
  vector[num_proteins] mu_chrom[num_tissue];
  
  vector<lower=0>[num_tissue] sigma_chrom;
  
  vector<lower=0>[num_proteins] chrom_intensity_protein;
    
  vector[num_proteins] chrom_intensity_cell_line_n_min_1[num_cell_lines - 1]; 
  }

transformed parameters{
  vector[num_proteins] prot_intensity_est[num_samples];
  
  vector[num_samples] enrichment;
  
  for(s in 1:(num_samples - 1) ){
    enrichment[s] = enrichment_n_min_1[s];
  }
  
  enrichment[num_samples] = 1 - sum(enrichment_n_min_1);
   
  for(s in 1:num_samples){
    
    if(cell_line[s] == num_cell_lines) prot_intensity_est[s] = enrichment[s] + chrom_intensity_protein;
    else prot_intensity_est[s] = enrichment[s] + chrom_intensity_cell_line_n_min_1[cell_line[s]] + chrom_intensity_protein;
  }
}


model {
  
  sigma ~ gamma(2,2);
  sigma_chrom ~ gamma(2,2);
  
  chrom_intensity_protein ~ normal(15,2);
  enrichment_n_min_1 ~ normal(0, 0.3);
  
  for(t in 1:num_tissue){
  mu_chrom[t] ~ normal(0,0.5); 
  }
  for(c in 1:(num_cell_lines-1)){
      
      chrom_intensity_cell_line_n_min_1[c] ~ normal(mu_chrom[tissue[c]], sigma_chrom[tissue[c]]);
  }

  for(s in 1:num_samples){
        // generally discouraged to have complex maths in call to sample from distribution
    prot_intensity[s] ~ lognormal(prot_intensity_est[s], sigma);
  }
} 

generated quantities {
  real gen_prot_intensity[num_samples, num_proteins];
  
  for(s in 1:num_samples){
    gen_prot_intensity[s] = lognormal_rng(prot_intensity_est[s], sigma);
  }
}
