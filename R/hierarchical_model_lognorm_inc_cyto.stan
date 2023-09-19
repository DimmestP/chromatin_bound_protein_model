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
  vector<lower=0,upper=1>[num_samples] enrichment;
  
  vector<lower=0>[num_proteins] sigma_total_intensity;
  
  vector<lower=0>[num_proteins] chrom_intensity_protein;
  
  vector<lower=0>[num_proteins] cyto_intensity_protein;
    
  vector<lower=0>[num_proteins] cyto_chrom_diff_cell_line_n_min_1[num_cell_lines - 1]; 
  }

transformed parameters{
  vector[num_proteins] prot_intensity_est[num_samples];
  vector[num_proteins] cyto_intensity_temp;
  vector[num_proteins] chrom_intensity_temp;
   
  for(s in 1:num_samples){
    if(cell_line[s] == num_cell_lines) {
      cyto_intensity_temp = cyto_intensity_protein;
      chrom_intensity_temp = chrom_intensity_protein;
      }
    else {
      cyto_intensity_temp = cyto_intensity_protein - cyto_chrom_diff_cell_line_n_min_1[cell_line[s]];
      chrom_intensity_temp = chrom_intensity_protein + cyto_chrom_diff_cell_line_n_min_1[cell_line[s]];
      }
    prot_intensity_est[s] = enrichment[s] * chrom_intensity_temp + (1 - enrichment[s]) * cyto_intensity_temp;
  }
}


model {
  
  sigma_total_intensity ~ lognormal(15,1);
  
  cyto_intensity_protein ~ lognormal(17,2);
  chrom_intensity_protein ~ lognormal(15,2);
  
  enrichment ~ normal(0.7,0.2);
  
  for(c in 1:(num_cell_lines-1)){
      
      cyto_chrom_diff_cell_line_n_min_1[c] ~ lognormal(15, 3);
  }

  for(s in 1:num_samples){
        // generally discouraged to have complex maths in call to sample from distribution
    prot_intensity[s] ~ normal(prot_intensity_est[s], sigma_total_intensity);
  }
} 

generated quantities {
  real gen_prot_intensity[num_samples, num_proteins];
  
  for(s in 1:num_samples){
    gen_prot_intensity[s] = normal_rng(prot_intensity_est[s], sigma_total_intensity);
  }
}
