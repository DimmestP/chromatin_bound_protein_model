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
  
  real<lower=0> sigma_total_intensity;
  
  vector<lower=0>[num_proteins] chrom_intensity_protein_base;
  
  vector<lower=0>[num_proteins] cyto_intensity_protein_base;
  
  vector[num_proteins] cyto_chrom_diff_protein[num_cell_lines - 1];

}

transformed parameters{
  vector[num_proteins] chrom_intensity_sum[num_cell_lines];
  vector[num_proteins] cyto_intensity_sum[num_cell_lines];
  
  for(c in 1:num_cell_lines){
    if(c == num_cell_lines){
      chrom_intensity_sum[c] = chrom_intensity_protein_base;
      cyto_intensity_sum[c] =  cyto_intensity_protein_base;
    }
    else{
      chrom_intensity_sum[c] = chrom_intensity_protein_base + cyto_chrom_diff_protein[c];
      cyto_intensity_sum[c] =  cyto_intensity_protein_base - cyto_chrom_diff_protein[c];
    }
  }
}

model {
  sigma_total_intensity ~ normal(1,0.5);
  
  for(c in 1:(num_cell_lines - 1)){
    cyto_chrom_diff_protein[c] ~ normal(0,3);
  }
  
  chrom_intensity_protein_base ~ normal(8,2);
  cyto_intensity_protein_base ~ normal(10,2);
  
  enrichment ~ normal(0.5,0.1);
  
  for (p in 1:num_proteins){
    for(s in 1:num_samples){
      target += log_sum_exp(log(enrichment[s]) + normal_lpdf( prot_intensity[s,p] | chrom_intensity_sum[cell_line[s],p], sigma_total_intensity), 
      log(1 - enrichment[s]) + normal_lpdf( prot_intensity[s,p] | cyto_intensity_sum[cell_line[s],p], sigma_total_intensity));
    }
  }
} 

generated quantities {
  real gen_prot_intensity[num_samples, num_proteins];
  
  for (p in 1:num_proteins){
    for(s in 1:num_samples){
      gen_prot_intensity[s,p] =  enrichment[s] * normal_rng(chrom_intensity_sum[cell_line[s],p], sigma_total_intensity) + 
       (1 - enrichment[s]) * normal_rng( cyto_intensity_sum[cell_line[s],p], sigma_total_intensity);
    }
  }
}
