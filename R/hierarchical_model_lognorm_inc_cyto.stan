data {
  int<lower=1> num_samples; // total of different samples (cell lines * rep)
  int<lower=1> num_cell_lines; // tnumber of cell lines
  int<lower=1> num_proteins; //number of proteins
  int<lower=1> cell_line[num_samples]; //  which cell line the sample comes from
  int<lower=1> num_tissue; // number of  tissues
  int<lower=1> tissue[num_cell_lines]; // which tissue the sample comes from
  vector<lower=0>[num_proteins] prot_intensity_obs[num_samples]; // array of matricies of protein intensity from MS data per cell line (across replicates)
  }
  
parameters {
  
  vector<lower=0,upper=1>[num_samples] enrichment;

  real<lower=0,upper=1> on_off_chrom_mix[num_proteins];
  
  real<lower=0> sigma_total_intensity;
  
  vector<lower=0>[num_proteins] total_intensity_protein_base[num_tissue];
  
  vector<lower=0,upper=1>[num_proteins] cyto_chrom_ratio [num_cell_lines];
}

transformed parameters {
  vector[num_proteins] chrom_intensity_sum[num_cell_lines];
  vector[num_proteins] cyto_intensity_sum[num_cell_lines];
  
  
  
  for(c in 1:num_cell_lines){
    int cur_tissue = tissue[c];
    chrom_intensity_sum[c] = log(cyto_chrom_ratio[c]) + total_intensity_protein_base[cur_tissue];
    cyto_intensity_sum[c] = log(1 - cyto_chrom_ratio[c]) + total_intensity_protein_base[cur_tissue];
  }
}

model {
  sigma_total_intensity ~ normal(1,0.5);
  
  for(c in 1:num_cell_lines){
    for(p in 1:num_proteins){
       target += log_sum_exp( log(on_off_chrom_mix[p]) + normal_lpdf(cyto_chrom_ratio[c,p] | 0.7,0.1), 
       log(1 - on_off_chrom_mix[p]) + normal_lpdf(cyto_chrom_ratio[c,p] | 0.3,0.1));
    }
  }
  
  for(t in 1:num_tissue){
    total_intensity_protein_base[t] ~ normal(17,2);
  }
  
  enrichment ~ normal(0.6,0.1);
  
  for (p in 1:num_proteins){
    for(s in 1:num_samples){
      int cur_cell_line = cell_line[s];
      target += log_sum_exp(log(enrichment[s]) + normal_lpdf( prot_intensity_obs[s,p] |  chrom_intensity_sum[cur_cell_line,p], sigma_total_intensity), 
      log(1 - enrichment[s]) + normal_lpdf( prot_intensity_obs[s,p] | cyto_intensity_sum[cur_cell_line,p], sigma_total_intensity));
    }
  }
} 

generated quantities {
  real gen_prot_intensity[num_samples, num_proteins];
  
  
  for (p in 1:num_proteins){
    for(s in 1:num_samples){
      int cur_cell_line = cell_line[s];
      gen_prot_intensity[s,p] =  enrichment[s] * normal_rng(chrom_intensity_sum[cur_cell_line,p], sigma_total_intensity) + 
       (1 - enrichment[s]) * normal_rng( cyto_intensity_sum[cur_cell_line,p], sigma_total_intensity);
    }
  }
}
