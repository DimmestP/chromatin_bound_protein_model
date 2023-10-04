data {
  int<lower=1> num_samples; // total of different samples (cell lines * rep)
  int<lower=1> num_cell_lines; // tnumber of cell lines
  int<lower=1> num_proteins; //number of proteins 
  int<lower=0> N_mis; // total number of missing values
  int<lower=1> cell_line[num_samples]; //  which cell line the sample comes from
  int<lower=1> num_tissue; // number of  tissues
  int<lower=1> tissue[num_cell_lines]; // which tissue the sample comes from
  vector<lower=0>[num_proteins] prot_intensity_obs[num_samples]; // array of matricies of protein intensity from MS data per cell line (across replicates)
  
  int<lower=1, upper=num_proteins> ii_mis[N_mis]; // protein number that is missing a value accross all samples
  int<lower=0, upper=num_proteins> N_mis_sample[num_samples]; // number of missing values per sample
  }
  
parameters {
  
  vector<lower=0>[N_mis] prot_intensity_mis;
  
  vector<lower=0,upper=1>[num_samples] enrichment;
  
  real<lower=0> sigma_total_intensity;
  
  vector<lower=0>[num_proteins] chrom_intensity_protein_base[num_tissue];
  
  vector[num_proteins] total_intensity_protein_cell_line[num_cell_lines - 1];
  
  vector<lower=0>[num_proteins] cyto_intensity_protein_base[num_tissue];

}

transformed parameters {
  vector[num_proteins] chrom_intensity_sum[num_cell_lines];
  
  for(c in 1:num_cell_lines){
    int cur_tissue = tissue[c];
    if(c == num_cell_lines){
      chrom_intensity_sum[c] = chrom_intensity_protein_base[cur_tissue];
    }
    else{
      chrom_intensity_sum[c] = chrom_intensity_protein_base[cur_tissue] + total_intensity_protein_cell_line[c];
      }
  }
}

model {
  // Define Temporary variables
  int pos_start;
  vector[num_proteins] temp_prot_intensity_comb;
  vector[num_proteins] prot_intensity_comb[num_samples];
  
  prot_intensity_comb=prot_intensity_obs;
  
  // combine inferred missing values and known values
  pos_start = 1;
  
  for(s in 1:num_samples){
    if (N_mis_sample[s] > 0){ 
      int cur_cell_line = cell_line[s];
      int cur_tissue = tissue[cur_cell_line];
      for(p in 1:N_mis_sample[s]){
        prot_intensity_mis[pos_start+p-1] ~ normal(log_sum_exp(chrom_intensity_protein_base[cur_tissue,ii_mis[(pos_start+p-1)]],cyto_intensity_protein_base[cur_tissue,ii_mis[(pos_start+p-1)]]),0.2);
      }
      temp_prot_intensity_comb = prot_intensity_comb[s,];
      
      temp_prot_intensity_comb[ii_mis[pos_start:(pos_start+N_mis_sample[s]-1)]] = prot_intensity_mis[pos_start:(pos_start+N_mis_sample[s]-1)];
      
      prot_intensity_comb[s,] = temp_prot_intensity_comb;
      
      pos_start = pos_start + N_mis_sample[s];
    }
  }
  
  sigma_total_intensity ~ normal(1,0.5);
  
  for(c in 1:(num_cell_lines - 1)){
    total_intensity_protein_cell_line[c] ~ normal(0,0.5);
  }
  
  for(t in 1:num_tissue){
    cyto_intensity_protein_base[t] ~ normal(17,3);
    chrom_intensity_protein_base[t] ~ normal(15,3);
  }
  
  enrichment ~ normal(0.6,0.1);
  
  for (p in 1:num_proteins){
    for(s in 1:num_samples){
      int cur_cell_line = cell_line[s];
      target += log_sum_exp(log(enrichment[s]) + normal_lpdf( prot_intensity_comb[s,p] | chrom_intensity_sum[cell_line[s],p], sigma_total_intensity), 
      log(1 - enrichment[s]) + normal_lpdf( prot_intensity_comb[s,p] | cyto_intensity_protein_base[tissue[cur_cell_line],p], sigma_total_intensity));
    }
  }
} 

generated quantities {
  real gen_prot_intensity[num_samples, num_proteins];
  
  
  for (p in 1:num_proteins){
    for(s in 1:num_samples){
      int cur_cell_line = cell_line[s];
      gen_prot_intensity[s,p] =  enrichment[s] * normal_rng(chrom_intensity_sum[cell_line[s],p], sigma_total_intensity) + 
       (1 - enrichment[s]) * normal_rng( cyto_intensity_protein_base[tissue[cur_cell_line],p], sigma_total_intensity);
    }
  }
}
