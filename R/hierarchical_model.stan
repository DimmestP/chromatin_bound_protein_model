data {
  int<lower=0> cell_lines; //number of different cell_lines
  int<lower=0> proteins; //number of proteins 
  int<lower=0> replicates; //number of replicates
  matrix[cell_lines,proteins] prot_intensity[replicates]; // array of matricies of protein intensity from MS data per cell line (across replicates)
  vector[cell_lines] input_enrichment[replicates]; // array of vectors of chromatin enrichment
  }
  
parameters {
  real<lower=0>  b; 
  real  a;
  vector[cell_lines] noiseless_input_enrichment[replicates];
  vector<lower=0,upper=1>[cell_lines] enrichment[replicates];
  
  
  
  real<lower=0> sigma;
  vector<lower=0>[proteins] chrom; 
  vector<lower=0>[proteins] cyto;
  
    
  // swapped dimensions to allow for easy vector addition (protein number is always column number of matrix)
  matrix[cell_lines,proteins] chrom_C; 
  matrix[cell_lines,proteins] cyto_C;
  
 
 
  }

transformed parameters{
  vector[proteins] chromatin_fraction; 
  matrix[cell_lines,proteins] prot_intensity_est[replicates];
  vector[cell_lines] enrichment_est[replicates];
   
  chromatin_fraction = chrom ./ (cyto+chrom); // as chrom and cyto are vectors you can do vector addition and division to remove for loop
  
  for(r in 1:replicates) {
    enrichment_est[r] = inv_logit(a + b * noiseless_input_enrichment[r]); // the way to use stan's inbuilt sigmoid function
  }
  for(r in 1:replicates){
    for(c in 1:cell_lines){
      prot_intensity_est[r,c,] = enrichment[r,c]*chrom_C[c,]+
        (1-enrichment[r,c])*cyto_C[c,];
    }
  }
}


model {
  a ~ normal(0,1);
  b ~ normal(2,1);
  // chromatin_fraction ~ beta(2.5, 2.5); removed for same issue with enrichment
  sigma ~ cauchy(1,1); // generally should define priors for all variables
  chrom ~ normal(20,4) ; 
  cyto ~ normal(20,4) ; 
  
  for(C in 1:cell_lines){ 
    // remove nested for loop by calls cols of matrix
   
    chrom_C[C,] ~ normal(chrom,1);
    cyto_C[C,] ~ normal(cyto,1); 
  }

  for(r in 1:replicates){
    noiseless_input_enrichment[r] ~ normal(input_enrichment[r], 1); // all input data should be assumed to have noise!
    enrichment[r] ~ normal(enrichment_est[r], 0.1); 

    for(p in 1:proteins){ 
      // generally discouraged to have complex maths in call to sample from distribution
      prot_intensity[r,,p] ~ normal(prot_intensity_est[r,,p], sigma);
    }
  }
} 
