data {
  int<lower=0> N; //number of samples
  int<lower=0> P; //number of proteins 
  matrix[N,P] prot_intensity; // matrix of protein intensity from MS data against samples, no missing values
  vector[N] input_enrichment; // each sample has been enriched for chromatin to a specific degree
  int cell_line[N];// type of cell line that each sample is (e.g. SKMEL-10), each cell line has 2 reps
  int<lower=0> cell_lines; //number of different cell_lines
  
  }
  
parameters {
    // my input_enrichment is my approximation of what percentage of my sample is pure chromatin and what % is pure non-chrom
    //so that chromatin proteins will be more abundant in samples that are higher enriched ie more % pure chromain
  
  //Im asking the model to convert my input enrichment values, to actual true percentages of pure chromatin in the sample
  real<lower=0>  b; // this forces a positive correlation in the sigmoidal transformation
  real  a;
  real  c; // constant for moving scaled enrichment along the x axis
  
  
  real<lower=0> sigma;
  vector<lower=0>[P] chrom; //  how much protein abundance each protein has, in pure chromatin samples e.g. global mean
  vector<lower=0>[P] cyto;// how much protein abundance each protein has, in pure non-chromatin samples
  // the chrom and cyto values characterise the general behavior of proteins, across samples
  // eg a chromatin protein would have high chrom value and low chrom value
  matrix[P,cell_lines] chrom_C; // here we allow the global chrom value of a protein to vary in a specific cell line, ideally in a hierachical manner
  matrix[P,cell_lines] cyto_C;// here we allow the global cyto value of a protein to vary in a specific cell line, ideally in a hierachical manner
 
  }

transformed parameters{
      vector[P] chromatin_fraction; // we take the ratio of chrom/(chrom+cyto) to determine if protein is mostly chrom or cytoplasmic
    for(i in 1:P){
    chromatin_fraction[i]= chrom[i]/(cyto[i]+chrom[i]);
    }
    
    vector[N] enrichment;// convert our input_enrichemnt values into actual % of pure chromatin in sample
    for(i in 1:N){
        enrichment[i] =c+ exp(a+input_enrichment[i]*b)/(1+exp(a+b*input_enrichment[i]));
    }
    

}


model {
 enrichment ~ beta(7.5, 7.5);// prior on enrichment of sample, since they are % of pure chrom, have to be bound between 0 and 1
 chromatin_fraction ~ beta(2.5, 2.5);// prior on chromatin_fraction of all proteins, since they are % protein on chrom, have to be bound between 0 and 1
   // enrichment ~ normal(input_enrichment*b+ a,0.01);
 chrom ~ normal(20,4) ;// prior of any protein abundance on pure chrom
 cyto ~ normal(20,4) ;// prior of any protein abundance on pure non-chrom
 
      for(C in 1:cell_lines){ // *trying to allow each protein abundance to vary in each cell line, from the global mean
     for(p in 1:P){
   
           chrom_C[p,C] ~ normal(chrom[P],1) ;
          cyto_C[p,C] ~ normal(cyto[P],1)  ;
          }
          
          }


     for(C in 1:cell_lines){ // building the matrix of protein abundance against samples,
    for(p in 1:P){// where each protein abundnace is an average of how much protein is on pure chrom of that cell line, how much of that sample is pure chrom (enrichment)

              for(n in 1:N){
   
            prot_intensity[n,p] ~ normal(enrichment[n]*(chrom_C[p,cell_line[n]])+
            (1-enrichment[n])*(cyto_C[p,cell_line[n]]), sigma); //multiplying enrichment of sample against the true chromatin and cytoplasm abundance of the protein
        }
    }
 }
 
    b ~ lognormal(0,0.5 );// prior for sigmoidal transformation
    a ~ normal(0,0.5);// prior for sigmoidal transformation

    
} 


