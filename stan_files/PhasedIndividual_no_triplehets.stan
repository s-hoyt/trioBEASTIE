functions {
real likelihood(int[,] count,int[] isPhased,int numSites,real p)
{
   real loglik=0;
   for(site in 1:numSites) {
      int mat=count[site,1];
      int pat=count[site,2];
      
      if(isPhased[site]) 
         loglik+=binomial_lpmf(mat | mat+pat,p);
      else {
         loglik += 0;
      //   real phase1=log(0.5) + binomial_lpmf(mat | mat+pat,p);
      //   real phase2=log(0.5) + binomial_lpmf(pat | mat+pat,p);
      //   loglik+=log_sum_exp(phase1,phase2);
      }
   }
   return loglik;
}}
data {
   int N_SITES;
   int<lower=0> count[N_SITES,2]; // [site,individual,haplotype]
   int<lower=0,upper=1> isPhased[N_SITES]; // triple hets are unphased
}
parameters 
{
   real<lower=0.000001,upper=1000000> theta; // amount of ASE
}
transformed parameters 
{
   real p = theta / (1 + theta);
}
model 
{
   // Priors:
   log2(theta) ~ normal(0, 1);
   target += -log(theta * log(2)); // Jacobian

   // Likelihoods:
   target+=likelihood(count,isPhased,N_SITES,p);
}


