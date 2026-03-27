//====================================================================
// TrioBEAST.stan
// (C)2022 W.H. Majoros (bmajoros@alumni.duke.edu)
// This is OPEN SOURCE software, released under the GPL 3.0 license.
//
// This version has been refactored to construct the likelihood
// terms automatically for each requested mode of inheritance.
//
// Indexing of arrays:
//   Individuals: 1=mother, 2=father, 3=child
//   Haplotypes:  1=maternal, 2=paternal
//
// For child, genotype XY means X is maternal and Y is paternal
// For a parent, gentype XY means X is transmitted to child, Y is not
//====================================================================

functions {

real binom_lpmf(int x,int n,int het,real p) 
{
   return het ? binomial_lpmf(x|n,p) : 0.0;
}

real computeElem(int[,,] count,int[,] het,int[] isPhased,int site,
   int MI,int FI,int CI,real MP,real FP,real CP) 
{
   int MN=count[site,1,1]+count[site,1,2]; // Mother
   int FN=count[site,2,1]+count[site,2,2]; // Father
   int CN=count[site,3,1]+count[site,3,2]; // Child
   int MC=count[site,1,MI];
   int FC=count[site,2,FI];
   int CC=count[site,3,CI];
   if(isPhased[site]) 
      return 
           binom_lpmf(MC | MN,het[site,1],MP) // Mother
         + binom_lpmf(FC | FN,het[site,2],FP) // Father
         + binom_lpmf(CC | CN,het[site,3],CP);// Child
   else {
      return 0;
   //   real phase1=log(0.5) +
   //        binom_lpmf(MC | MN,het[site,1],MP) // Mother
   //      + binom_lpmf(FC | FN,het[site,2],FP) // Father
   //      + binom_lpmf(CC | CN,het[site,3],CP);// Child
   //   real phase2=log(0.5) +
   //        binom_lpmf(MC | MN,het[site,1],1-MP) // Mother
   //      + binom_lpmf(FC | FN,het[site,2],1-FP) // Father
   //      + binom_lpmf(CC | CN,het[site,3],1-CP);// Child
   //   return log_sum_exp(phase1,phase2);
   }  
}

int countDenovos(int parent,int[,] V) {
   // Count as de novo only if parent has no affected copies and child does
   if(V[parent,1]==0 && V[parent,2]==0 && V[3,parent]==1) return 1;
   return 0;
}

int isHet(int parent,int[,] V) { return V[parent,1]!=V[parent,2]; }

real recombTerm(int parent,int[,] V,real logRecomb,real logNoRecomb) {
   // There is no recombination term if the parent is homozygous, because
   // the terms for recombination and no recombination cancel
   //Edit Stephanie Oct 2025 - 'homozygous' here doesn't refer to genotype, but to mode affected status
   real s=0;
   if(isHet(parent,V))
      if(V[3,parent]==V[parent,1]) s+=logNoRecomb;
      else s+=logRecomb;
   return s;
}

int affectedCopy(int indiv,int[,] V) {
   if(V[indiv,1]==1) return 1;
   if(V[indiv,2]==1) return 2;
   return 1;
}

real getP(int indiv,int[,] V,real p) {
   if(V[indiv,1]==V[indiv,2]) return 0.5; //no affected copy
   return p; //one affected copy
}

real modeLik(int site,int[,] mode,int[,,] count,int[,] het,
   real p,int[] isPhased) 
{
   // Count some things and initialize variables
   real s=0; 
   
   int mCopy=affectedCopy(1,mode); int fCopy=affectedCopy(2,mode);
   int cCopy=affectedCopy(3,mode); real mP=getP(1,mode,p);
   real fP=getP(2,mode,p); real cP=getP(3,mode,p);
   
   // Compute binomial terms
   s+=computeElem(count,het,isPhased,site,mCopy,fCopy,cCopy,mP,fP,cP);
   return s;
}

real likelihoods(int[,] mode,int[,,] count,int[,] het,
   real logAffected,real logUnaffected,real logDenovo,real logNoDenovo,
   real logRecomb,real logNoRecomb,real p,int N_SITES,int[] isPhased) 
{
   real array = 0.0;
      
   // Priors for recombinations (parents) and de novos (child)
   int numDenovos=countDenovos(1,mode)+countDenovos(2,mode);
   real recomb=recombTerm(1,mode,logRecomb,logNoRecomb)+
       recombTerm(2,mode,logRecomb,logNoRecomb);
   array+=recomb;
   array+=numDenovos*logDenovo + (2-numDenovos)*logNoDenovo;

   // Priors for parent affected status
   int numParentAffected=mode[1,1]+mode[1,2]+mode[2,1]+mode[2,2];
   int numParentUnaffected=4-numParentAffected;
   array+=numParentAffected*logAffected+numParentUnaffected*logUnaffected;

   for(site in 1:N_SITES) {
      array += modeLik(site,mode,count,het,p,isPhased);
   }
   return array;
}
}

data {
   int N_SITES; // number of sites (variants)
   int mode[3,2]; // [individual,haplotype]
   int<lower=0,upper=1> het[N_SITES,3]; // [site,indvidual]
   int<lower=0> count[N_SITES,3,2]; // [site,individual,haplotype]
   int<lower=0,upper=1> isPhased[N_SITES]; // triple hets are unphased
   real<lower=0,upper=1> probAffected; // prior prob of 1 parent copy affected
}

transformed data {
   real logAffected=log(probAffected);
   real probRecomb = 0.001;
   real probDenovo = 0.0001;
}

parameters 
{
   real<lower=0.001,upper=1000> theta; // amount of ASE
}

transformed parameters 
{
   real p = theta / (1 + theta);
   real logUnaffected=log(1-probAffected);
   real logDenovo=log(probDenovo);
   real logNoDenovo=log(1-probDenovo);
   real logRecomb=log(probRecomb);
   real logNoRecomb=log(1-probRecomb);
}

model 
{
   // Priors:
   log2(theta) ~ normal(0, 1);
   target += -log(theta * log(2)); // Jacobian
   
   // Likelihoods:
   target+=likelihoods(mode,count,het,logAffected,logUnaffected,
      logDenovo,logNoDenovo,logRecomb,logNoRecomb,p,N_SITES,isPhased);
}

generated quantities 
{
   real numerator = likelihoods(mode,count,het,logAffected,logUnaffected,
      logDenovo,logNoDenovo,logRecomb,logNoRecomb,p,N_SITES,isPhased);
}
