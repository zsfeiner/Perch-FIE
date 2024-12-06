//Data block
data {
  int<lower=1> N;   //Sample size
  int<lower=1> K;   //Number of continuous parameters
  vector[N] Age;   //Age, explanatory variable
  vector[N] TL;    //Length, response variable
  vector[N] RW;   //relative weight
  vector[N] GDD5; //annual GDD5
  vector[N] TP; //annual mean TP
  vector[N] Fished; //Fishing indicator
  int<lower=1> nCohorts;   //Number of cohorts
  int<lower=1,upper=nCohorts> Cohort[N];  //Cohort ID
  int<lower=0,upper=1> Mat[N];   // Maturation status
  real mean_TL;    //mean total length for scaling
  real<lower=0> sd_TL;  //sd total length for scaling
  real mean_RW;   //mean relative weight for scaling
  real<lower=0> sd_RW;  //sd relative weith for scaling
  real mean_GDD5;   //mean GDD5 for scaling
  real<lower=0> sd_GDD5;  //sd GDD5 for scaling
  real mean_TP;   //mean TP for scaling
  real<lower=0> sd_TP;  //sd TP for scaling
  int<lower=1> nAbiotics; //number of years with abiotic vars
  matrix[nAbiotics,5] abiotics; //table of Year, CohortYear, GDD5Annual, TP, Fishing
  int<lower=1> Year[N]; //Year sampled
//  int<lower=0> nFishing; //Cohorts fishing
//  matrix[nFishing, 3] fishing; //table of cohortyear, cohort number, and fishing indicator (1 before 1997, 0 after)

}


transformed data {
  vector[N] sc_TL;
  vector[N] sc_RW;
  vector[N] sc_GDD5;
  vector[N] sc_TP;
  matrix[nAbiotics,5] sc_abiotics;
  sc_TL = (TL - mean_TL)/sd_TL;
  sc_RW = (RW - mean_RW)/sd_RW;
  sc_GDD5 = (GDD5 - mean_GDD5)/sd_GDD5;
  sc_TP = (TP - mean_TP)/sd_TP;
  sc_abiotics[,1] = abiotics[,1]; //Year 1983-2018
  sc_abiotics[,2] = abiotics[,2]; //CohortYear 5-40
  sc_abiotics[,3] = (abiotics[,3] - mean_GDD5) / sd_GDD5; //GDD scaled
  sc_abiotics[,4] = (abiotics[,4] - mean_TP) / sd_TP; //TP scaled
  sc_abiotics[,5] = abiotics[,5]; //Fishing indicator linked to year
}


//Parameter block
parameters {
  //Maturation parameters
  vector[K+3] beta;  //Mean coefficient vector of length K+2 (intercept + age + TL + age:TL + age:RW)
  matrix[K+2,nCohorts] z_u;  //Random effects (intercept, age effect, TL effect, interaction effect)
  vector<lower=0>[K+2] sigma_u;  //Variance for random effects (intercept, age effect, TL effect, interaction effect)
  cholesky_factor_corr[K+2] L_u;
  
  //Growth parameters
  vector[nCohorts] phi;  //Coefficient vector of length nCohorts - intercept
  vector[nCohorts] gamma;  //Coefficient vector of length nCohorts - slope
  real phi_mu;
  real<lower=0> phi_sigma;
  real gamma_mu;
  real<lower=0> gamma_sigma;
  real<lower=0> tau;    //Error, lower limit zero
  
  //Condition parameters
  vector[nCohorts] omega1;  //Coefficient vector of length nCohorts - intercept
  vector[nCohorts] omega2;  //Coefficient vector of length nCohorts - age effect
  vector[nCohorts] omega3;  //Coefficient vector of length nCohorts - length effect
  vector[nCohorts] omega4;  //Coefficient vector of length nCohorts - age*length effect
  real <lower=0> zeta;  //Error, lower limit zero
  
}

transformed parameters {
  matrix[K+2,nCohorts] u;
  u = diag_pre_multiply(sigma_u, L_u) * z_u;  //Cohort random effects
}

//Model block
//Maturation model
model {
  real mu[N];  //Maturation
  vector[N] y;  //Growth
  vector[N] w;  //Condition
  
  //Priors
  //Maturation
  beta ~ normal(0,10);
  L_u ~ lkj_corr_cholesky(2.0);
  to_vector(z_u) ~ normal(0,1);
  sigma_u ~ cauchy(0,5);
  
  //Growth
  phi ~ normal(phi_mu, phi_sigma);
  phi_mu ~ normal(100, 20);
  phi_sigma ~ cauchy(0,5);
  
  gamma ~ normal(gamma_mu, gamma_sigma);
  gamma_mu ~ normal(100,20);
  gamma_sigma ~ cauchy(0,5);
  
  //Condition
  omega1 ~ normal(0,1); 
  omega2 ~ normal(0,1); 
  omega3 ~ normal(0,1); 
  omega4 ~ normal(0,1); 
  zeta ~ cauchy(0,5);
  
  tau ~ cauchy(0,5);
  
  for (i in 1:N) {
    //Maturation
    mu[i] = (beta[1] + u[1,Cohort[i]]) + 
    (beta[2] + u[2,Cohort[i]]) * Age[i] +
    (beta[3] + u[3,Cohort[i]]) * sc_TL[i] +
    (beta[4] + u[4,Cohort[i]]) * sc_RW[i] +    // relative weight add elsewhere
    (beta[5] + u[5,Cohort[i]]) * Age[i] * sc_TL[i] + 
    (beta[6] + u[6,Cohort[i]]) * Age[i] * sc_RW[i] + //new line add elsewhere
    (beta[7] + u[7,Cohort[i]]) * sc_GDD5[i] +  //add annual GDD5 - need to have interaction with age?  I think no?;  
    (beta[8] + u[8,Cohort[i]]) * sc_TP[i] +  //add TP - need to have interaction with age?  I think no? - should these have cohort random effects?
    beta[9] * Fished[i];  //Fishing effect - no cohort random slopes with fishing
 
    //Growth
    y[i] = phi[Cohort[i]] + gamma[Cohort[i]] * log(Age[i]);
    
    //Condition
    w[i] = omega1[Cohort[i]] + omega2[Cohort[i]]*Age[i] + omega3[Cohort[i]]*TL[i] + omega4[Cohort[i]]*Age[i]*TL[i];
  
  }
  
  Mat ~ bernoulli_logit(mu);  //Maturation
  TL ~ normal(y, tau);  //Growth
  RW ~ normal(w, zeta); //Condition
  

}


generated quantities {
  vector[N] s;
  vector[N] inc;
  //vector[N] prev_wr;
  vector[N] wr_inc;
  vector[N] prev_sc_wr;
  vector<lower=0,upper=1>[N] p;
  vector<lower=0,upper=1>[N] prev_p;
  vector[N] m;

  for (i in 1:N) {
    //if(Cohort[i] == 1){ //if Cohort year is 1983 (first cohort year) then skip, we don't have TP for 1982
    //    s[i] = 0;
    //    p[i] = 0;
    //    prev_p[i] = 0;
    //    m[i] = 0;
    //}
    //else{
      
      if(Age[i] == 1) {
        inc[i] = 0;
        s[i] = (0 - mean_TL) / sd_TL;
        wr_inc[i] = 0;
        //prev_wr[i] = 0;
        prev_sc_wr[i] = (0 - mean_RW) / sd_RW;
      }
      else {
        inc[i] = (phi[Cohort[i]] + gamma[Cohort[i]] * log(Age[i])) - (phi[Cohort[i]] + gamma[Cohort[i]] * log(Age[i]-1));
        s[i] = ((TL[i] - inc[i]) - mean_TL) / sd_TL;
        
        wr_inc[i] = (omega1[Cohort[i]] + omega2[Cohort[i]]*Age[i] + omega3[Cohort[i]]*TL[i] + omega4[Cohort[i]]*Age[i]*TL[i]) - (omega1[Cohort[i]] + omega2[Cohort[i]]*(Age[i]-1) + omega3[Cohort[i]]*(TL[i]-inc[i]) + omega4[Cohort[i]]*(Age[i]-1)*(TL[i]-inc[i]));
        //prev_wr[i] = omega1[Cohort[i]] + omega2[Cohort[i]]*(Age[i]-1) + omega3[Cohort[i]]*(TL[i]-inc[i]) + omega4[Cohort[i]]*(Age[i]-1)*(TL[i]-inc[i]);
        prev_sc_wr[i] = ((RW[i] - wr_inc[i]) - mean_RW) / sd_RW;
      }
      
      p[i] = inv_logit((beta[1] + u[1,Cohort[i]]) + 
      (beta[2] + u[2,Cohort[i]]) * Age[i] +
      (beta[3] + u[3,Cohort[i]]) * sc_TL[i] +
      (beta[4] + u[4,Cohort[i]]) * sc_RW[i] + 
      (beta[5] + u[5,Cohort[i]]) * Age[i] * sc_TL[i] +
      (beta[6] + u[6,Cohort[i]]) * Age[i] * sc_RW[i] +
      (beta[7] + u[7,Cohort[i]]) * sc_GDD5[i] +
      (beta[8] + u[8,Cohort[i]]) * sc_TP[i] +
      beta[9] * Fished[i]);
      
      prev_p[i] = inv_logit((beta[1] + u[1,Cohort[i]]) + 
      (beta[2] + u[2,Cohort[i]]) * (Age[i] - 1) +
      (beta[3] + u[3,Cohort[i]]) * s[i] +
      (beta[4] + u[4,Cohort[i]]) * prev_sc_wr[i] +               //Use modeled Wr from previous year and length
      (beta[5] + u[5,Cohort[i]]) * (Age[i] - 1) * s[i] +
      (beta[6] + u[6,Cohort[i]]) * (Age[i] - 1) * prev_sc_wr[i] +
      (beta[7] + u[7,Cohort[i]]) * sc_abiotics[Year[i],3] +    // GDD from previous year
      (beta[8] + u[8,Cohort[i]]) * sc_abiotics[Year[i],4] +    // TP from previous year
      beta[9] * sc_abiotics[Year[i],5]);    // Fishing from previous year
      
      m[i] = (p[i] - prev_p[i]) / (1 - prev_p[i] + 0.00001); //add small offset, prevents division by 0
    //}
  } 
}
