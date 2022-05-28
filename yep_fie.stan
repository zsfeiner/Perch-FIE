//Data block
data {
  int<lower=1> N;   //Sample size
  int<lower=1> K;   //Number of continuous parameters
  vector[N] Age;   //Age, explanatory variable
  vector[N] TL;    //Length, response variable
  int<lower=1> nCohorts;   //Number of cohorts
  int<lower=1,upper=nCohorts> Cohort[N];  //Cohort ID
  int<lower=0,upper=1> Mat[N];   // Maturation status
  real mean_TL;    //mean total length for scaling
  real<lower=0> sd_TL;  //sd total length for scaling
}


transformed data {
  vector[N] sc_TL;
  sc_TL = (TL - mean_TL)/sd_TL;
}


//Parameter block
parameters {
  //Maturation parameters
  vector[K+2] beta;  //Mean coefficient vector of length K+2 (intercept + age + TL + age:TL)
  matrix[4,nCohorts] z_u;  //Random effects (intercept, age effect, TL effect, interaction effect)
  vector<lower=0>[4] sigma_u;  //Variance for random effects (intercept, age effect, TL effect, interaction effect)
  cholesky_factor_corr[4] L_u;
  
  //Growth parameters
  vector[nCohorts] phi;  //Coefficient vector of length nCohorts - intercept
  vector[nCohorts] gamma;  //Coefficient vector of length nCohorts - slope
  real phi_mu;
  real<lower=0> phi_sigma;
  real gamma_mu;
  real<lower=0> gamma_sigma;
  real<lower=0> sigma;    //Error, lower limit zero
}

transformed parameters {
  matrix[4,nCohorts] u;
  u = diag_pre_multiply(sigma_u, L_u) * z_u;  //Cohort random effects
}

//Model block
//Maturation model
model {
  real mu[N];  //Maturation
  vector[N] y;  //Growth
  
  //Priors
  //Maturation
  beta ~ normal(0,10);
  L_u ~ lkj_corr_cholesky(2.0);
  to_vector(z_u) ~ normal(0,1);
  
  //Growth
  phi ~ normal(phi_mu, phi_sigma);
  phi_mu ~ normal(200, 20);
  phi_sigma ~ cauchy(0,5);
  
  gamma ~ normal(gamma_mu, gamma_sigma);
  gamma_mu ~ normal(0,10);
  gamma_sigma ~ cauchy(0,5);
  
  sigma ~ cauchy(0,5);
  
  for (i in 1:N) {
    //Maturation
    mu[i] = (beta[1] + u[1,Cohort[i]]) + 
    (beta[2] + u[2,Cohort[i]]) * Age[i] +
    (beta[3] + u[3,Cohort[i]]) * sc_TL[i] +
    (beta[4] + u[4,Cohort[i]]) * Age[i] * sc_TL[i];
    
    //Growth
    y[i] = phi[Cohort[i]] + gamma[Cohort[i]] * Age[i];
  
  }
  
  Mat ~ bernoulli_logit(mu);  //Maturation
  TL ~ normal(y, sigma);  //Growth
  

}


generated quantities {
  vector[N] s;
  vector<lower=0,upper=1>[N] p;
  vector<lower=0,upper=1>[N] prev_p;
  vector[N] m;

  for (i in 1:N) {
    s[i] = ((TL[i] - gamma[Cohort[i]]) - mean_TL) / sd_TL;
    
    p[i] = inv_logit((beta[1] + u[1,Cohort[i]]) + 
    (beta[2] + u[2,Cohort[i]]) * Age[i] +
    (beta[3] + u[3,Cohort[i]]) * sc_TL[i] +
    (beta[4] + u[4,Cohort[i]]) * Age[i] * sc_TL[i]);
    
    prev_p[i] = inv_logit((beta[1] + u[1,Cohort[i]]) + 
    (beta[2] + u[2,Cohort[i]]) * (Age[i] - 1) +
    (beta[3] + u[3,Cohort[i]]) * s[i] +
    (beta[4] + u[4,Cohort[i]]) * (Age[i] - 1) * s[i]);
    
    m[i] = (p[i] - prev_p[i]) / (1 - prev_p[i]);
  
  } 
}
