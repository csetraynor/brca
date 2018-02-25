functions {
  vector sqrt_vec(vector x) {
    vector[dims(x)[1]] res;

    for (m in 1:dims(x)[1]){
      res[m] = sqrt(x[m]);
    }

    return res;
  }
//horseshoe prior
  vector hs_prior_lp(real r1_global, real r2_global, vector r1_local, vector r2_local, real nu) {
    r1_global ~ normal(0.0, 1.0);
    r2_global ~ inv_gamma(0.5, 0.5);

    r1_local ~ normal(0.0, 1.0);
    r2_local ~ inv_gamma(0.5 * nu, 0.5 * nu);

    return (r1_global * sqrt(r2_global)) * r1_local .* sqrt_vec(r2_local);
  }
}

data {
    int<lower=2> K;
    int<lower=0> N;
    int<lower=1> D;
    int<lower=1,upper=K> y[N];
    matrix[N,K] x; 
  }
  parameters {
    vector beta;
  }
  model {
     vector[N] beta_x; 
     beta_x = beta_raw * (5 * x); 
     for (k in 1:K) 
       beta_raw[k] ~ normal(0,1); 
    
      y[n] ~ categorical(softmax(beta * x[n]));
  }
  

















// 
data {
  int<lower=1> N; //number of samples (nrow dataset)
  int<lower=0> Y; // number of metastasis
  int<lower = 0> M; //number of mutations
  int<lower=1, upper=N> s[N];     // sample id
  real y[Y]; // metastasis
  matrix[N, M] x;                 // mutation matrix
  real<lower=0> sigma[J]; // s.e. of effect estimates 
}
parameters {
  real mu; 
  real<lower=0> tau;
  real eta[J];
}
transformed parameters {
  real theta[J];
  for (j in 1:J)
    theta[j] = mu + tau * eta[j];
}
model {
  target += normal_lpdf(eta | 0, 1);
  target += normal_lpdf(y | theta, sigma);
}