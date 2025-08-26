functions{
  vector diagSPD_EQ(real alpha, real rho, real L, int M)
    {
      return alpha * sqrt(sqrt(2*pi()) * rho) * exp(-0.25*(rho*pi()/2/L)^2 * linspaced_vector(M, 1, M)^2);
    }

  // Eigenfunctions
  matrix PHI(int M, real L, vector x)
  {
    return sin(diag_post_multiply(rep_matrix(pi()/(2*L) * (x+L), M), linspaced_vector(M, 1, M)))/sqrt(L);
  }

  matrix kron_mvprod(matrix A, matrix B, matrix V)
  {
    return (B*V) * transpose(A);
  }

  // Transform restructured matrix back to a squared matrix
  matrix inverse_restruct(matrix B, array[] int nn_idx)
  {
    return to_matrix(to_vector(B)[nn_idx], cols(B), cols(B), 0);
  }

matrix hsgp(int A, real alpha, real rho1, real rho2, real L1, real L2, int M1, int M2,
            matrix PHI1, matrix PHI2, matrix z)
{
  vector[M1] sqrt_spd_1 = diagSPD_EQ(alpha, rho1, L1, M1);
  vector[M2] sqrt_spd_2 = diagSPD_EQ(alpha, rho2, L2, M2);

  matrix[A,A] f = kron_mvprod(
    diag_post_multiply( PHI1, sqrt_spd_1 ),
    diag_post_multiply( PHI2, sqrt_spd_2 ),
    z
  );

  return(f);
}

}

data{
  int<lower=1> N; //No. observations or no.pairs
  int<lower=1> P; //number of unique recipients
  int<lower=1> S; //number of source categories
  vector<lower=0>[N] y; //Genetic Distances
  vector<lower=0>[N] x; //Time elapsed

  int<lower=1> D; //Number of Dimensions
  vector[D] L; //Boundaries for each dimension
  array[D] int M; //Number of basis functions in each dimension
  int<lower=1> M_nD; //M1*M2

  matrix[M_nD, D] indices; //Matrix of D-tuples
  int<lower=0> A; // number of age groups
  int<lower=0> n; // number of 1yr ages

  // parameters from clock model
  real log_alpha1;
  real log_alpha1_pair_sd;
  real log_phi;
  real log_phi_pair_sd;

  //Indexing
  int<lower=0> Nk_max; //Max obs per age
  array[N, D] int<lower=0, upper=N> coordinates; // indexing for ages of sources and recipients
  int<lower=0> Nj_max; //Max obs per recip 1yr age
  int<lower=0> Nj_max_coarse; //Max obs per recip age band
  array[Nj_max_coarse, A, A] int<lower=0> age_recip_coarse; 
  array[N] int<lower=1> pt_idx; // unique recipient ID for pair ij
  array[P, N] int pt_map; // idx of ID for pair ij
  matrix[S,N] idx_src; // idx of source categories
  matrix[S,N] idx_rec; // idx of recipient categories
  //vector<lower=0>[N] idx_true_pairs; // indicator for true pairs


  int<lower=1> Nu_pairs; //number of unique observed time elapsed
  array[N] int<lower=1, upper=Nu_pairs> IDX_UNIQUE_PAIR; // index of unique values
  matrix[n, D] ages; //Age matrix (Source-Recip Columns)
  real<lower=0> sd1;
  real<lower=0> sd2;
}

transformed data{
  real<lower=0> unif_pdf = 1/max(y);
  real unif_lpdf = log(unif_pdf);
  matrix[n, M[1]] PHI1;
  matrix[n, M[2]] PHI2;
  PHI1 = PHI(M[1], L[1], ages[,1]);
  PHI2 = PHI(M[2], L[2], ages[,2]);
}

parameters{
  vector[N] log_alpha1_pair;
  vector[N] log_phi_pair;
  real logit_y_mix_0; //Mixing probability parameter (on logit scale)

  vector<lower=0>[D] lscale; //Lengthscale
  real<lower=0> gpscale; //Amplitude
  matrix[M[1],M[2]] z1;
}

transformed parameters{
  matrix[n,n] f; //Latent function
  {
        f = hsgp(n, gpscale, lscale[1], lscale[2],
              L[1], L[2], M[1], M[2], PHI1, PHI2, z1);
  }
  vector<lower=0, upper=1>[N] y_mix; //Mixing probability parameter (on natural scale)
  vector[N] tpair_mean; //Mean of transmission pair
  vector[N] tpair_alpha; //Shape param of signal component
  vector[N] tpair_beta; //Rate param of signal component

  for(i in 1:N){
      y_mix[i] = inv_logit(logit_y_mix_0 + f[coordinates[i,1],coordinates[i,2]]);
  }
  tpair_beta = exp(-(log_phi + log_phi_pair)); //Beta = e^(-phi0_log + phipair_log) where phipair is a vector i.e log(beta) = -log(phi0) + log(phipair)
  tpair_mean = exp(log_alpha1 + log_alpha1_pair) .* x; //mu = e^(log_alpha1 + log_alpha_pair) so log(mu) = log_alpha1 + log_alpha_pair
  tpair_alpha = tpair_mean .* tpair_beta; //For gamma dist, mu = alpha / beta -> alpha = mu .* beta
}

model{
  target += normal_lpdf(log_alpha1_pair | 0, log_alpha1_pair_sd);
  target += normal_lpdf(log_phi_pair | 0, log_phi_pair_sd);
  target += normal_lpdf(logit_y_mix_0 | 0, 2); //mean = 0.005 = 200/200^2, (was 0 for 0.5) -5.293305, 5
  // GP priors
  lscale[1] ~ inv_gamma(5, 5);
  lscale[2] ~ inv_gamma(5, 5);
  gpscale ~ normal(0,0.15);

 for(i in 1:M[1]){
    for(j in 1:M[2]){
      z1[i,j] ~ normal(0,1);
    }
  }

  for (i in 1:N){
    target += log_mix(y_mix[i], gamma_lpdf(y[i] | tpair_alpha[i], tpair_beta[i]), unif_lpdf); //Background pairs given uniform
  }
}

generated quantities{
  vector[N] d_pred;
  int npairs;
  int k;
  real num;     // create vector of appropriate size for all npairs
  vector[2] den_one;
  vector<lower=0,upper=1>[N] tpair_prob;
  vector<lower=0,upper=1>[N] tpair_prob_w;
  row_vector<lower=0,upper=1>[N] idx;
  matrix<lower=0>[S,S] flows; //rows=sources
  matrix<lower=0>[S,S] pflows; // prop. total flows between groups
  matrix<lower=0>[S,S] pflows_to; // prop. flows within recipients B from each source
  vector<lower=0>[S] pflows_from; // prop. flows from sources A
  //vector<lower=0>[S] true_flows; // prop. flows from sources A
  //vector<lower=0>[S]  AE_from; // AE for flows_from sources A
  //real<lower=0> MAE_from; // MAE for flows_from sources A

  for (i in 1:N)
  {

    // generate distances using parameters for posterior predictive check
    d_pred[i] = log_sum_exp(log(y_mix[i]) + gamma_rng(tpair_alpha[i], tpair_beta[i]), log1m(y_mix[i]) + uniform_rng(0,max(y)));

    npairs = sum(pt_map[pt_idx[i],]); //count number of pairs for recipient in pair i

    vector[npairs+1] den;     // create vector of appropriate size for all npairs

    // numerator of pr(transm pair)
    num = gamma_lpdf(y[i] | tpair_alpha[i], tpair_beta[i]) +
    log(y_mix[i]); // first add the llk of pair i being a true pair

    k = 1;
    if(npairs>1){
      for(j in 1:N)
      { // add llk of other pairs for same recipient being non-pairs
          if(i!=j && pt_map[pt_idx[i],j]==1){
            num += unif_lpdf + log1m(y_mix[i]);
            k += 1;
          }
       }
    }

    // denominator of pr(transm pair)
    k = 1;
    if(npairs>1){
      for(j in 1:N)
      { // all combinations of pairs with other non-pairs
          if(pt_map[pt_idx[i],j]==1){
              den[k] = gamma_lpdf(y[j] | tpair_alpha[j], tpair_beta[j]) + log(y_mix[i]);
              den[npairs+1] = unif_lpdf + log1m(y_mix[i]);
            for(l in 1:N)
            {
              if(j!=l && pt_map[pt_idx[i],l]==1){
                den[k] += unif_lpdf + log1m(y_mix[i]);
                den[npairs+1] += unif_lpdf + log1m(y_mix[i]);
              }
           }
            k += 1;
          }
       }
    }else{ // if only one possible pair for the recipient, only add pr(pair | non-pair)
      den_one[1] = gamma_lpdf(y[i] | tpair_alpha[i], tpair_beta[i]) + log(y_mix[i]);
      den_one[2] = unif_lpdf + log1m(y_mix[i]);
    }

    if(npairs==1){
        tpair_prob_w[i] = exp(num-log_sum_exp(den_one));
    }else{
        tpair_prob_w[i] = exp(num-log_sum_exp(den));
    }

    tpair_prob[i] = exp(
          gamma_lpdf(y[i] | tpair_alpha[i], tpair_beta[i]) +
          log(y_mix[i]) -
          (
            log_mix( y_mix[i],
              gamma_lpdf(y[i] | tpair_alpha[i], tpair_beta[i]),
              unif_lpdf
              )
          )
        );

  }
  for(a in 1:S){
    for(b in 1:S){
        idx = (idx_src[a,] .* idx_rec[b,]);
        flows[a,b] = idx * tpair_prob_w;
      }
  }
    pflows = flows/sum(flows[:,:]);
  for(a in 1:S){
    for(b in 1:S){
        pflows_to[a,b] = flows[a,b]/sum(flows[:,b]);
    }
  }
  for(a in 1:S){
    pflows_from[a] = sum(flows[a,:])/sum(flows);
    //true_flows[a] = (idx_src[a,]*idx_true_pairs)/sum(idx_src*idx_true_pairs);
    //AE_from[a] = abs(pflows_from[a]-true_flows[a]);
  }
  //MAE_from = sum(AE_from)/S;
}
