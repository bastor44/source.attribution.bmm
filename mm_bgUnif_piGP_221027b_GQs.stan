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
  //int M[D]; //Number of basis functions in each dimension
  array[D] int M; 
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
  //int<lower=0,upper=N> coordinates[N, D]; // indexing for ages of sources and recipients
  array[N,D] int<lower=0, upper=N> coordinates; 
  int<lower=0> Nj_max; //Max obs per recip 1yr age
  int<lower=0> Nj_max_coarse; //Max obs per recip age band
  //int<lower=0> age_recip_coarse[Nj_max_coarse, A, A];
  array[Nj_max_coarse, A, A] int<lower=0> age_recip_coarse; 
  //int<lower=1> pt_idx[N]; // unique recipient ID for pair ij
  array[N] int<lower=1> pt_idx; 
  //int pt_map[P,N]; // idx of ID for pair ij
  array[P,N] int pt_map; 
  matrix[S,N] idx_src; // idx of source categories
  matrix[S,N] idx_rec; // idx of recipient categories


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

generated quantities
{
  vector[N] log_lik;

  for (i in 1:N){
      log_lik[i] = log_mix(y_mix[i], gamma_lpdf(y[i] | tpair_alpha[i], tpair_beta[i]), unif_lpdf);
  }

}
