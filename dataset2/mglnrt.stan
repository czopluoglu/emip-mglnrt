data{

	int  J;                      //  number of items

	int  i1;                      //  number of individuals in condition1
	int  n_obs1;                  //  number of observed responses in condition 1
	int  istatus1[n_obs1];        //  item indicator of a response in condition 1
  int  item1[n_obs1];           //  item id in condition 1
  int  id1[n_obs1];             //  person id in condition 1
  real Y1[n_obs1];               //  vector of log of response times in condition 1

	int  i2;                      //  number of individuals in condition 2
	int  n_obs2;                  //  number of observed responses in condition 2
	int  istatus2[n_obs2];        //  item indicator of a response in condition 2
  int  item2[n_obs2];           //  item id in condition 2
  int  id2[n_obs2];             //  person id in condition 2
  real Y2[n_obs2];              //  vector of log of response times in condition 2

  vector[2] mu_tau1;            // vector of average latent speed in reference group
}

parameters {

	vector[J] beta;               //item specific time intensity parameter
	real mu_beta;
	real<lower=0> sigma_beta;

	vector <lower=0> [J]  alpha;  // item specific residual standard deviation,
	real mu_alpha;
	real<lower=0> sigma_alpha;

	vector[2] tau1[i1];         // person specific latent working speed for uncompromised items Condition 1

	vector[2] tau2[i2];         // person specific latent working speed for uncompromised items Condition 1

  vector[2] mu_tau2;          // vector of average latent speed in the second group

	corr_matrix[2] omega1;
	corr_matrix[2] omega2;

	vector<lower=0>[2] sigma1;
	vector<lower=0>[2] sigma2;
}


model{

	sigma1   ~ cauchy(0,2.5);
	omega1   ~ lkj_corr(1);

	tau1    ~ multi_normal(mu_tau1,quad_form_diag(omega1, sigma1));
                                                                       

	mu_tau2  ~ normal(0,5);
	sigma2   ~ cauchy(0,2.5);
	omega2   ~ lkj_corr(1);
	tau2    ~ multi_normal(mu_tau2,quad_form_diag(omega2, sigma2));

	mu_beta      ~ normal(0,5);
	sigma_beta   ~ cauchy(0,2.5);
	beta     ~ normal(mu_beta,sigma_beta);

	mu_alpha     ~ normal(0,5);
	sigma_alpha  ~ cauchy(0,2.5);
	alpha       ~ normal(mu_alpha,sigma_alpha);
                  

	for (i in 1:n_obs1) {    

		real p_t = beta[item1[i]]-tau1[id1[i],1];
		real p_c = beta[item1[i]]-tau1[id1[i],2];

		real p2 = (1-istatus1[i])*p_t +(istatus1[i])*p_c;

		Y1[i] ~ normal(p2,1/(alpha[item1[i]]));
  }

	for (i in 1:n_obs2) {        

		real p_t = beta[item2[i]]-tau2[id2[i],1];
		real p_c = beta[item2[i]]-tau2[id2[i],2];

		real p2 = (1-istatus2[i])*p_t +(istatus2[i])*p_c;

		Y2[i] ~ normal(p2,1/(alpha[item2[i]]));
  }
}

generated quantities{

  real rt1_rep[n_obs1];
  real rt2_rep[n_obs2];

	for (i in 1:n_obs1) {    

		real p_t = beta[item1[i]]-tau1[id1[i],1];
		real p_c = beta[item1[i]]-tau1[id1[i],1];

		real p2 = (1-istatus1[i])*p_t +(istatus1[i])*p_c;

		rt1_rep[i] = normal_rng(p2,1/(alpha[item1[i]]));
  }

  for (i in 1:n_obs2) {        

		real p_t = beta[item2[i]]-tau2[id2[i],1];
		real p_c = beta[item2[i]]-tau2[id2[i],2];

		real p2 = (1-istatus2[i])*p_t +(istatus2[i])*p_c;

		rt2_rep[i] = normal_rng(p2,1/(alpha[item2[i]]));
  }
}
