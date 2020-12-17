require(MASS)
require(MBESS)
require(rstan)
require(matrixStats)

####################################################################
#            MODEL PARAMETER GENERATION
#
# These parameters are reported in the paper
# and taken from the real data analysis
####################################################################

# Item parameters

   beta  <- rnorm(25,4.19,0.38)
   alpha <- rnorm(25,1.42,0.37)

# Tau for control group

	cor0 <- matrix(c(1,.84,.84,1),2,2)
	tau0 <- mvrnorm(33,c(0,0),cor2cov(cor0,c(0.18,0.17)))

# Tau for experimental group (Items Disclosed)

	cor1 <- matrix(c(1,.56,.56,1),2,2)
	tau1 <- mvrnorm(30,c(0.05,1.34),cor2cov(cor1,c(0.55,0.47)))

# Tau for experimental group (Items and Answers Disclosed)

	cor2 <- matrix(c(1,.43,.43,1),2,2)
	tau2 <- mvrnorm(30,c(-0.18,1.47),cor2cov(cor2,c(0.43,0.46)))

# A vector for item status (0: not disclosed, 1:disclosed)

	C    <- c(0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0)

####################################################################
#
#         RESPONSE TIME GENERATION
#
####################################################################

# Note that the response time data generated is
# already on the log scale
	
# Control Group
	
  rt0 <- matrix(nrow=33,ncol=25)

  for(i in 1:33){
    for(j in 1:25){
      p_t = beta[j]-tau0[i,1]
      p_c = beta[j]-tau0[i,2]
      p   = p_t*(1-C[j]) + p_c*C[j]
      rt0[i,j] = rnorm(1,p,1/alpha[j])   
    }
  }

# Experimental Group 1 (Item Disclosed)
	
	rt1 <- matrix(nrow=30,ncol=25)
	
	for(i in 1:30){
	  for(j in 1:25){
	    p_t = beta[j]-tau1[i,1]
	    p_c = beta[j]-tau1[i,2]
	    p   = p_t*(1-C[j]) + p_c*C[j]
	    rt1[i,j] = rnorm(1,p,1/alpha[j])   
	  }
	}

# Experimental Group 2 (Item and Answers Disclosed)
	
	rt2 <- matrix(nrow=30,ncol=25)
	
	for(i in 1:30){
	  for(j in 1:25){
	    p_t = beta[j]-tau2[i,1]
	    p_c = beta[j]-tau2[i,2]
	    p   = p_t*(1-C[j]) + p_c*C[j]
	    rt2[i,j] = rnorm(1,p,1/alpha[j])   
	  }
	}

# Combine the groups 

rt <- rbind(cbind(data.frame(exp(rt0)),gr=1),
            cbind(data.frame(exp(rt1)),gr=2),
            cbind(data.frame(exp(rt1)),gr=3))

#describeBy(rt[,1:25],rt$gr)

####################################################################
#
#         FIT THE MODEL
#
####################################################################

# Reshaping data for each group into long format
# and prepare data objects for Stan modeling

  rt0 <- as.data.frame(rt0)
  rt0$id <- 1:33
  rt0.long <- reshape(data        = rt0,
                      idvar       = "id",
                      varying     = list(colnames(rt0)[1:25]),
                      timevar     = "Item",
                      times       = 1:25,
                      v.names      = "RT",
                      direction   = "long")

  rt1 <- as.data.frame(rt1)
  rt1$id <- 1:30
  rt1.long <- reshape(data        = rt1,
                      idvar       = "id",
                      varying     = list(colnames(rt1)[1:25]),
                      timevar     = "Item",
                      times       = 1:25,
                      v.names      = "RT",
                      direction   = "long")
  
  rt2 <- as.data.frame(rt2)
  rt2$id <- 1:30
  rt2.long <- reshape(data        = rt2,
                      idvar       = "id",
                      varying     = list(colnames(rt2)[1:25]),
                      timevar     = "Item",
                      times       = 1:25,
                      v.names      = "RT",
                      direction   = "long")
  
  
  rt0.long$i.status  <- ifelse((rt0.long$Item)%%2 == 1,0,1)
  rt1.long$i.status  <- ifelse((rt1.long$Item)%%2 == 1,0,1)
  rt2.long$i.status  <- ifelse((rt2.long$Item)%%2 == 1,0,1)

  id0            <- rt0.long$id
  item0          <- rt0.long$Item
  istatus0       <- rt0.long$i.status
  rt0            <- rt0.long$RT
  n_obs0         <- length(rt0.long$RT)

  id1            <- rt1.long$id
  item1          <- rt1.long$Item
  istatus1       <- rt1.long$i.status
  rt1            <- rt1.long$RT
  n_obs1         <- length(rt1.long$RT)
  
  id2            <- rt2.long$id
  item2          <- rt2.long$Item
  istatus2       <- rt2.long$i.status
  rt2            <- rt2.long$RT
  n_obs2         <- length(rt2.long$RT)
  
  
# Stan Model Syntax

code_rt<-'

data{

	int  J;               //  number of items

	int  i1;              //  number of individuals in control group
	int  n_obs1;          //  number of observed responses in control group
	int  istatus1[n_obs1];//  item indicator of a response in control group
  int  item1[n_obs1];   //  item id in control group
  int  id1[n_obs1];     //  person id in control group
  real Y1[n_obs1];      //  vector of log response times in control group

	int  i2;              //  number of individuals in experimental group1
	int  n_obs2;          //  number of observed responses in experimental group 1
	int  istatus2[n_obs2];//  item indicator of a response in experimental group 1
  int  item2[n_obs2];   //  item id in experimental group 1
  int  id2[n_obs2];     //  person id in experimental group 1
  real Y2[n_obs2];      //  vector of log response times in experimental group 1

	int  i3;              //  number of individuals in experimental group 2
	int  n_obs3;          //  number of observed responses in experimental group 2 
	int  istatus3[n_obs3];//  item indicator of a response in experimental group 2
  int  item3[n_obs3];   //  item id in experimental group 2
  int  id3[n_obs3];     //  person id in experimental group 2
  real Y3[n_obs3];      //  vector of log response times in experimental group 2

  vector[2] mu_tau1;    // vector for mean latent speed ability 
                        // for compromised and uncompromised items 
                        // in control group
}

parameters {

	vector[J] beta;               //item specific time intensity parameter
	  real mu_beta;               // mean of beta
	  real<lower=0> sigma_beta;   // sd of beta

	vector <lower=0> [J]  alpha;  // item specific residual standard deviation,
  	real mu_alpha;              // mean alpha   
  	real<lower=0> sigma_alpha;  // sd alpha

	vector[2] tau1[i1];  // latent working speed for compromised and uncompromised items 
	                     // in control group
	                     
	vector[2] tau2[i2];  // latent working speed for compromised and uncompromised items 
	                     // in experimental group 1

	vector[2] tau3[i3];  // latent working speed for compromised and uncompromised items 
	                     // in experimental group 2

  vector[2] mu_tau2;   // vector for mean latent speed ability 
                       // for compromised and uncompromised items 
                       // in experimental group 1
  
  vector[2] mu_tau3;   // vector for mean latent speed ability 
                       // for compromised and uncompromised items 
                       // in experimental group 2

	vector<lower=0>[2] sigma1;// std.dev. of latent speed abilities in control group
	vector<lower=0>[2] sigma2;// std.dev. of latent speed abilities in experimental group 1
	vector<lower=0>[2] sigma3;// std.dev. of latent speed abilities in experimental group 2

	corr_matrix[2] omega1; // correlation matrix between latent speed parameters
	                       // in control group 1
	                       
	corr_matrix[2] omega2; // correlation matrix between latent speed parameters
	                       // in experimental group 1
	
	corr_matrix[2] omega3; // correlation matrix between latent speed parameters
	                       // in experimental group 2
	                     
}


model{

	sigma1 ~ cauchy(0,2.5);
	omega1 ~ lkj_corr(1);
	tau1   ~ multi_normal(mu_tau1,quad_form_diag(omega1, sigma1)); 
	      
	      // http://stla.github.io/stlapblog/posts/StanLKJprior.html
	      // http://modernstatisticalworkflow.blogspot.com/2016/09/what-are-covariance-matrices-anyway.html

	mu_tau2 ~ normal(0,5);
	sigma2  ~ cauchy(0,2.5);
	omega2  ~ lkj_corr(1);
	tau2    ~ multi_normal(mu_tau2,quad_form_diag(omega2, sigma2));

	mu_tau3 ~ normal(0,5);
	sigma3  ~ cauchy(0,2.5);
	omega3  ~ lkj_corr(1);
	tau3    ~ multi_normal(mu_tau3,quad_form_diag(omega3, sigma3));

	mu_beta  ~ normal(0,5);
	  sigma_beta ~ cauchy(0,2.5);
	  beta       ~ normal(mu_beta,sigma_beta);

	mu_alpha     ~ normal(0,5);
	  sigma_alpha ~ cauchy(0,2.5);
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

	for (i in 1:n_obs3) {        
		real p_t = beta[item3[i]]-tau3[id3[i],1];
		real p_c = beta[item3[i]]-tau3[id3[i],2];
		real p2 = (1-istatus3[i])*p_t +(istatus3[i])*p_c;
		Y3[i] ~ normal(p2,1/(alpha[item3[i]]));
	}
}

'

# Estimation

data_rt<-list(J        = 25,
              i1       = 33,
              n_obs1   = n_obs0,
              item1    = item0,
              id1      = id0,
              istatus1 = istatus0,
              Y1       = rt0,
              i2       = 30,
              n_obs2   = n_obs1,
              item2    = item1,
              id2      = id1,
              istatus2 = istatus1,
              Y2       = rt1,
              i3       = 30,
              n_obs3   = n_obs2,
              item3    = item2,
              id3      = id2,
              istatus3 = istatus2,
              Y3       = rt2,
              mu_tau1  = c(0,0))


rt_stan <- stan(model_code = code_rt, 
                data       = data_rt, 
                iter       = 10000,
                chains     = 4,
                control    = list(max_treedepth = 15,adapt_delta=0.99))



print(rt_stan, 
      pars = c("beta","mu_beta","sigma_beta",
               "alpha","mu_alpha","sigma_alpha"),
      probs = c(0.025,0.975), 
      digits = 3)

print(rt_stan, 
      pars = c("mu_tau2","mu_tau3",
               "sigma1","sigma2","sigma3",
               "omega1","omega2","omega3"),
      probs = c(0.025,0.975), 
      digits = 3)

save.image(paste0("C:\\Users\\Dr Zopluoglu\\Desktop\\simdata\\rep",xxx,".Rdata"))
rm(list=ls(all=TRUE))


#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################


other        <- as.data.frame(matrix(nrow=17,ncol=100))
rownames(other)<- c("mu_beta","sigma_beta","mu_alpha","sigma_alpha",
                    "mu_tau2_1","mu_tau2_2","mu_tau3_1","mu_tau3_2",
                    "sigma1_1","sigma1_2","sigma2_1","sigma2_2","sigma3_1","sigma3_2",
                    "omega1","omega2","omega3")

Rhat       <- as.data.frame(matrix(nrow=17,ncol=100))
rownames(Rhat)<- c("mu_beta","sigma_beta","mu_alpha","sigma_alpha",
                    "mu_tau2_1","mu_tau2_2","mu_tau3_1","mu_tau3_2",
                    "sigma1_1","sigma1_2","sigma2_1","sigma2_2","sigma3_1","sigma3_2",
                    "omega1","omega2","omega3")

Lower       <- as.data.frame(matrix(nrow=17,ncol=100))
rownames(Lower)<- c("mu_beta","sigma_beta","mu_alpha","sigma_alpha",
                    "mu_tau2_1","mu_tau2_2","mu_tau3_1","mu_tau3_2",
                    "sigma1_1","sigma1_2","sigma2_1","sigma2_2","sigma3_1","sigma3_2",
                    "omega1","omega2","omega3")

Higher       <- as.data.frame(matrix(nrow=17,ncol=100))
rownames(Higher)<- c("mu_beta","sigma_beta","mu_alpha","sigma_alpha",
                    "mu_tau2_1","mu_tau2_2","mu_tau3_1","mu_tau3_2",
                    "sigma1_1","sigma1_2","sigma2_1","sigma2_2","sigma3_1","sigma3_2",
                    "omega1","omega2","omega3")

write.table(other,"Parameter Estimates.csv",sep=",")
write.table(Rhat,"Rhat.csv",sep=",")
write.table(Lower,"LowerBoundary.csv",sep=",")
write.table(Higher,"UpperBoundary.csv",sep=",")

for(x in 1:100){

	load(paste0("rep",x,".Rdata"))
	
	other[,x] <- as.numeric(summary(rt_stan)$summary[c(26,27,53,54,241:250,252,256,260),1])
	Rhat[,x] <- as.numeric(summary(rt_stan)$summary[c(26,27,53,54,241:250,252,256,260),10])
	Lower[,x] <- as.numeric(summary(rt_stan)$summary[c(26,27,53,54,241:250,252,256,260),4])
	Higher[,x] <- as.numeric(summary(rt_stan)$summary[c(26,27,53,54,241:250,252,256,260),8])

	print(x)
	rm(list=setdiff(ls(),c("Higher","Lower","Rhat","other"))) 
}



par.summary <- as.data.frame(matrix(
					c(4.19,0.38,1.18,0.16,0.05,1.34,-0.18,1.47,0.18,0.17,0.55,0.47,0.43,0.46,0.84,0.56,0.43),
					nrow=17,ncol=1))
rownames(par.summary)<- c("mu_beta","sigma_beta","mu_alpha","sigma_alpha",
                    "mu_tau2_1","mu_tau2_2","mu_tau3_1","mu_tau3_2",
                    "sigma1_1","sigma1_2","sigma2_1","sigma2_2","sigma3_1","sigma3_2",
                    "omega1","omega2","omega3")


par.summary$Rhat <- rowMeans(Rhat < matrix(1.1,17,100))

t = matrix(c(4.19,0.38,1.18,0.16,0.05,1.34,-0.18,1.47,0.18,0.17,0.55,0.47,0.43,0.46,0.84,0.56,0.43),17,100,byrow=F)
par.summary$Coverage <- rowMeans(t < Higher & t > Lower)

par.summary$Mean <- rowMeans(other)

par.summary$SDs <- rowSds(as.matrix(other))

round(matrixStats::rowRanges(as.matrix(Rhat)),3)





