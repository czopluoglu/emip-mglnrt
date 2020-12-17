
require(rstan)
options(mc.cores = parallel::detectCores())
Sys.setenv(LOCAL_CPPFLAGS = '-march=native')

require(psych)
require(here)

###########################################################################################################################

d <- read.csv(here("data/Real Data/UVA Experiment Data.csv"))

###########################################################################################################################

d.sub <- d[,c("ID","COND","Q1RT","Q2RT","Q3RT","Q4RT","Q5RT","Q6RT","Q7RT","Q8RT","Q9RT","Q10RT",
              "Q11RT","Q12RT","Q13RT","Q14RT","Q15RT","Q16RT","Q17RT","Q18RT","Q19RT","Q20RT",
              "Q21RT","Q22RT","Q23RT","Q24RT","Q25RT")]



subset <- log(d.sub[which(d.sub$COND==3),3:25])

corr <- cor(subset,use="pairwise.complete.obs")

eigen(corr)$values

plot(eigen(corr)$values)

fa(corr,3,rotate="promax")


##############################################################################################################

d.long <- reshape(data        = d.sub,
                  idvar       = "ID",
                  varying     = list(colnames(d.sub)[3:27]),
                  timevar     = "Item",
                  times       = 1:25,
                  v.names      = "RT",
                  direction   = "long")

d.long1 <- d.long[which(d.long$COND==1),]
d.long2 <- d.long[which(d.long$COND==2),]
d.long3 <- d.long[which(d.long$COND==3),]

d.long1$id <- 1:33
d.long2$id <- 1:30
d.long3$id <- 1:30


d.long1$i.status  <- ifelse((d.long1$Item)%%2 == 1,0,1)
d.long2$i.status  <- ifelse((d.long2$Item)%%2 == 1,0,1)
d.long3$i.status  <- ifelse((d.long3$Item)%%2 == 1,0,1)

d.long1 <- na.omit(d.long1)
d.long2 <- na.omit(d.long2)
d.long3 <- na.omit(d.long3)


id1            <- d.long1$id
item1          <- d.long1$Item
istatus1       <- d.long1$i.status
rt1            <- log(d.long1$RT)
n_obs1         <- length(rt1)

id2            <- d.long2$id
item2          <- d.long2$Item
istatus2       <- d.long2$i.status
rt2            <- log(d.long2$RT)
n_obs2         <- length(rt2)

id3            <- d.long3$id
item3          <- d.long3$Item
istatus3       <- d.long3$i.status
rt3            <- log(d.long3$RT)
n_obs3         <- length(rt3)

###########################################################################################################################

code_rt<-'

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

	int  i3;                      //  number of individuals in condition 3
	int  n_obs3;                  //  number of observed responses in condition 3
	int  istatus3[n_obs3];        //  item indicator of a response in condition 3
      int  item3[n_obs3];           //  item id in condition 3
      int  id3[n_obs3];             //  person id in condition 3
      real Y3[n_obs3];              //  vector of log of response times in condition 3

      vector[2] mu_tau1;  // item specific residual standard deviation,
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

	vector[2] tau3[i3];         // person specific latent working speed for uncompromised items Condition 1

    	vector[2] mu_tau2;
    	vector[2] mu_tau3;

	corr_matrix[2] omega1;
	corr_matrix[2] omega2;
	corr_matrix[2] omega3;

	vector<lower=0>[2] sigma1;
	vector<lower=0>[2] sigma2;
	vector<lower=0>[2] sigma3;
}


model{

	sigma1   ~ cauchy(0,2.5);
	omega1   ~ lkj_corr(1);

	tau1    ~ multi_normal(mu_tau1,quad_form_diag(omega1, sigma1)); // http://stla.github.io/stlapblog/posts/StanLKJprior.html
                                                                            // http://modernstatisticalworkflow.blogspot.com/2016/09/what-are-covariance-matrices-anyway.html

	mu_tau2  ~ normal(0,5);
	sigma2   ~ cauchy(0,2.5);
	omega2   ~ lkj_corr(1);
	tau2    ~ multi_normal(mu_tau2,quad_form_diag(omega2, sigma2));

	mu_tau3  ~ normal(0,5);
	sigma3   ~ cauchy(0,2.5);
	omega3   ~ lkj_corr(1);
	tau3    ~ multi_normal(mu_tau3,quad_form_diag(omega3, sigma3));

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

		Y1[i] ~ normal(p2,1/(alpha[item1[i]]^2));
      }



	for (i in 1:n_obs2) {        

		real p_t = beta[item2[i]]-tau2[id2[i],1];
		real p_c = beta[item2[i]]-tau2[id2[i],2];

		real p2 = (1-istatus2[i])*p_t +(istatus2[i])*p_c;

		Y2[i] ~ normal(p2,1/(alpha[item2[i]]^2));
      }

	for (i in 1:n_obs3) {        

		real p_t = beta[item3[i]]-tau3[id3[i],1];
		real p_c = beta[item3[i]]-tau3[id3[i],2];

		real p2 = (1-istatus3[i])*p_t +(istatus3[i])*p_c;

		Y3[i] ~ normal(p2,1/(alpha[item3[i]]^2));
      }
}

generated quantities{

  real rt1_rep[n_obs1];
  real rt2_rep[n_obs2];
  real rt3_rep[n_obs3];


	for (i in 1:n_obs1) {    

		real p_t = beta[item1[i]]-tau1[id1[i],1];
		real p_c = beta[item1[i]]-tau1[id1[i],1];

		real p2 = (1-istatus1[i])*p_t +(istatus1[i])*p_c;

		rt1_rep[i] = normal_rng(p2,1/(alpha[item1[i]]^2));
      }

  for (i in 1:n_obs2) {        

		real p_t = beta[item2[i]]-tau2[id2[i],1];
		real p_c = beta[item2[i]]-tau2[id2[i],2];

		real p2 = (1-istatus2[i])*p_t +(istatus2[i])*p_c;

		rt2_rep[i] = normal_rng(p2,1/(alpha[item2[i]]^2));
      }

  for (i in 1:n_obs3) {        

		real p_t = beta[item3[i]]-tau3[id3[i],1];
		real p_c = beta[item3[i]]-tau3[id3[i],2];

		real p2 = (1-istatus3[i])*p_t +(istatus3[i])*p_c;

		rt3_rep[i] = normal_rng(p2,1/(alpha[item3[i]]^2));
  }
}
'


data_rt<-list(J        = 25,
              i1       = 33,
              n_obs1   = length(rt1),
              item1    = item1,
              id1      = id1,
              istatus1 = istatus1,
              Y1       = rt1,
              i2       = 30,
              n_obs2   = length(rt2),
              item2    = item2,
              id2      = id2,
              istatus2 = istatus2,
              Y2       = rt2,
              i3       = 30,
              n_obs3   = length(rt3),
              item3    = item3,
              id3      = id3,
              istatus3 = istatus3,
              Y3       = rt3,
              mu_tau1  = c(0,0))


rt_stan <- stan(model_code = code_rt, 
                data       = data_rt, 
                iter       = 10000,
                chains     = 4,
                control    = list(max_treedepth = 15,adapt_delta=0.99))

################################################################################
# Load the workspace with fitted model

load(here("data/Real Data/run1_10000_4chain_withcorrelation"))

################################################################################

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


plot(rt_stan, 
     plotfun = "hist", 
     pars = c("beta[1]","beta[2]","beta[3]","beta[4]","beta[5]","beta[6]","beta[7]","beta[8]","beta[9]","beta[10]",
              "beta[11]","beta[12]","beta[13]","beta[14]","beta[15]","beta[16]","beta[17]","beta[18]","beta[19]","beta[20]",
              "beta[21]","beta[22]","beta[23]","beta[24]","beta[25]"), 
     inc_warmup = FALSE,
     col="black")


plot(rt_stan, 
     plotfun = "hist", 
     pars = c("alpha[1]","alpha[2]","alpha[3]","alpha[4]","alpha[5]","alpha[6]","alpha[7]","alpha[8]","alpha[9]","alpha[10]",
              "alpha[11]","alpha[12]","alpha[13]","alpha[14]","alpha[15]","alpha[16]","alpha[17]","alpha[18]","alpha[19]","alpha[20]",
              "alpha[21]","alpha[22]","alpha[23]","alpha[24]","alpha[25]"), 
     inc_warmup = FALSE)



alphas <- extract(rt_stan,pars = c("alpha[1]","alpha[2]","alpha[3]","alpha[4]","alpha[5]","alpha[6]","alpha[7]","alpha[8]","alpha[9]","alpha[10]",
                                    "alpha[11]","alpha[12]","alpha[13]","alpha[14]","alpha[15]","alpha[16]","alpha[17]","alpha[18]","alpha[19]","alpha[20]",
                                    "alpha[21]","alpha[22]","alpha[23]","alpha[24]","alpha[25]"),permuted=TRUE)


i= 25

round(mean(alphas[[i]]^2),2)
round(sd(alphas[[i]]^2),2)
round(quantile(alphas[[i]]^2,c(.025,.975)),2)

a <- c()
for(i in 1:25){
  a[i] <- mean(alphas[[i]]^2)
}



posterior <- as.array(rt_stan)

dimnames(posterior)[[3]][1:25] <- c("Item 1","Item 2","Item 3","Item 4","Item 5","Item 6","Item 7","Item 8","Item 9","Item 10",
                                    "Item 11","Item 12","Item 13","Item 14","Item 15","Item 16","Item 17","Item 18","Item 19",
                                    "Item 20","Item 21","Item 22","Item 23","Item 24","Item 25")


dimnames(posterior)[[3]][1:25] <- c("beta[1]","beta[2]","beta[3]","beta[4]","beta[5]","beta[6]","beta[7]","beta[8]","beta[9]","beta[10]",
              "beta[11]","beta[12]","beta[13]","beta[14]","beta[15]","beta[16]","beta[17]","beta[18]","beta[19]","beta[20]",
              "beta[21]","beta[22]","beta[23]","beta[24]","beta[25]")


dimnames(posterior)[[3]][28:52] <- c("Item 1","Item 2","Item 3","Item 4","Item 5","Item 6","Item 7","Item 8","Item 9","Item 10",
                                    "Item 11","Item 12","Item 13","Item 14","Item 15","Item 16","Item 17","Item 18","Item 19",
                                    "Item 20","Item 21","Item 22","Item 23","Item 24","Item 25")



color_scheme_set("darkgray")

p <- mcmc_hist(posterior, 
     pars = c("Item 1","Item 2","Item 3","Item 4","Item 5","Item 6","Item 7","Item 8","Item 9","Item 10",
                                    "Item 11","Item 12","Item 13","Item 14","Item 15","Item 16","Item 17","Item 18","Item 19",
                                    "Item 20","Item 21","Item 22","Item 23","Item 24","Item 25"))

p + theme_bw() +
  theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()) 


tau1  <- summary(rt_stan)$summary[55:120,]
tau2  <- summary(rt_stan)$summary[121:180,]
tau3  <- summary(rt_stan)$summary[181:240,]

par(mfrow=c(3,2))

hist(tau1[,10])
hist(tau2[,10])
hist(tau3[,10])

tau1_ <- cbind(tau1[seq(1,66,2),1],tau1[seq(2,66,2),1])
tau2_ <- cbind(tau2[seq(1,60,2),1],tau2[seq(2,60,2),1])
tau3_ <- cbind(tau3[seq(1,60,2),1],tau3[seq(2,60,2),1])


describe(tau1_)

tau1_[,1] = (tau1_[,1]-mean(tau1_[,1]))*(0.176/sd(tau1_[,1]))
tau1_[,2] = (tau1_[,2]-mean(tau1_[,2]))*(0.174/sd(tau1_[,2]))
describe(tau1_)

cor(tau1_)

describe(tau2_)
tau2_[,1] = ((tau2_[,1]-mean(tau2_[,1]))*(0.549/sd(tau2_[,1])))+0.05
tau2_[,2] = ((tau2_[,2]-mean(tau2_[,2]))*(0.474/sd(tau2_[,2])))+1.34
describe(tau2_)
cor(tau2_)


describe(tau3_)
tau3_[,1] = ((tau3_[,1]-mean(tau3_[,1]))*(0.430/sd(tau3_[,1])))-0.18
tau3_[,2] = ((tau3_[,2]-mean(tau3_[,2]))*(0.457/sd(tau3_[,2])))+1.47
describe(tau3_)
cor(tau3_)

tau <- data.frame(rbind(tau1_,tau2_,tau3_))

plot(tau[1:33,1],tau[1:33,2],xlim=c(-.75,2.75),ylim=c(-.75,2.75),cex=.8,cex.lab=.8,
     xlab="Latent Speed Parameter Estimate - Undisclosed Items",
     ylab="Latent Speed Parameter Estimate - Disclosed Items")

points(tau[34:63,1],tau[34:63,2],pch=2,cex=.8)
points(tau[64:93,1],tau[64:93,2],pch=3,cex=.8)
abline(coef=c(0,1),lty=2,col="gray")

legend("bottomright",c("Control Group","Items Disclosed","Items and Answers Disclosed"),
       pch=c(1,2,3),cex=.7)


########################################################################################################################

# Model Fit Assessment

# Model generated posterior predictions

  # Each has 100,000 rows (25,000 iterations x 4 chains)
  # Each row indicates a separate generated dataset
  
	rt1_rep <- extract(rt_stan,pars="rt1_rep",permuted=TRUE)$rt1_rep
	dim(rt1_rep)

	length(rt1)  # observed data, rt1 vs rt1_rep to be compared


	rt2_rep <- extract(rt_stan,pars="rt2_rep",permuted=TRUE)$rt2_rep
	dim(rt2_rep)

	length(rt2)  # observed data, rt2 vs rt2_rep to be compared

	rt3_rep <- extract(rt_stan,pars="rt3_rep",permuted=TRUE)$rt3_rep
	dim(rt3_rep)

	length(rt3)  # observed data, rt3 vs rt3_rep to be compared


layout(matrix(c(1,1,1,1,1,0,2,2,2,2,2,0,0,0,3,3,3,3,3,0,0,0), 
                nrow=2, 
                ncol=11, 
                byrow = TRUE))


# The lower-tail posterior predictive p-values of the observed logtimes

  # for all items and test takers in Condition 1

		pis1 <- c()

		for(i in 1:length(rt1)){
		   pis1[i]=sum(rt1_rep[,i] < rt1[i])/nrow(rt1_rep)
		}

		pis1 <- round(pis1,3)

		myTable <- data.frame( table(pis1))
		myTable$Prop <- prop.table( myTable$Freq )
		myTable$CumProp <-  cumsum( myTable$Prop )
	 
            myTable
	
		hist(pis1)

		plot(as.numeric(as.character(myTable$pis1)),myTable$CumProp,type="l")
		abline(coef=c(0,1))


  # for all items and test takers in Condition 2

		pis2 <- c()

		for(i in 1:length(rt2)){
		   pis2[i]=sum(rt2_rep[,i] < rt2[i])/nrow(rt2_rep)
		}

		pis2 <- round(pis2,3)

		myTable <- data.frame( table(pis2))
		myTable$Prop <- prop.table( myTable$Freq )
		myTable$CumProp <-  cumsum( myTable$Prop )
	 
            myTable
	
		hist(pis2)

		plot(as.numeric(as.character(myTable$pis2)),myTable$CumProp,type="l")
		abline(coef=c(0,1))

  # for all items and test takers in Condition 3

		pis3 <- c()

		for(i in 1:length(rt3)){
		   pis3[i]=sum(rt3_rep[,i] < rt3[i])/nrow(rt3_rep)
		}

		pis3 <- round(pis3,3)

		myTable <- data.frame( table(pis3))
		myTable$Prop <- prop.table( myTable$Freq )
		myTable$CumProp <-  cumsum( myTable$Prop )
	 
            myTable
	
		hist(pis3)

		plot(as.numeric(as.character(myTable$pis3)),myTable$CumProp,type="l")
		abline(coef=c(0,1))


pis <- c(pis1,pis2,pis3)

myTable <- data.frame( table(pis))
		myTable$Prop <- prop.table( myTable$Freq )
		myTable$CumProp <-  cumsum( myTable$Prop )
	 
plot(as.numeric(as.character(myTable$pis)),myTable$CumProp,type="l",
     main="All 25 items and 93 test takers",xlab=expression(paste(pi,"-value")),
     ylab="Cumulative proportion",cex.main=0.8)
abline(coef=c(0,1),lty=2,col="gray")

	      temp1    = data.frame(cbind(id1,item1,pis1))
             colnames(temp1) <- c("id","item","pis")

	      temp2    = data.frame(cbind(id2+33,item2,pis2))
             colnames(temp2) <- c("id","item","pis")

	      temp3    = data.frame(cbind(id3+63,item3,pis3))
             colnames(temp3) <- c("id","item","pis")


temp <- rbind(temp1,temp2,temp3)
temp$pis <- round(temp$pis,2)

# Each item across test takers


	i = 1 
		temp.tab         = data.frame(table(temp[which(temp$item==i),]$pis))
		temp.tab$prop    = prop.table(temp.tab$Freq)
		temp.tab$Cumprop = cumsum(temp.tab$prop)
		plot(as.numeric(as.character(temp.tab[,1])),temp.tab[,4],
                  type="l",xlab=expression(paste(pi,"-value")),ylab="Cumulative proportion",lty=2,col="gray",
			main="For each of the 25 items across all 93 test takers",cex.main=0.75,,xlim=c(0,1),ylim=c(0,1))
		
	for(i in 2:25){
		temp.tab         = data.frame(table(temp[which(temp$item==i),]$pis))
		temp.tab$prop    = prop.table(temp.tab$Freq)
		temp.tab$Cumprop = cumsum(temp.tab$prop)
		points(as.numeric(as.character(temp.tab[,1])),temp.tab[,4],type="l",lty=2,col="gray")
	} 

	abline(coef=c(0,1))

# Each test takers across items 

	i = 1 
		temp.tab         = data.frame(table(temp[which(temp$id==i),]$pis))
		temp.tab$prop    = prop.table(temp.tab$Freq)
		temp.tab$Cumprop = cumsum(temp.tab$prop)
		plot(as.numeric(as.character(temp.tab[,1])),temp.tab[,4],
                  type="l",xlab=expression(paste(pi,"-value")),ylab="Cumulative proportion",lty=2,col="gray",
			main="For each of the 93 test taker across all 25 items",cex.main=0.75,xlim=c(0,1),ylim=c(0,1))
		
	for(i in 2:93){
		temp.tab         = data.frame(table(temp[which(temp$id==i),]$pis))
		temp.tab$prop    = prop.table(temp.tab$Freq)
		temp.tab$Cumprop = cumsum(temp.tab$prop)
		points(as.numeric(as.character(temp.tab[,1])),temp.tab[,4],type="l",lty=2,col="gray")
	} 

	abline(coef=c(0,1))



	i = 10
		temp.tab         = data.frame(table(temp[which(temp$id==i),]$pis))
		temp.tab$prop    = prop.table(temp.tab$Freq)
		temp.tab$Cumprop = cumsum(temp.tab$prop)
		plot(as.numeric(as.character(temp.tab[,1])),temp.tab[,4],
                  type="l",xlab="",ylab="",xlim=c(0,1),ylim=c(0,1))

		abline(coef=c(0,1))
		
	


################################################################################################

source("https://raw.githubusercontent.com/betanalpha/knitr_case_studies/master/qr_regression/stan_utility.R")


hist(summary(rt_stan)$summary[1:262,10],main="",xlab="Rhat")
describe(summary(rt_stan)$summary[1:262,10])


check_treedepth(rt_stan)

check_energy(rt_stan)

check_div(rt_stan)


################################################################################################

beta  <- summary(rt_stan)$summary[1:25,1]
alpha <- summary(rt_stan)$summary[28:52,1]

tau10  <- c(tau1_[,1])
tau11  <- c(tau1_[,2])

tau20  <- tau2_[,1]
tau21  <- tau2_[,2]

tau30  <- tau3_[,1]
tau31  <- tau3_[,2]


pred      <- data.frame(cbind(id1,item1,rt1))
pred$pred <- NA
pred$sres <- NA

for(i in 1:nrow(pred)){

	if(pred[i,2]%%2==1){
		pred[i,]$pred = beta[pred[i,2]] - tau10[pred[i,1]]
		pred[i,]$sres = (pred[i,3] - (beta[pred[i,2]] - tau10[pred[i,1]]))/(1/(alpha[pred[i,2]]))
	} else {
		pred[i,]$pred = beta[pred[i,2]] - tau11[pred[i,1]]
		pred[i,]$sres = (pred[i,3] - (beta[pred[i,2]] - tau11[pred[i,1]]))/(1/(alpha[pred[i,2]]))

	}
} 

pred1 = pred
colnames(pred1) <- c("id","item","rt","pred","sres")
cor(pred1[,3],pred1[,4])


pred      <- data.frame(cbind(id2,item2,rt2))
pred$pred <- NA
pred$sres <- NA

for(i in 1:nrow(pred)){

	if(pred[i,2]%%2==1){
		pred[i,]$pred = beta[pred[i,2]] - tau20[pred[i,1]]
		pred[i,]$sres = (pred[i,3] - (beta[pred[i,2]] - tau20[pred[i,1]]))/(1/(alpha[pred[i,2]]))
	} else {
		pred[i,]$pred = beta[pred[i,2]] - tau21[pred[i,1]]
		pred[i,]$sres = (pred[i,3] - (beta[pred[i,2]] - tau21[pred[i,1]]))/(1/(alpha[pred[i,2]]))

	}
} 

pred2 = pred
colnames(pred2) <- c("id","item","rt","pred","sres")
pred2$id = pred2$id + 33

cor(pred2[,3],pred2[,4])



pred      <- data.frame(cbind(id3,item3,rt3))
pred$pred <- NA
pred$sres <- NA

for(i in 1:nrow(pred)){

	if(pred[i,2]%%2==1){
		pred[i,]$pred = beta[pred[i,2]] - tau30[pred[i,1]]
		pred[i,]$sres = (pred[i,3] - (beta[pred[i,2]] - tau30[pred[i,1]]))/(1/(alpha[pred[i,2]]))
	} else {
		pred[i,]$pred = beta[pred[i,2]] - tau31[pred[i,1]]
		pred[i,]$sres = (pred[i,3] - (beta[pred[i,2]] - tau31[pred[i,1]]))/(1/(alpha[pred[i,2]]))

	}
} 

pred3 = pred
colnames(pred3) <- c("id","item","rt","pred","sres")
pred3$id = pred3$id + 63

cor(pred3[,3],pred3[,4])

predict <- rbind(pred1,pred2,pred3)

cor(predict$rt,predict$pred)

plot(predict$rt,predict$pred)

plot(predict$rt,predict$sred)

describe(predict$sres)


predict2 <- predict[,c(-3,-4)]

	res.wide <- reshape(data=predict2,
                          idvar="id",
                          v.names="sres",
                          timevar="item",
                          direction="wide")
	
      res.wide <- res.wide[,-1]

	boxplot(res.wide,col="gray",xaxt="n",ylab="Standardized Residual",xlab="Item Location")
	axis(side=1,at=1:25)

	plot(1:25,colMeans(res.wide,na.rm=TRUE),ylim=c(-.1,.1),type="l",
          ylab="Mean Standardized Residual",xlab="Item Location",xaxt="n")

	axis(side=1,at=1:25)

      abline(h=0,lty=2,col="gray")
      abline(h=-0.5,lty=2,col="gray")

	plot(1:25,colMeans(res.wide,na.rm=TRUE),ylim=c(-.1,.1),type="l")

      abline(h=0,lty=2)


	res.cor <- cor(res.wide,use="pairwise.complete.obs")

    res.cor[1,2]
    res.cor[2,3]
    res.cor[3,4]
    res.cor[4,5]
    res.cor[5,6]
    res.cor[6,7]
    res.cor[7,8]
    res.cor[8,9]
    res.cor[9,10]
    res.cor[10,11]
    res.cor[11,12]
    res.cor[12,13]
    res.cor[13,14]
    res.cor[14,15]
    res.cor[15,16]
    res.cor[16,17]
    res.cor[17,18]
    res.cor[18,19]
    res.cor[19,20]
    res.cor[20,21]
    res.cor[21,22]
    res.cor[22,23]
    res.cor[23,24]
    res.cor[24,25]

	diag(res.cor) <- NA
      describe(res.cor)

	# Below -0.203 or above .203 is the significant correlation

      length(which(abs(res.cor)>.203))/600


require(MASS)

corr = c()

for(i in 1:100000){

     corr[i]=cor(mvrnorm(93,rep(0,2),diag(2)))[1,2]
}

describe(corr)

hist(corr)


qt(0.025/300,91)

r=0.38
(r*sqrt(91))/sqrt(1-r^2)

# .38 critical value for significant after controlling for
# multiple testing


res.cor2 = (res.cor>.38)*1
library(reshape2)
melted_cormat <- melt(res.cor)

ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()



# Get lower triangle of the correlation matrix
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }

upper_tri <- get_upper_tri(res.cor)
upper_tri
melted_cormat <- melt(upper_tri, na.rm=TRUE)

melted_cormat[,1] = as.numeric(substring(melted_cormat[,1],first=6))
melted_cormat[,2] = as.numeric(substring(melted_cormat[,2],first=6))

melted_cormat = melted_cormat[which(abs(melted_cormat[,3])>.38),]

ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile(color="black")+theme_bw()+xlab("Item Location")+ylab("Item Location")+
  scale_x_continuous(breaks=1:25,limits=c(.5,25.5))+
  scale_y_continuous(breaks=1:25,limits=c(.5,25.5))+
  scale_fill_continuous(name = "")+
  scale_fill_gradientn(colours = c("white", "black"), values = c(0,1))





hist(upper_tri,main="",lty=2,xlab="Residual Correlations")
abline(v=.38)
abline(v=-.38)

describe(upper_tri)


length(which(melted_cormat[,1]%%2==0 & melted_cormat[,2]%%2==0))


###############################################################################################

require(MASS)

beta  <- rnorm(25,4.19,0.38)
alpha <- rnorm(25,1.18,0.16)

cor1 <- matrix(c(1,.84,.84,1),2,2)
tau1 <- mvrnorm(33,c(0,0),cor2cov(cor1,c(0.18,0.17)))

cor2 <- matrix(c(1,.56,.56,1),2,2)
tau2 <- mvrnorm(33,c(0.05,1.34),cor2cov(cor2,c(0.55,0.47)))

C    <- c(0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0)

rt1 <- matrix(nrow=33,ncol=25)

for(i in 1:33){
 for(j in 1:25){
  
   p_t = beta[j]-tau1[i,1]
   p_c = beta[j]-tau1[i,2]
 
   p   = p_t*(1-C[j]) + p_c*C[j]

   rt1[i,j] = rnorm(1,p,1/alpha[j]^2)   
 }
}

rt2 <- matrix(nrow=30,ncol=25)

for(i in 1:30){
 for(j in 1:25){
  
   p_t = beta[j]-tau2[i,1]
   p_c = beta[j]-tau2[i,2]
 
   p   = p_t*(1-C[j]) + p_c*C[j]

   rt2[i,j] = rnorm(1,p,1/alpha[j]^2)   
 }
}

rt <- rbind(cbind(data.frame(exp(rt1)),gr=1),
            cbind(data.frame(exp(rt2)),gr=2))

require(psych)

describeBy(rt[,1:25],rt$gr)

























