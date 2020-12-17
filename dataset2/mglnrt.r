
d.long <- read.csv('/gpfs/projects/edquant/cengiz/mglnrt3/rt_long.csv')

d.long <- na.omit(d.long)

d.long$item <- as.numeric(substring(d.long$Item,6))



describeBy(exp(d.long[which(d.long$i_flag=='Flagged'),]$RT),
           list(d.long[which(d.long$i_flag=='Flagged'),]$item),mat=TRUE)[,c(2,4,5,6)]

describeBy(exp(d.long[which(d.long$i_flag=='Unflagged'),]$RT),
           list(d.long[which(d.long$i_flag=='Unflagged'),]$item),mat=TRUE)[,c(2,4,5,6)]


icode = c(201,204,205,206,208,210,211,214,216,218,220,222,223,224,226,227,228,
          229,230,231,232,234,236,237,239,241,242,243,244,247,248,249,253,255,
          257,258,261,263,264,265,267,268,269,272,279,280,281,282,283,286,287,
          288,292,296,297,298,302,303,305,307,308,314,316,318,322,323,327,329,
          331,338,340,341,343,345,346,347,348,350,353,359,366,367,368)

icode2 = 171:253


for(i in 1:83){
  d.long[which(d.long$item == icode[i]),]$item = icode2[i]
}

unique(d.long$item)

tab <- as.data.frame(table(d.long$Item,d.long$item))
tab = tab[which(tab[,3]!=0),]



d.long$id <- NA

ind = which(substring(d.long$EID,1,3)=='e10')
d.long[ind,]$id = as.numeric(substring(d.long[ind,]$EID,4,7))

ind = which(substring(d.long$EID,1,3)=='e20')
d.long[ind,]$id = as.numeric(substring(d.long[ind,]$EID,4,7))+1636

tab <- as.data.frame(table(d.long$EID,d.long$id))
tab = tab[which(tab[,3]!=0),]

########################################################################

d.long1 <- d.long[which(d.long$p_flag==0),]
d.long2 <- d.long[which(d.long$p_flag==1),]


d.long1$id2 <- NA

icode = unique(d.long1$id)
for(i in 1:3186){
  d.long1[which(d.long1$id == icode[i]),]$id2 = i
}

tab <- as.data.frame(table(d.long1$id,d.long1$id2))
tab = tab[which(tab[,3]!=0),]


d.long2$id2 <- NA

icode = unique(d.long2$id)
for(i in 1:94){
  d.long2[which(d.long2$id == icode[i]),]$id2 = i
}

tab <- as.data.frame(table(d.long2$id,d.long2$id2))
tab = tab[which(tab[,3]!=0),]


id1            <- d.long1$id2
item1          <- d.long1$item
istatus1       <- ifelse(d.long1$i_flag=='Flagged',1,0)
rt1            <- d.long1$RT
n_obs1         <- length(rt1)

id2            <- d.long2$id2
item2          <- d.long2$item
istatus2       <- ifelse(d.long2$i_flag=='Flagged',1,0)
rt2            <- d.long2$RT
n_obs2         <- length(rt2)

###########################################################################################################################

require(cmdstanr)

set_cmdstan_path('/gpfs/home/cengiz/.cmdstanr/cmdstan-2.24.1')

mod <- cmdstan_model('/gpfs/projects/edquant/cengiz/mglnrt3/mglnrt.stan')


mod <- cmdstan_model('B:/Ongoing_Research/Murat/ResponseTime/EMIP/Real Data 2/mglnrt.stan')


data_rt<-list(J        = 253,
              i1       = 3186,
              n_obs1   = length(rt1),
              item1    = item1,
              id1      = id1,
              istatus1 = istatus1,
              Y1       = rt1,
              i2       = 94,
              n_obs2   = length(rt2),
              item2    = item2,
              id2      = id2,
              istatus2 = istatus2,
              Y2       = rt2,
              mu_tau1  = c(0,0))

fit <- mod$sample(
  data = data_rt,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  iter_warmup   = 500,
  iter_sampling = 1500,
  refresh = 100
)


fit$save_output_files('/gpfs/projects/edquant/cengiz/mglnrt3')

stanfit <- rstan::read_stan_csv(fit$output_files())

save.image("/gpfs/projects/edquant/cengiz/mglnrt3/mglnrt.RData")

######################################################################


print(stanfit, 
      pars = c("beta","mu_beta","sigma_beta",
               "alpha","mu_alpha","sigma_alpha"),
      probs = c(0.025,0.975), 
      digits = 3)


print(stanfit, 
      pars = c("mu_tau2",
               "sigma1","sigma2",
               "omega1","omega2"),
      probs = c(0.025,0.975), 
      digits = 3)



alphas <- extract(stanfit,pars = 'alpha',permuted=TRUE)


a <- matrix(nrow=253,ncol=4)

for(i in 1:253){
  a[i,1]   <- mean(alphas[[1]][,i]^2)
  a[i,2]   <- sd(alphas[[1]][,i]^2)
  a[i,3:4] <- quantile(alphas[[1]][,i]^2,c(.025,.975))
}


#######################################################################


tau1  <- summary(stanfit)$summary[511:6882,]
tau2  <- summary(stanfit)$summary[6883:7070,]


tau1_ <- cbind(tau1[seq(1,6372,2),1],tau1[seq(2,6372,2),1])
tau2_ <- cbind(tau2[seq(1,188,2),1],tau2[seq(2,188,2),1])

plot(tau1_[,1],tau1_[,2],xlim=c(-1,1),ylim=c(-1,1),cex=.8,cex.lab=.8,
     xlab="Latent Speed Parameter Estimate - Unflagged Items",
     ylab="Latent Speed Parameter Estimate - Flagged Items",
     main='Unflagged Test Takers')

plot(tau2_[,1],tau2_[,2],xlim=c(-1,1),ylim=c(-1,1),cex=.8,cex.lab=.8,
     xlab="Latent Speed Parameter Estimate - Unflagged Items",
     ylab="Latent Speed Parameter Estimate - Flagged Items",
     main='Flagged Test Takers')

abline(coef=c(0,1),lty=2,col="gray")

legend("bottomright",c("Unflagged Test Takers","Flagged Test Takers"),
       pch=c(1,2),cex=.7)


t.test(tau2_[,1],tau2_[,2],paired=TRUE)


effsize::cohen.d(tau2_[,1],tau2_[,2],paired=T,within=F)


describe(summary(stanfit)$summary[,10])

hist(summary(stanfit)$summary[,10],col='white',xlab='Rhat',main='')

########################################################################

r = 1

a = extract(stanfit,pars=paste0("tau2[",r,",1]"),permuted=TRUE)[[1]]
b = extract(stanfit,pars=paste0("tau2[",r,",2]"),permuted=TRUE)[[1]]


tau1  <- summary(stanfit)$summary[511:6882,]
tau2  <- summary(stanfit)$summary[6883:7070,]


a = extract(stanfit,pars='tau2',permuted=TRUE)[[1]]

a1 = a[,,1]
a2 = a[,,2]

aa = rowMeans(a1) - rowMeans(a2)
hist(aa)

round(quantile(aa,c(.025,.975)),3)
mean(aa)

########################################################################

# Model Fit Assessment

# Model generated posterior predictions

  # Each has 8000 rows (2000 iterations x 4 chains)
  # Each row indicates a separate generated dataset

  # summary(stanfit)$summary[7085:548614,] #rt1_rep
  # summary(stanfit)$summary[548615:564446,] #rt1_rep

  # Group 1
	# observed data, rt1 vs rt1_rep to be compared

	# rt1[1]  -- the first observation
	# rt1_rep[1] -- 6000 posterior draws for this observation
                      # 1500 iter per chain, 2000 total iter - 500 warmup iter
                      # 1500 iter, 4 chains

	pis <- c()

	for(r in 1:541530){
		kkk = extract(stanfit,pars=paste0("rt1_rep[",r,"]"),permuted=TRUE)
		pis[r] = sum(kkk[[1]] < rt1[r])/length(kkk[[1]])
		print(r)
	}

  # Group 2

	# observed data, rt2 vs rt2_rep to be compared

	# rt2[1]  -- the first observation for the second group
	# rt2_rep[1] -- 6000 posterior draws for this observation
                      # 1500 iter per chain, 2000 total iter - 500 warmup iter
                      # 1500 iter, 4 chains

	pis2 <- c()

	for(r in 1:15832){
		kkk = extract(stanfit,pars=paste0("rt2_rep[",r,"]"),permuted=TRUE)
		pis2[r] = sum(kkk[[1]] < rt2[r])/length(kkk[[1]])
		print(r)
	}


layout(matrix(c(1,1,1,1,1,0,2,2,2,2,2,0,0,0,3,3,3,3,3,0,0,0), 
                nrow=2, 
                ncol=11, 
                byrow = TRUE))


# The lower-tail posterior predictive p-values of the observed logtimes

  # for all items and test takers in Condition 1

		pis1 <- round(pis1,3)

		myTable <- data.frame( table(pis1))
		myTable$Prop <- prop.table( myTable$Freq )
		myTable$CumProp <-  cumsum( myTable$Prop )
	 
            myTable
	
		hist(pis1)

		plot(as.numeric(as.character(myTable$pis1)),myTable$CumProp,type="l")
		abline(coef=c(0,1))


  # for all items and test takers in Condition 2

		pis2 <- round(pis2,3)

		myTable <- data.frame( table(pis2))
		myTable$Prop <- prop.table( myTable$Freq )
		myTable$CumProp <-  cumsum( myTable$Prop )
	 
            myTable
	
		hist(pis2)

		plot(as.numeric(as.character(myTable$pis2)),myTable$CumProp,type="l")
		abline(coef=c(0,1))



		pis <- c(pis1,pis2)
		myTable         <- data.frame( table(pis))
		myTable$Prop    <- prop.table( myTable$Freq )
		myTable$CumProp <-  cumsum( myTable$Prop )

		plot(as.numeric(as.character(myTable$pis)),myTable$CumProp,type="l",
                 xlab=expression(paste(pi,"-value")),
                 ylab="Cumulative proportion",
                 main="All 253 items and 3280 test takers",
                 cex.main=0.8)

		abline(coef=c(0,1),lty=2)






  # head(id1)
  # head(item1)	


  # Individuals

	ind.pis1 <- vector('list',length(unique(id1)))	
	
	for(i in 1:length(unique(id1))){
		ind.pis1[[i]] = pis1[which(id1==i)]
		print(i)
	}
	

	a <- round(ind.pis1[[1]],3)

	myTable         <- data.frame( table(a))
	myTable$Prop    <- prop.table( myTable$Freq )
	myTable$CumProp <-  cumsum( myTable$Prop )
	 
	plot(as.numeric(as.character(myTable$a)),myTable$CumProp,type="l",col='gray',
           main = 'For each of the 3280 test-takers across 253 items',
	     xlab=expression(paste(pi,"-value")),
	     ylab="Cumulative proportion",)
	
	for(i in 2:3186){
		a <- round(ind.pis1[[i]],3)
		myTable         <- data.frame( table(a))
		myTable$Prop    <- prop.table( myTable$Freq )
		myTable$CumProp <-  cumsum( myTable$Prop )
	 
		points(as.numeric(as.character(myTable$a)),myTable$CumProp,type="l",col='gray')
	}

	abline(coef=c(0,1))


  # Items

	ind.pis1 <- vector('list',length(unique(item1)))	
	
	for(i in 1:length(unique(item1))){
		ind.pis1[[i]] = pis1[which(item1==i)]
		print(i)
	}
	

	a <- round(ind.pis1[[1]],3)

	myTable         <- data.frame( table(a))
	myTable$Prop    <- prop.table( myTable$Freq )
	myTable$CumProp <-  cumsum( myTable$Prop )
	 
	plot(as.numeric(as.character(myTable$a)),myTable$CumProp,type="l",col='gray',
           main = 'For each of the 253 items across 3280 test-takers',
	     xlab=expression(paste(pi,"-value")),
	     ylab="Cumulative proportion",)
	
	for(i in 2:253){
		a <- round(ind.pis1[[i]],3)
		myTable         <- data.frame( table(a))
		myTable$Prop    <- prop.table( myTable$Freq )
		myTable$CumProp <-  cumsum( myTable$Prop )
	 
		points(as.numeric(as.character(myTable$a)),myTable$CumProp,type="l",col='gray')
	}

	abline(coef=c(0,1))

########################################################################

# Residuals


beta  <- summary(stanfit)$summary[1:253,1]
alpha <- summary(stanfit)$summary[256:508,1]

tau10  <- c(tau1_[,1])
tau11  <- c(tau1_[,2])

tau20  <- tau2_[,1]
tau21  <- tau2_[,2]


pred      <- data.frame(cbind(id1,item1,istatus1,rt1))
pred$pred <- NA
pred$sres <- NA

for(i in 1:nrow(pred)){

	if(pred[i,3]==0){
		pred[i,]$pred = beta[pred[i,2]] - tau10[pred[i,1]]
		pred[i,]$sres = (pred[i,4] - (beta[pred[i,2]] - tau10[pred[i,1]]))/(1/(alpha[pred[i,2]]))
	} else {
		pred[i,]$pred = beta[pred[i,2]] - tau11[pred[i,1]]
		pred[i,]$sres = (pred[i,4] - (beta[pred[i,2]] - tau11[pred[i,1]]))/(1/(alpha[pred[i,2]]))

	}

	print(i)
} 


pred1 = pred
colnames(pred1) <- c("id","item","istatus","rt","pred","sres")
cor(pred1[,4],pred1[,5])


pred      <- data.frame(cbind(id2,item2,istatus2,rt2))
pred$pred <- NA
pred$sres <- NA

for(i in 1:nrow(pred)){

	if(pred[i,3]==0){
		pred[i,]$pred = beta[pred[i,2]] - tau20[pred[i,1]]
		pred[i,]$sres = (pred[i,4] - (beta[pred[i,2]] - tau20[pred[i,1]]))/(1/(alpha[pred[i,2]]))
	} else {
		pred[i,]$pred = beta[pred[i,2]] - tau21[pred[i,1]]
		pred[i,]$sres = (pred[i,4] - (beta[pred[i,2]] - tau21[pred[i,1]]))/(1/(alpha[pred[i,2]]))

	}

	print(i)
} 



pred2 = pred
colnames(pred2) <- c("id","item","istatus","rt","pred","sres")
pred2$id = pred2$id + 3186

cor(pred2[,4],pred2[,5])



predict <- rbind(pred1,pred2)

cor(predict$rt,predict$pred)

plot(predict$rt,predict$pred)

plot(predict$rt,predict$sred)

describe(predict$sres)


predict2 <- predict[,c(-3,-4,-5)]

predict2 <- predict2[order(predict2$item),]

	res.wide <- reshape(data=predict2,
                          idvar="id",
                          v.names="sres",
                          timevar="item",
                          direction="wide")
	
      res.wide <- res.wide[,-1]

	plot(1:253,colMeans(res.wide,na.rm=TRUE),ylim=c(-.01,.01),type="l",
          ylab="Mean Standardized Residual",xlab="Item Location",xaxt="n",col='gray')

	axis(side=1,at=c(1,seq(20,250,20)),cex.axis=0.8)

      abline(h=0,lty=2)
  
	boxplot(res.wide,col="gray",xaxt="n",ylab="Standardized Residual",xlab="Item Location")
	axis(side=1,at=1:25)


max(colMeans(res.wide,na.rm=TRUE))
min(colMeans(res.wide,na.rm=TRUE))


#######################################################################

res.cor <- cor(res.wide,use="pairwise.complete.obs")

diag(res.cor) <- NA
describe(res.cor)

	# Below -0.048 or above 0.048 is the significant correlation

      length(which(abs(res.cor)>.048))/600



qt(0.025/24989,1650)

r=0.115
(r*sqrt(1650))/sqrt(1-r^2)

# .11 critical value for significant after controlling for
# multiple testing

res.cor2 = (res.cor>.11)*1
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

melted_cormat$ind = as.factor(ifelse(abs(melted_cormat[,3])>.115,'Significant','Not Significant'))

ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=ind)) + 
  geom_tile()+
  theme_bw()+
  scale_fill_manual(name="",values = c("white","black"))+
  xlab("Item Location")+
  ylab("Item Location")

dim(melted_cormat)

table(melted_cormat$ind)


hist(upper_tri,main="",lty=2,xlab="Residual Correlations",col='white')
abline(v=.115)
abline(v=-.115)

describe(upper_tri)


length(which(melted_cormat[,1]%%2==0 & melted_cormat[,2]%%2==0))

head(melted_cormat)

describe(melted_cormat$value)

################################################################################################

source("https://raw.githubusercontent.com/betanalpha/knitr_case_studies/master/qr_regression/stan_utility.R")


check_treedepth(stanfit)

check_energy(stanfit)

check_div(stanfit)


################################################################################################
















