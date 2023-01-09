seed_id = 1982
set.seed(seed_id)

install.packages("pacman")
pacman::p_load(R2jags, parallel, polspline, ggplot2, glue, hesim, VGAM)

# defining a function for calculating the maximum of the posterior density (not exactly the same as the mode)
MPD <- function(x) {
  density(x)$x[which(density(x)$y==max(density(x)$y))]
}


##### Preprocessing #####

groupSize <- 4
ntrials <- 10
pi <- 1.6 # used to be 1.4, but the original paper (and Josh' preprint both say 1.6)
ntokens <- 20
vals <- seq(0,ntokens,1)
#vals <- seq(1,21,1) #possible values to contribute - from 0 to 20 tokens

rawDat <- read.csv("data/HerrmannThoeniGaechterDATA.csv", skip = 3) # Public goods game

# extract every third line - data file has lines representing others responses and we don't need that
redDat <- rawDat[seq(1,length(rawDat$sessionid),3),]

group_names <- unique(redDat$groupid)
ngroups <- length(group_names)

# THIS WILL REMOVE SUBJECTS WITH MISSING DATA IN NO PUNISHMENT CONDITION
ngroups <- 269

subject_names <- unique(redDat$subjectid)
nsubjects <- length(subject_names)

# data for no punishment condition #
membership_no_punish <- array(0,c(groupSize,ngroups))
c_no_punish <- array(0,c(groupSize,ntrials,ngroups)) # choices
Gga_no_punish <- array(0,c(ntrials,ngroups)) # group-averaged contribution (only 1 entry per group - hence the G + ga)
Ggas_no_punish <- array(0,c(groupSize,ntrials,ngroups)) # same as Gga, but specified for each subject (cuz that's how the JAGS-code wants it) - hence the s
Gc_no_punish <- array(0,c(groupSize,ntrials,ngroups)) # summed ("cumulated" hence the c) group contribution not including oneself (therefore specified for each subject) - we don't use this - this refers to the sum-command in the loop-structure
Ga_no_punish <- array(0,c(groupSize,ntrials,ngroups)) # group-averaged contribution without oneself - we don't use this cuz the participants don't see this (hence, they have to do quite a bit of mental arithmetics to represent this), and thus we're modeling their conditional preference relative to the averaged group contribution (including their own)
missing <- array(0,ngroups)
member_NA_no_punish <- array(FALSE,ngroups)

for (g in 1:ngroups) {
  c_no_punish[,,g] <- rbind(redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="N-experiment"][1:10],
                            redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="N-experiment"][11:20],
                            redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="N-experiment"][21:30],
                            redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="N-experiment"][31:40])
  
  Gga_no_punish[,g] <- colMeans(c_no_punish[,,g])
  
  missing[g] <- is.na(c_no_punish[1,1,g])
  
  for (s in 1:groupSize) {
    membership_no_punish[s,g] <- redDat$membership[redDat$groupid==group_names[g]&redDat$p=="N-experiment"][s+(s-1)*10]
    Gc_no_punish[s,,g] <- colSums(c_no_punish[-s,,g])
    Ga_no_punish[s,,g] <- colMeans(c_no_punish[-s,,g])
    Ggas_no_punish[s,,g] <- colMeans(c_no_punish[,,g])
  }
  member_NA_no_punish[g] <- any(is.na(membership_no_punish[,g]))
}

# data for punishment condition #
membership_punish <- array(0,c(groupSize,ngroups))
c_punish <- array(0,c(groupSize,ntrials,ngroups)) # choices
Gga_punish <- array(0,c(ntrials,ngroups)) # group-averaged contribution (only 1 entry per group - hence the G + ga)
Ggas_punish <- array(0,c(groupSize,ntrials,ngroups)) # same as Gga, but specified for each subject (cuz that's how the JAGS-code wants it) - hence the s
Gc_punish <- array(0,c(groupSize,ntrials,ngroups)) # summed ("cumulated" hence the c) group contribution not including oneself (therefore specified for each subject) - we don't use this - this refers to the sum-command in the loop-structure
Ga_punish <- array(0,c(groupSize,ntrials,ngroups)) # group-averaged contribution without oneself - we don't use this cuz the participants don't see this (hence, they have to do quite a bit of mental arithmetics to represent this), and thus we're modeling their conditional preference relative to the averaged group contribution (including their own)
member_NA_punish <- array(FALSE,ngroups)

for (g in 1:ngroups) {
  c_punish[,,g] <- rbind(redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="P-experiment"][1:10],
                         redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="P-experiment"][11:20],
                         redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="P-experiment"][21:30],
                         redDat$senderscontribution[redDat$groupid==group_names[g]&redDat$p=="P-experiment"][31:40])
  
  Gga_punish[,g] <- colMeans(c_punish[,,g])
  
  for (s in 1:groupSize) {
    membership_punish[s,g] <- redDat$membership[redDat$groupid==group_names[g]&redDat$p=="P-experiment"][s+(s-1)*10]
    Gc_punish[s,,g] <- colSums(c_punish[-s,,g])
    Ga_punish[s,,g] <- colMeans(c_punish[-s,,g])
    Ggas_punish[s,,g] <- colMeans(c_no_punish[,,g])
  }
  member_NA_punish[g] <- any(is.na(membership_punish[,g]))
}

# compile data from each condition into 4D matrix
c <- array(0,c(groupSize,ntrials,ngroups,2))
c[,,,1] <- c_no_punish
c[,,,2] <- c_punish

Gga <- array(0,c(ntrials,ngroups,2))
Gga[,,1] <- Gga_no_punish
Gga[,,2] <- Gga_punish

Ggas <- array(0,c(groupSize,ntrials,ngroups,2))
Ggas[,,,1] <- Ggas_no_punish
Ggas[,,,2] <- Ggas_punish


Gc <- array(0,c(groupSize,ntrials,ngroups,2))
Gc[,,,1] <- Gc_no_punish
Gc[,,,2] <- Gc_punish

Ga <- array(0,c(groupSize,ntrials,ngroups,2))
Ga[,,,1] <- Ga_no_punish
Ga[,,,2] <- Ga_punish

membership <- array(0,c(groupSize,ngroups,2))
membership[,,1] <- membership_no_punish
membership[,,2] <- membership_punish

c <- c[,,!member_NA_no_punish,]
Gga <- Gga[,!member_NA_no_punish,]
Ggas <- Ggas[,,!member_NA_no_punish,]
Gc <- Gc[,,!member_NA_no_punish,]
Ga <- Ga[,,!member_NA_no_punish,]
membership <- membership[,!member_NA_no_punish,]

c_choice_index <- c

ngroups <- dim(membership)[2]

c_win <- c_no_punish

# calculate the winnings (i.e. apply the multiplication-factor to the sum of each groups contributions)
# winnings per group
winnings <-  array(0, ngroups)
for (g in 1:ngroups) {
  winnings[g] <- sum(colSums(c_win[,,g])*pi)
}

member_winnings <- array(NA, c(groupSize, ngroups))
for (g in 1:ngroups) {
  for (s in 1:groupSize) {
    member_winnings[s,g] <- sum((colSums(c_win[,,g])*pi)/4)-sum(c[s,,g,1])
  }
}

member_c <- array(NA, c(groupSize, ngroups))
for (g in 1:ngroups) {
  for (s in 1:groupSize) {
    member_c[s,g] <- sum(c[s,,g,1])
  }
}

first_c <- array(NA, c(groupSize, ngroups))
for (g in 1:ngroups) {
  for (s in 1:groupSize) {
    first_c[s,g] <- sum(c[s,1,g,1])
  }
}

#################################################################
#------------------ CC model analysis ---------------------------
#################################################################

#-------------------  compare belief weights and slope of prefs in CC model for membership vs. no-membership ---------------

# standardise covariate

Ga_old <- Ga
Ga <- Ggas

data <- list("groupSize", "ngroups", "ntrials","c","Ga","membership") 
params <- c("mu_alpha_log_member","mu_alpha_log","mu_rho_probit_member","mu_rho_probit",
            "mu_omega_probit_member","mu_omega_probit") 

# - run jags code
start_time = Sys.time()
CC.samples <- jags.parallel(data, inits=NULL, params,
                   model.file ="CC_corr_member_split.txt",
                   n.chains=3, n.iter=15000, n.burnin=5000, n.thin=1, n.cluster=4)

end_time = Sys.time()
end_time - start_time


#################################################################
#------------------ Plotting/Bayes Factors-----------------------
#################################################################

#### GROUP AVG: ALPHA for membership (= belief about initial group contribution)

BF_effect <- NULL
BF_null <- NULL

# savage dickey plot
plot(density(rnorm(30000,0,1/sqrt(1))),ylim=c(0,1),main=" ")
lines(density(CC.samples$BUGSoutput$sims.list$mu_alpha_log_member),col="red")

fit.posterior <- logspline(CC.samples$BUGSoutput$sims.list$mu_alpha_log_member)
null.posterior <- dlogspline(0, fit.posterior)

null.prior     <- dnorm(0,0,(1/sqrt(1)))                   

BF_effect$alpha_member <- null.prior/null.posterior
BF_null$alpha_member <- null.posterior/null.prior

# estimating and plotting Bayesian credible interval (95%)
MAP_alpha_member <- MPD(CC.samples$BUGSoutput$sims.list$mu_alpha_log_member)
MAP_alpha_member

alpha_member_cred = quantile(CC.samples$BUGSoutput$sims.list$mu_alpha_log_member,c(0.025,0.975))
alpha_member_cred

# Plot density of posterior + credible interval (95%) incl. MAP (maximum a posteori probability)
plot(density(CC.samples$BUGSoutput$sims.list$mu_alpha_log_member), main=" ")
lines(alpha_member_cred, c(0,0), col="red")
points(MAP_alpha_member, 0, pch=19, col="red")


#### GROUP AVG: ALPHA for no-membership (= belief about initial group contribution)

BF_effect <- NULL
BF_null <- NULL

# savage dickey plot
plot(density(rnorm(30000,0,1/sqrt(1))),ylim=c(0,1),main=" ")
lines(density(CC.samples$BUGSoutput$sims.list$mu_alpha_log),col="red")

fit.posterior <- logspline(CC.samples$BUGSoutput$sims.list$mu_alpha_log)
null.posterior <- dlogspline(0, fit.posterior)

null.prior     <- dnorm(0,0,(1/sqrt(1)))                   

BF_effect$mu_alpha_log <- null.prior/null.posterior
BF_null$mu_alpha_log <- null.posterior/null.prior

# estimating and plotting Bayesian credible interval (95%)
MAP_mu_alpha_log <- MPD(CC.samples$BUGSoutput$sims.list$mu_alpha_log)
MAP_mu_alpha_log

mu_alpha_log_cred = quantile(CC.samples$BUGSoutput$sims.list$mu_alpha_log,c(0.025,0.975))
mu_alpha_log_cred

# Plot density of posterior + credible interval (95%) incl. MAP (maximum a posteori probability)
plot(density(CC.samples$BUGSoutput$sims.list$mu_alpha_log), main=" ")
lines(mu_alpha_log_cred, c(0,0), col="red")
points(MAP_mu_alpha_log, 0, pch=19, col="red")




#### GROUP AVG: RHO for membership (= degree of undermatching, aka. conditional cooperation)

# savage dickey plot
plot(density(rnorm(30000,0,1/sqrt(1))),ylim=c(0,1),main=" ")
lines(density(CC.samples$BUGSoutput$sims.list$mu_rho_probit_member),col="red")

fit.posterior <- logspline(CC.samples$BUGSoutput$sims.list$mu_rho_probit_member)
null.posterior <- dlogspline(0, fit.posterior)

null.prior     <- dnorm(0,0,(1/sqrt(1)))                   

BF_effect$rho_member <- null.prior/null.posterior
BF_null$rho_member <- null.posterior/null.prior

# estimating and plotting Bayesian credible interval (95%)
MAP_rho_member <- MPD(CC.samples$BUGSoutput$sims.list$mu_rho_probit_member)
MAP_rho_member

rho_member_cred = quantile(CC.samples$BUGSoutput$sims.list$mu_rho_probit_member,c(0.025,0.975))
rho_member_cred

# Plot density of posterior + credible interval (95%) incl. MAP (maximum a posteori probability)
plot(density(CC.samples$BUGSoutput$sims.list$mu_rho_probit_member), main=" ")
lines(rho_member_cred, c(0,0), col="red")
points(MAP_rho_member, 0, pch=19, col="red")

#### GROUP AVG: RHO for no-membership (= belief about initial group contribution)

BF_effect <- NULL
BF_null <- NULL

# savage dickey plot
plot(density(rnorm(30000,0,1/sqrt(1))),ylim=c(0,1),main=" ")
lines(density(CC.samples$BUGSoutput$sims.list$mu_rho_probit),col="red")

fit.posterior <- logspline(CC.samples$BUGSoutput$sims.list$mu_rho_probit)
null.posterior <- dlogspline(0, fit.posterior)

null.prior     <- dnorm(0,0,(1/sqrt(1)))                   

BF_effect$mu_rho_probit <- null.prior/null.posterior
BF_null$mu_rho_probit <- null.posterior/null.prior

# estimating and plotting Bayesian credible interval (95%)
MAP_mu_rho_probit <- MPD(CC.samples$BUGSoutput$sims.list$mu_rho_probit)
MAP_mu_rho_probit

mu_rho_probit_cred = quantile(CC.samples$BUGSoutput$sims.list$mu_rho_probit,c(0.025,0.975))
mu_rho_probit_cred

# Plot density of posterior + credible interval (95%) incl. MAP (maximum a posteori probability)
plot(density(CC.samples$BUGSoutput$sims.list$mu_rho_probit), main=" ")
lines(mu_rho_probit_cred, c(0,0), col="red")
points(MAP_mu_rho_probit, 0, pch=19, col="red")




#### GROUP AVG: OMEGA for membership (= attention to others, i.e. how much weight on people's previous contribution compared to initial belief)

# savage dickey plot
plot(density(rnorm(30000,0,1/sqrt(1))),ylim=c(0,1),main=" ")
lines(density(CC.samples$BUGSoutput$sims.list$mu_omega_probit_member),col="red")

fit.posterior <- logspline(CC.samples$BUGSoutput$sims.list$mu_omega_probit_member)
null.posterior <- dlogspline(0, fit.posterior)

null.prior     <- dnorm(0,0,(1/sqrt(1)))                   

BF_effect$omega_member <- null.prior/null.posterior
BF_null$omega_member <- null.posterior/null.prior

# estimating and plotting Bayesian credible interval (95%)
MAP_omega_member <- MPD(CC.samples$BUGSoutput$sims.list$mu_omega_probit_member)
MAP_omega_member

omega_member_cred = quantile(CC.samples$BUGSoutput$sims.list$mu_omega_probit_member,c(0.025,0.975))
omega_member_cred

# Plot density of posterior + credible interval (95%) incl. MAP (maximum a posteori probability)
plot(density(CC.samples$BUGSoutput$sims.list$mu_omega_probit_member), main=" ")
lines(omega_member_cred, c(0,0), col="red")
points(MAP_omega_member, 0, pch=19, col="red")

#### GROUP AVG: OMEGA for no-membership (= belief about initial group contribution)

BF_effect <- NULL
BF_null <- NULL

# savage dickey plot
plot(density(rnorm(30000,0,1/sqrt(1))),ylim=c(0,1),main=" ")
lines(density(CC.samples$BUGSoutput$sims.list$mu_omega_probit),col="red")

fit.posterior <- logspline(CC.samples$BUGSoutput$sims.list$mu_omega_probit)
null.posterior <- dlogspline(0, fit.posterior)

null.prior     <- dnorm(0,0,(1/sqrt(1)))                   

BF_effect$mu_omega_probit <- null.prior/null.posterior
BF_null$mu_omega_probit <- null.posterior/null.prior

# estimating and plotting Bayesian credible interval (95%)
MAP_mu_omega_probit <- MPD(CC.samples$BUGSoutput$sims.list$mu_omega_probit)
MAP_mu_omega_probit

mu_omega_probit_cred = quantile(CC.samples$BUGSoutput$sims.list$mu_omega_probit,c(0.025,0.975))
mu_omega_probit_cred

# Plot density of posterior + credible interval (95%) incl. MAP (maximum a posteori probability)
plot(density(CC.samples$BUGSoutput$sims.list$mu_omega_probit), main=" ")
lines(mu_omega_probit_cred, c(0,0), col="red")
points(MAP_mu_omega_probit, 0, pch=19, col="red")


#################################################################
# Converting output to interpretable values (aka. reparameterizing)
#################################################################

##### GROUP AVG: ALPHA for membership (from normal to gamma distr - via loglink( inverse=TRUE)
# estimating and plotting Bayesian credible interval (95%)
MAP_alpha_member <- MPD(loglink(CC.samples$BUGSoutput$sims.list$mu_alpha_log_member, inverse=TRUE))
MAP_alpha_member

alpha_member_cred = quantile(loglink(CC.samples$BUGSoutput$sims.list$mu_alpha_log_member, inverse=TRUE),c(0.025,0.975))
alpha_member_cred

# Plot density of posterior + credible interval (95%) incl. MAP (maximum a posteori probability)
plot(density(loglink(CC.samples$BUGSoutput$sims.list$mu_alpha_log_member, inverse=TRUE)), 
     main = "Group AVG alpha membership (reparameterized) - belief about initial contrib")
lines(alpha_member_cred, c(0,0), col="red")
points(MAP_alpha_member, 0, pch=19, col="red")

##### GROUP AVG: ALPHA for no-membership (from normal to gamma distr - via loglink( inverse=TRUE)
# estimating and plotting Bayesian credible interval (95%)
MAP_alpha <- MPD(loglink(CC.samples$BUGSoutput$sims.list$mu_alpha_log, inverse=TRUE))
MAP_alpha

alpha_cred = quantile(loglink(CC.samples$BUGSoutput$sims.list$mu_alpha_log, inverse=TRUE),c(0.025,0.975))
alpha_cred

# Plot density of posterior + credible interval (95%) incl. MAP (maximum a posteori probability)
plot(density(loglink(CC.samples$BUGSoutput$sims.list$mu_alpha_log, inverse=TRUE)), 
     main = "Group AVG alpha no-membership (reparameterized) - belief about initial contrib")
lines(alpha_cred, c(0,0), col="red")
points(MAP_alpha, 0, pch=19, col="red")


#### GROUP AVG: RHO fpr membership (= degree of undermatching, aka. conditional cooperation)
# estimating and plotting Bayesian credible interval (95%)
# MAP_rho_member <- MPD(pnorm(CC.samples$BUGSoutput$sims.list$mu_rho_probit_member))
MAP_rho_member <- MPD(probitlink(CC.samples$BUGSoutput$sims.list$mu_rho_probit_member, inverse=TRUE))
MAP_rho_member

rho_member_cred = quantile(probitlink(CC.samples$BUGSoutput$sims.list$mu_rho_probit_member, inverse=TRUE),c(0.025,0.975))
rho_member_cred

# Plot density of posterior + credible interval (95%) incl. MAP (maximum a posteori probability)
plot(density(probitlink(CC.samples$BUGSoutput$sims.list$mu_rho_probit_member, inverse=TRUE)), 
     main = "Group AVG rho membership (reparameterized) - conditional cooperation")
lines(rho_member_cred, c(0,0), col="red")
points(MAP_rho_member, 0, pch=19, col="red")

#### GROUP AVG: RHO no-membership (= degree of undermatching, aka. conditional cooperation)
# estimating and plotting Bayesian credible interval (95%)
MAP_mu_rho_probit <- MPD(probitlink(CC.samples$BUGSoutput$sims.list$mu_rho_probit, inverse=TRUE))
MAP_mu_rho_probit

mu_rho_probit_cred = quantile(probitlink(CC.samples$BUGSoutput$sims.list$mu_rho_probit, inverse=TRUE),c(0.025,0.975))
mu_rho_probit_cred

# Plot density of posterior + credible interval (95%) incl. MAP (maximum a posteori probability)
plot(density(probitlink(CC.samples$BUGSoutput$sims.list$mu_rho_probit, inverse=TRUE)), 
     main = "Group AVG rho no-membership (reparameterized) - conditional cooperation")
lines(mu_rho_probit_cred, c(0,0), col="red")
points(MAP_mu_rho_probit, 0, pch=19, col="red")


#### GROUP AVG: OMEGA membership (= attention to others, i.e. how much weight on people's previous contribution compared to initial belief)
# estimating and plotting Bayesian credible interval (95%)
MAP_omega_member <- MPD(probitlink(CC.samples$BUGSoutput$sims.list$mu_omega_probit_member, inverse=TRUE))
MAP_omega_member

omega_member_cred = quantile(pnorm(CC.samples$BUGSoutput$sims.list$mu_omega_probit_member),c(0.025,0.975))
omega_member_cred

# Plot density of posterior + credible interval (95%) incl. MAP (maximum a posteori probability)
plot(density(probitlink(CC.samples$BUGSoutput$sims.list$mu_omega_probit_member, inverse=TRUE)), 
     main = "Group AVG omega membership (reparameterized) - attention to others")
lines(omega_member_cred, c(0,0), col="red")
points(MAP_omega_member, 0, pch=19, col="red")

#### GROUP AVG: OMEGA no-membership (= degree of undermatching, aka. conditional cooperation)
# estimating and plotting Bayesian credible interval (95%)
MAP_mu_omega_probit <- MPD(probitlink(CC.samples$BUGSoutput$sims.list$mu_omega_probit, inverse=TRUE))
MAP_mu_omega_probit

mu_omega_probit_cred = quantile(probitlink(CC.samples$BUGSoutput$sims.list$mu_omega_probit, inverse=TRUE),c(0.025,0.975))
mu_omega_probit_cred

# Plot density of posterior + credible interval (95%) incl. MAP (maximum a posteori probability)
plot(density(probitlink(CC.samples$BUGSoutput$sims.list$mu_omega_probit, inverse=TRUE)), 
     main = "Group AVG omega no-membership (reparameterized) - attention to others")
lines(mu_omega_probit_cred, c(0,0), col="red")
points(MAP_mu_omega_probit, 0, pch=19, col="red")

