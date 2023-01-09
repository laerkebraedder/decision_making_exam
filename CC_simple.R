CC_member <- function(groupSize, ntrials, ngroups, rate_alpha, shape_alpha,
                      shape1_rho, shape2_rho, shape1_omega, shape2_omega
                      ) {
  
  # arrays to populate for simulation
  Gb <- array(NA,c(groupSize,ntrials,ngroups))
  Ga <- array(NA,c(groupSize,ntrials,ngroups,2))
  p <- array(NA,c(groupSize,ntrials,ngroups))
  c <- array(NA,c(groupSize,ntrials,ngroups,2))
  alpha_s <- array(NA,c(groupSize,ngroups))
  rho_s <- array(NA,c(groupSize,ngroups))
  omega_s <- array(NA,c(groupSize,ngroups))
  
  #-------------------------------------------------------------------
  #-------------------  Individual level model -----------------------
  #-------------------------------------------------------------------
  
  for (g in 1:ngroups) {
    
    for (s in 1:groupSize) {
      #--------------- Model priors ------------------------------------------------------
      
      alpha_s[s,g] <- rgamma(1,shape_alpha,rate_alpha)
      rho_s[s,g] <- rbeta(1,shape1_rho,shape2_rho)
      omega_s[s,g]  <- rbeta(1,shape1_omega,shape2_omega)
      
      #beliefs about others on first trial - gamma-poisson distribution
      Gb[s,1,g] <- rpois(1,alpha_s[s,g])
      
      # modelled preference and first contribution - see below
      p[s,1,g] <- (rho_s[s,g]*Gb[s,1,g])
      c[s,1,g,1] <- rpois(1,p[s,1,g])
      
    }
    
    # updating the actual group contribution
    Ga[,1,g,1] <- mean(c[,1,g,1])
    
    #--------------- Implementation of CC model --------------------------------
    
    for (t in 2:ntrials) {
      
      for (s in 1:groupSize) {
        
        #- Belief about group contribution
        Gb[s,t,g] <- ((1-omega_s[s,g])*(Gb[s,t-1,g]))+(omega_s[s,g]*(Ga[s,t-1,g,1]))
        
        #- Contribution preference, given belief and matching preference rho  
        p[s,t,g] <- rho_s[s,g]*Gb[s,t,g]
        
        #- Contribution as discrete sample from preferences
        c[s,t,g,1] <- rpois(1,p[s,t,g])
      }
      
      # updating the actual group contribution
      Ga[,t,g,1] <- mean(c[,t,g,1])
      
    }
  }
  
  result <- list(c=c,
                 Ga=Ga)
  
  return(result)

}