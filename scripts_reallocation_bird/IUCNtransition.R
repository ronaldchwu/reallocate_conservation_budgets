###############################################################
#
# IUCN transition model: for computing the probability of transitioning between IUCN status based on investment level
#   - Based on historica data of how Australian birds' IUCN status transition correlate with budget investment (see McCarthy et al. 2008)
#   - We follow the orginal studies in assuming:
#       - Status can improve for at most one level at a time (e.g. CR->EN; while CR->VU is not possible)
#           - All one-level improvement are modeled together using a exponential regression (exp(f + g*x))
#       - Status can deterioate for at most two levels at a time; yet NTLC cannot deterioate for more than one level (i.e. only allow NTLC->VU/NTLC)
#           - All one-level deterioation, except 'NTLC->VU', are modeled together using a exponential regression (exp(c + d*x))
#               - 'NTLC->VU' is modelled using an independent exponential regression (exp(y + z*x))
#           - All two-level deterioation are modeled together using a exponential regression (exp(a + b*x))
#
#
#   - Inputs
#       - i: bird species ID
#       - x: per-time-step monetary investment level
#       - t_now: current time step
#       - t_end: target end time step to assess transition probabiilty
#       - prob_scale: scaling factor (x1,x2,....x10) to make state transition more likely
#                     
#   - Outputs
#       - IUCN_end_prob: prob of the bird species ending up in [EX, CR, EN, VU, NTLC] status at end time
#
###############################################################

IUCNtransition <- function(i, x, t_now, t_end, prob_scale){
  
  #debug
  # i <- 1
  # x <- 0
  #t_now <-0
  # t_end <-2
  # IUCN_now[i]
  
  # unpack regression model parameters (regression coefficients)
  #  - from McCarthy et al., 2008
  #  - formula of regression: 

  # use mean estimate of the parameters
  a <- par_mean[1]
  b <- par_mean[2]
  c <- par_mean[3]
  d <- par_mean[4]
  f <- par_mean[5]
  g <- par_mean[6]
  y <- par_mean[7]
  z <- par_mean[8]

  # retrieve species status
  status <- IUCN_now[i]
  
  # variable to store results  
  #   - matrix: row - current status, col - future status
  IUCN_end_prob_all <- matrix(0, nrow = n_status, ncol = n_status)
  
  # Estimate the probability of IUCN status transition
  # if the species is already extinct, it remains extinct
    IUCN_end_prob_all[1,1] <- 1
    
  # if the species is in CR status 
    IUCN_end_prob_all[2,1] <- exp(c + d*x)*prob_scale  # 'CR->EX'
    IUCN_end_prob_all[2,3] <- exp(f + g*x)*prob_scale  # 'CR->EN'       
    sum_trans <- IUCN_end_prob_all[2,1]+IUCN_end_prob_all[2,3]
    if(sum_trans>=1){   # if transition prob out of the current status >=1, adjust them to sum to one
      IUCN_end_prob_all[2,1] <- round(IUCN_end_prob_all[2,1]/sum_trans,3)
      IUCN_end_prob_all[2,3] <- round(IUCN_end_prob_all[2,3]/sum_trans,3)
      IUCN_end_prob_all[2,2] <- 0
    }else{
      IUCN_end_prob_all[2,2] <- 1-sum_trans  #'CR->CR'
    }
    
    
  # if the species is in EN status
    
    IUCN_end_prob_all[3,1] <- exp(a + b*x)*prob_scale  # 'EN->EX'
    IUCN_end_prob_all[3,2] <- exp(c + d*x)*prob_scale  # 'EN->CR'
    IUCN_end_prob_all[3,4] <- exp(f + g*x)*prob_scale  # 'EN->VU'
    sum_trans <- IUCN_end_prob_all[3,1]+IUCN_end_prob_all[3,2]+IUCN_end_prob_all[3,4]
    if(sum_trans>=1){   # if transition prob out of the current status >=1, adjust them to sum to one
      IUCN_end_prob_all[3,1] <- round(IUCN_end_prob_all[3,1]/sum_trans,3)
      IUCN_end_prob_all[3,2] <- round(IUCN_end_prob_all[3,2]/sum_trans,3)
      IUCN_end_prob_all[3,4] <- round(IUCN_end_prob_all[3,4]/sum_trans,3)
      IUCN_end_prob_all[3,3] <- 0
    }else{
      IUCN_end_prob_all[3,3] <- 1- sum_trans  # 'EN->EN'                   
    }
    
    
    
  # if the species is in VU status
    
    IUCN_end_prob_all[4,2] <- exp(a + b*x)*prob_scale  # 'VU->CR'
    IUCN_end_prob_all[4,3] <- exp(c + d*x)*prob_scale  # 'VU->EN'
    IUCN_end_prob_all[4,5] <- exp(f + g*x)*prob_scale  # 'VU->NTLC'
    sum_trans <- IUCN_end_prob_all[4,2]+IUCN_end_prob_all[4,3]+IUCN_end_prob_all[4,5]
    if(sum_trans>=1){   # if transition prob out of the current status >=1, adjust them to sum to one
      IUCN_end_prob_all[4,2] <- round(IUCN_end_prob_all[4,2]/sum_trans,3)
      IUCN_end_prob_all[4,3] <- round(IUCN_end_prob_all[4,3]/sum_trans,3)
      IUCN_end_prob_all[4,5] <- round(IUCN_end_prob_all[4,5]/sum_trans,3)
      IUCN_end_prob_all[4,4] <- 0
    }else{
      IUCN_end_prob_all[4,4] <- 1- sum_trans # 'VU->VU'                   
    }
    
    
    
    
  # if the species is in NTLC status
    
    IUCN_end_prob_all[5,4] <- exp(y + z*x)*prob_scale  # 'NTLC->VU'
    sum_trans <- IUCN_end_prob_all[5,4]
    if(sum_trans>=1){   # if transition prob out of the current status >=1, adjust them to sum to one
      IUCN_end_prob_all[5,4] <- 1
      IUCN_end_prob_all[5,5] <- 0
    }else{
      IUCN_end_prob_all[5,5] <- 1-sum_trans  # 'NTLC->NTLC'
    }
    

  
  # Check how many time periods to assess
  if(t_end-t_now == 1){        # if assessing one-step (8-yr) transition, simply extract the corresponding row
    IUCN_end_prob <- IUCN_end_prob_all[status+1,]
  }else if(t_end-t_now == 2){  # if assessing two-step (16-yr) transition, need to compute probability of transitioning twice 
    IUCN_end_prob <- rowSums(mapply(function(i){
      IUCN_end_prob_all[status+1,i]*IUCN_end_prob_all[i,]
    }, i=c(1:n_status)))
  }else if(t_end-t_now == 3){  # if assessing two-step (16-yr) transition, need to compute probability of transitioning three times 
    temp <- rowSums(mapply(function(i){
      IUCN_end_prob_all[status+1,i]*IUCN_end_prob_all[i,]
    }, i=c(1:n_status)))
    
    IUCN_end_prob <- rowSums(mapply(function(i){
      temp[status+1,i]*IUCN_end_prob_all[i,]
    }, i=c(1:n_status)))
    
  }else if(t_end-t_now == 4){  # if assessing two-step (16-yr) transition, need to compute probability of transitioning four times 
    temp <- rowSums(mapply(function(i){
      IUCN_end_prob_all[status+1,i]*IUCN_end_prob_all[i,]
    }, i=c(1:n_status)))
    
    temp <- rowSums(mapply(function(i){
      temp[status+1,i]*IUCN_end_prob_all[i,]
    }, i=c(1:n_status)))
    
    IUCN_end_prob <- rowSums(mapply(function(i){
      temp[status+1,i]*IUCN_end_prob_all[i,]
    }, i=c(1:n_status)))
  }

  return(IUCN_end_prob)
}



