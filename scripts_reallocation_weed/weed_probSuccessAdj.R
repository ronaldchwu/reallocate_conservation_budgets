##############################################################
#
#  Compute probability of weed eradication success using Hazard rate model
#    - adjust for past management efforts if there has been past management efforts that continue till now
# 
#  Input:
#    - pp: management strategy being assessed
#    - i: weed species id
#    - H: total infestation area at t = 0 for the species
#    - xi: allocated freq of survey/treatment from now to an specified end time
#          (note: if the species has not emerged, xi would be -1)
#    - t_now: time step at which the current continuous management efforts started (freq of survey/treatment>0)
#    - t_end: end time where prob of eradication success is assessed
#
#  Output:
#    - outputR: probability of successsfult weed eradication by the specified end time
#

probSuccessAdj <- function(pp, i, H, xi, t_now, t_end){
  
  # debug
  # i=1
  # H <- dat_weed[i,"H1"]
  # xi <- 3
  # t_now <- 0
  # t_end <- T_h
  
  
  # retrieve history of mgmt efforts
  pf_history <- simu_weed_iter[c(2*np_all+i)+(pp-1)*np_all*3,5:ncol(simu_weed_iter)]        
  mgmt_period <- ifelse(pf_history>0,1,0)      # in-between which reallocation events the species has been actively managed (with efforts > 0)
  past_avg <- numeric()                        # average efforts during the continuous active management till now 
  
  # identify the start of the current continuous active management
  if(sum(mgmt_period>0)==0){  # if there has been no active mgmt so far, set t_s = -1, past_avg = -1
    t_s <- -1
    past_avg <- -1
  }else{  # if there has been active mgmet, find the start of the current continuous active management
    past_es <- numeric()       # past efforts at each period of the continuous active management till now
    kk <- which(T_Rk == t_now) # current reallocation event id
    t_s <- -1
    for(tk in (kk-1):1){
      if(mgmt_period[tk]==1){
        t_s <- T_Rk[tk]
        past_es <- c(past_es, unlist(pf_history[tk]))
      }else{
        break
      }
    }
    past_avg <- mean(past_es)
  }
  
  if(t_now > 0 & t_s >= 0 & xi!=-1){  # if the weed had been managed for some time (xi!=1 means the sp has emerged)
                                        # need adjusted prob of eradication
    
    # compute averaged efforts from t_s to T_h, under given future effort xi
    avg_pf <- (past_avg*(t_now-t_s) + xi*(T_h-t_now))/(T_h-t_s)  # average efforts
    
    # compute 'lower-bound' prlb 
    #  - = the prob of succesful eradication from time = t_s till now under efforts = avg_pf
    #  - note that infestation size at time = t_s is the initial H
    prlb <- HazardRateMod(H = H, cs[i], DP[i], Wa[i], dto[i], PL[i], avg_E = avg_pf, tm[i], se[i], HRMpars, t_start = t_s, t_end = t_now,  new_E = xi, neighbor_dist)
    
    # compute unadjusted prob. of eradication by T_h
    #  - = the prob of successful eradication from time = t_s till T_h under efforts = avg_pf
    #  - note that infestation size at time = t_s is the initial H
    prt <- HazardRateMod(H = H, cs[i], DP[i], Wa[i], dto[i], PL[i], avg_E = avg_pf, tm[i], se[i], HRMpars, t_start = t_s, t_end = t_end ,  new_E = xi, neighbor_dist)
    

    # get adjusted prob. of eradication by T_h
    # - = the prob of 'NOT' eradicated by now but SUCCEED by T_h
    outputR <- (prt-prlb)/(1-prlb)
    
  }else if(xi==-1){  # if the project has not emerged, outputR is set to zero

    outputR <- 0
    
  }else{  # if currernt time is 0, or if t_s == -1
    # it means the weed project has not been managed before or active mgmt has been paused prior to this reallocation event
    # the avg effort would simply be xi
    
    # compute prob. of eradication by T_h
    #   - note that infestation size start from initial H
    outputR <- HazardRateMod(H = H, cs[i], DP[i], Wa[i], dto[i], PL[i], avg_E = xi, tm[i], se[i], HRMpars, t_start = t_now, t_end = t_end, new_E = xi, neighbor_dist)

  }
  
  
  return(outputR)
}