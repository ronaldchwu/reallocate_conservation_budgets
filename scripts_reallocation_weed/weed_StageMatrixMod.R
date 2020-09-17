#######################################################################################################
#
# Stage matrix model: compute changes in infested/managed/monitoring area during weed infestation
#                     under allocated survey/treatment freqyuency plan from now till end time
#  Type of weed infestation area
#    - A: current area under management (weed plants still found)
#    - Mj: current area having been in monitoring for j time steps (no plants found)
#          - if no plant is found after JJ time steps (monitoring time needed to confirm eradication), 
#            declare eradication and remove the area
#          - anytime plants are found, the area move to Type A and start being actively managed.
#    - H: total infestation area; H = A + sum(Mj)
#
#  Inputs
#    - pid: weed project ID
#    - x: allocated survey efforts (number of survey per time step) to a weed project (vector of t_end-t_now elements for each remaining time steps)
#    - t_now: current time step
#    - t_end: end time step to assess system status
#    - A: current area under management (weed plants still found)
#    - Mj: current area having been in monitoring for j time steps (no plants found)
#    - H: current total infested area   
#    - neighbor_dist: see main script file
#    - He_threshold: threshold of total infested area below which eradication is declared
#
#
#  Outputs
#    - H_end, A_end, Mj_end: expected infestation status by the specified end time
#    - TC: expected total monetary costs till end time
#    - xm: Monetary costs at each time step associated with the allocated survey efforts over time (x)

#                            
StageMatrixMod <- function(pid, x, t_now, t_end,
                           initH,H, A, Mj, neighbor_dist, He_threshold){
  
  # debug
  # pid <- i
  # x <- pf_next[[1]][i]
  # t_now <- t_now
  # t_end <- t_end
  # initH <-dat_weed[i,"H1"]
  # H <- H_now[i]
  # A <- A_now[i]
  # Mj <- Mj_now[[i]]

  # if allocated survey effort is > 0, proceed to estimate management costs
  if(x > 0){
    
    # Parameters of management efforts and effectiveness, based on investment level (number of survey/treatment per year)
    u <- ((10000 * ca[pid]) * mean(x)) / (S[pid]* Wa[pid])  # search effort (hr/ha)
    P_FA <- 1 - (exp(-ca[pid]))^mean(x)                     # proportion of weed plants found in active area
    P_FM <- 1 - (exp(-cm[pid]))^mean(x)                     # proportion of weed plants found in monitored area before they reach maturity
    
    
    # variables to store monetary costs at each time step associated with the allocated survey efforts over time
    cost_t <- rep(0, t_end-t_now)
    
    
    # variables to track weed infestation progress
    #    - H: delimited total infestation area
    #    - A: area actively managed
    #    - Mj: area under monitoring for j time steps after being actively managed (i.e. ensuring no remaining mature plants exist and no seed germination occurs)
    
    Ht <- rep(0, t_end-t_now)  # delineated infestation area (ha)
    At <- rep(0, t_end-t_now)  # area actively managed (ha)
    Mjt <- matrix(0, nrow = JJ[pid], ncol = t_end-t_now)    # matrix of area under monitoring (ha); row = years monitored; col = time step
    Rt <- rep(0, t_end-t_now)   # area reverted back to active (ha); to be added to A_t at next time step
    
    
    
    # discounting factors at each time step, for computing the present value of all future costs 
    #  - convert to monetary value at t = 0  
    discount <- rep(1, t_end-t_now)
    for(tt in 1:(t_end-t_now)){
      discount[tt] <- (1 + d[pid])^(-(tt + t_now - 1))
    }
    
    
    # Deterministic simulation of stage transition under allocated survey efforts (for details, see Dodd et al., 2017, Hester et al., 2013)
    #
    # Step 1: Active area turn into m1 after survey and treatment
    #   - for proportion of area with successful detection and treatment
    # Step 2: mj area is monitored and controlled.
    #   - 1) area can make transition into next stage (m(j+1)) if plant found and killed before maturity
    #        - if 
    #             1) no germination occur (1-G) 
    #             2) newly germinated plant found before maturation and killed (G*P_FM*P_KM)
    #   - 2) area can on the other hand revert back to infested status (R)
    #        - if newly germinated plants are not found & killed before maturity  G * (1- P_FM*P_KM)
    # Step 3: merge R into A ---- to be actively managed at next time step
    #
    # Step 4: update H (= A + all mj)
    #
    
    
    
    # compute weed infestation progress at time t_now+1 based on initial conditions at t_now = 0
    #   - using the input intial status H, A, Mj
    Mjt[1,1] <- (P_FA * P_KA[pid]) * A   # proportion of actively managed area with plants successfully found & killed will advance to M(j=1) stage
    for(j in 2:JJ[pid]){
      Mjt[j,1] <- (1 - G[pid] +  G[pid]* P_FM * P_KM[pid]) * Mj[j-1]  # area advances along Mj classes if they match the conditions in the Step 2.1 description
    }
    Rt[1] <-  G[pid] * (1 - P_FM * P_KM[pid])* sum(Mj)   # Step 2.2: area summed over all Mj classes that revert back to active managed
    At[1] <- A - Mjt[1,1] + Rt[1]
    Ht[1] <- At[1] + sum(Mjt[,1])
    
    
    if(t_end - t_now >= 2){
      # compute weed infestation progress from time t_now+2 to the end of time horizon
      for(tt in 2:(t_end-t_now)){
        
        Mjt[1,tt] <- (P_FA * P_KA[pid]) * At[tt - 1] # proportion of actively managed area with plants successfully found & killed will advance to M(j=1) stage
        for(j in 2:JJ[pid]){
          Mjt[j,tt] <- (1 - G[pid] + G[pid]* P_FM * P_KM[pid]) * Mjt[j-1, tt-1]  # area advances along Mj classes if they match the conditions in the Step 1 description
        }
        Rt[tt] <- (1 - P_FM * P_KM[pid])* G[pid] * sum(Mjt[, tt-1])   # area summed over all Mj classes that revert back to active managed
        At[tt] <- At[tt-1] - Mjt[1,tt] + Rt[tt]
        Ht[tt] <- At[tt] + sum(Mjt[,tt])
        
      } # end of tt loop
    }

    
    # check if infested area fall below the He_threhold at some time step; 
    # if so, eradication will be declared and projec will be terminated at that time step, and there is no more cost after that. 
    t_below_He <- which(Ht < He_threshold)
    if(length(t_below_He > 0)){
      discount[min(t_below_He):(t_end-t_now)] <- 0  # force all future costs = 0 (no more surveys needed) after eradication declaration; for later calculating total monetary costs over time
    }
    
    # sum search cost, control cost, administrative and communication costs over time and compute their present value (with discounting)
    SC <- sum(u*L[pid]*initH*discount)   # total cost of search in present-value; note that the entire landspace (of area size of initial infestation)
                                     # is searched at each time step until eradication is declared
    Mjt_sums <- apply(Mjt, 2, sum)   # compute total mj area in each time step
    if(t_end-t_now >=2){
      CC <- c(((za[pid]*P_FA)*c(A, At[1:(t_end-t_now-1)]) + (zm[pid]*G[pid]*P_FM)* c(sum(Mj), Mjt_sums[1:(t_end-t_now-1)])) %*% discount)   # total cost of control in present-value
    }else{
      CC <- c(((za[pid]*P_FA)*A + (zm[pid]*G[pid]*P_FM)* sum(Mj)) %*% discount)   # total cost of control in present-value
    }
    
    TC <- SC + CC + (AC[pid] + RC[pid]) * sum(discount)  # total cost of eradication management over the time horizon
    
    
    
    # calculate monetary costs at each time step (xm)
    #  - if t_now > 0, xm would be zero before t_now
    xm <- rep(0,T_h)  
    for(tt in 1:(t_end-t_now)){
      if(tt == 1){
        SCt <- u*L[pid]*initH*discount[tt]
        CCt <- (((za[pid]*P_FA)* A) + (zm[pid]*G[pid]*P_FM)* sum(Mj))*discount[tt]
        ARC <- (AC[pid] + RC[pid])*discount[tt]
      }else{
        SCt <- u*L[pid]*initH*discount[tt]
        CCt <- (((za[pid]*P_FA)* At[tt - 1]) + (zm[pid]*G[pid]*P_FM)* Mjt_sums[tt - 1]) * discount[tt]
        ARC <- (AC[pid] + RC[pid])* discount[tt]
      }
      TCt <- SCt + CCt + ARC
      xm[t_now+tt] <- TCt
    }
    
    # prepare outputs
    H_new <- Ht[t_end-t_now]
    A_new <- At[t_end-t_now]
    Mj_new <- Mjt[,t_end-t_now]
  
  }else if(x == 0 | x==-1){   # if allocated survey budget/effort = 0,
                              # or the weed project has not emerged (x=-1) --> TC = 0, it = 0
    # prepare outputs
    H_new <- H
    A_new <- A
    Mj_new <- Mj
    TC <- 0
    xm <- rep(0,T_h)
  }
  
  
  
  outputs <- list(H_new, A_new, Mj_new, TC, xm)
  
  return(outputs)
}