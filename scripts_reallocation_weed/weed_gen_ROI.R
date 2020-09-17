####################################################################
#
#  Generate return on investment (ROI) for weed eradication projects
#  - Steps
#     - 1. Identify investment level points to assess their returns
#     - 2. Compute the monetary costs of each investment level
#     - 3. Compute corresponding returns using hazard rate model (as probability of successful eradication 
#          by the end of time horizon X species' risk score)
#           - If a weed project has just started, efforts input to HRM is the allocated survey/treatment frequency
#           - If a weed has been actively managed for some time, need to consider previous efforts and the fact that
#             the weed is still persistent now, to compute an adjusted prob. of eradication by terminal time
#               - First, identify when the ongoing continuous management started (t_s),
#                 compute an time-averged efforts from t_s to now to terminal time (avg_efforts)
#               - Then, use HRM to estimate the prob. of successfult eradication by NOW (t_now) under avg_efforts starting from t_s;
#                 the prob is called 'lower-bound' (prlb)
#               - Then, use HRM to estimate the prob. of successfult eradication by terminal time under avg_efforts starting from t_s;
#                 the prob is called 'unadjusted terminal prob' (prt)
#               - Finally, we get the adjusted prob. of eradication by terminal time, prt_adj = (prt-prlb)/(1-prlb)
#           
#  - Input:
#     - pp: management strategy considerd (S1, S2, S3)
#     - kk: current reallocation event id
#  - Output:
#     - ROI for each weed species project
#         - each species has a data frame of 
#           'Is', 'Ms', 'Rs' (investment level, monetary costs and return, respectively)
#           and a matrix to store monetary costs at each time step under each investment level
#

gen_ROI <- function(pp, kk){  
  
  # debug

  
  # current time step
  t_now <- T_Rk[kk]
  
  
  # ROIs to return
  # - list of 50 elements, one for each weed species project
  #      - each elment is a list of two
  #           - 1: a data.frame of key ROI information
  #           - 2: a matrix to store monetary costs at each time step under each investment level
  ROIs <- rep(list(list(data.frame("Is"= seq(min_xi, max_xi, by = bin_xi),  # raw investment level
                                   "Ms"= 0,                                 # montetary costs 
                                   "Rs"= 0),                                # return
                        matrix(0, nrow = n_xi, ncol = K))), np_all)
  
  
  # iterate through species
  for(i in 1:np_all){
    
    # debug
    #i <- 1
    
    # check if species project has emerged and is persisent now (successs = 0)
    if(success_now[i] == 0){
      
      # - 1.Identify investment level points to assess their returns
      xi <- xi_all[[i]]
      
      # - 2. Compute the monetary costs of each investment level using stage-matrix model
      #   - use information of current project status
      out <- mapply(function(xi){
        temp<-StageMatrixMod(i, xi, t_now, T_h, dat_weed[i,"H1"], H_now[i], A_now[i], Mj_now[[i]], neighbor_dist, He_threshold)
        return(list(temp[[4]],temp[[5]]))
      }, xi = xi_all[[i]], SIMPLIFY = F)
      
      Mi <- unlist(lapply(out, `[[`, 1))                # extract the total monetary costs
      Mit <- do.call("rbind", lapply(out, `[[`, 2))     # extract the monetary costs at each time step
      
      # - 3. Compute corresponding returns using hazard rate model (as probability of successful eradication at terminal times X species' risk score)
      # check if the weed has been managed for sometime
      #     - if so, compute adjusted prob. of eradicatio based on how long and how much efforts have been invested
      #          - for details see probSuccessAdj.R
      Ri <- mapply(function(xi){
        outputR <- probSuccessAdj(pp, i, dat_weed[i,"H1"], xi, t_now, T_h)
        Ri_temp <- outputR*dat_weed[i,"Risk"]
      }, xi = xi_all[[i]])
      
      
      # store into ROIs
      ROIs[[i]][[1]][,"Is"] <- xi_all[[i]]
      ROIs[[i]][[1]][,"Ms"] <- Mi
      ROIs[[i]][[1]][,"Rs"] <- Ri
      ROIs[[i]][[2]] <- Mit
      
    }else if(success_now[i]==-1){   # if this weed species has not emerged, set its Rs and Ms to -1
      
      # store into ROIs
      ROIs[[i]][[1]][,"Is"] <- xi_all[[i]]
      ROIs[[i]][[1]][,"Ms"] <- -1
      ROIs[[i]][[1]][,"Rs"] <- -1
      ROIs[[i]][[2]] <- matrix(-1, nrow = n_xi, ncol = T_h)
    }else if(success_now[i]==1){    # if this weed species has been successfully eradicated, 
      # set its Rs to zero and cost = 0
      # store into ROIs
      ROIs[[i]][[1]][,"Is"] <- xi_all[[i]]
      ROIs[[i]][[1]][,"Ms"] <- 0
      ROIs[[i]][[1]][,"Rs"] <- 0
      ROIs[[i]][[2]] <- matrix(0, nrow = n_xi, ncol = T_h)
    }
  } # end of i species loop
  
  return(ROIs)
}
