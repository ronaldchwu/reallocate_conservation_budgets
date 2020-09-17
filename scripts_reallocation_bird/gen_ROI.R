####################################################################
#
#  Generate return on investment (ROI) for bird conservation projects
#  - Steps
#     - 1. Identify investment level points to assess their returns
#     - 2. Compute the monetary costs of each investment level
#     - 3. Compute corresponding returns using IUCN transition model (as expected penalty scores based on
#          species' IUCN status at end time)
#           - Transition model estimate per-8-yr (per-unit-time in our case) transitions. So 16-yr transition probabiilties 
#             are computed by applying the transition model twice  
#           - We assume past investment level (e.g. during the 1st 8-yr period) does not affect future transitions
#             (e.g. in the 2nd 8-yr period)
#           
#  - Input:
#     - pp: management strategy considerd (S1, S2, S3, S4)
#     - kk: current reallocation event id
#  - Output:
#     - ROI for each brid species project
#         - each species has a data frame of 
#           'Is', 'Ms', 'Rs' (investment level, monetary costs and return, respectively)
#           and a matrix to store monetary costs at each time step under each investment level
#         - in the bird case, Is = Ms.

gen_ROI <- function(pp, kk, prob_scale){  
  
  # debug

  
  # current time step
  t_now <- T_Rk[kk]
  
  
  # ROIs to return
  # - list of 270 elements, one for each bird species project
  #      - each elment is a list of two
  #           - 1: a data.frame of key ROI information
  #           - 2: a matrix to store monetary costs at each time step under each investment level
  ROIs <- rep(list(list(data.frame("Is"= numeric(n_xi),                    # total investment level till end time
                                   "Ms"= numeric(n_xi),                    # montetary costs 
                                   "Rs"= numeric(n_xi)),                   # return: penality score
                        matrix(0, nrow = n_xi, ncol = K))), np_all)
  
  
  # iterate through species
  for(i in 1:np_all){
    
    #np_all
    
    # debug
    #i <- 109
    
    # check if the bird species has gone extinct or not
    if(IUCN_now[i] != 0){  # if species still present, proceed to compute ROI
      
      # - 1.Identify investment level points to assess their returns
      #    - note: xi is per-time-step monetary investment
      xi <- xi_all[[i]]
      
      # - 2. Compute the monetary costs of each investment level
      Mi <- xi                                                           # investment levels already are monetary
      Mit <- matrix(rep(xi,2), nrow = n_xi, ncol=K, byrow = F)           # extract the monetary costs at each remaining time step
                                                                         # assume allocated budget evenly distribute among time periods
      
      # - 3. Compute corresponding returns using IUCN transition model
      #       - for details see IUCNtransitio.R
      Ri <- mapply(function(x){
        IUCN_end_prob <-IUCNtransition(i, x, t_now, t_end=T_h, prob_scale)  # prob of ending up in [EX, CR, EN, VU, NTLC] status at end time
        return(sum(IUCN_end_prob*penalty_weights))
      }, x = xi) 
      

      # store into ROIs
      #  - need to convert Is/Ms to total budget investment till end of time horizon
      ROIs[[i]][[1]][,"Is"] <- xi_all[[i]] * (K+1-kk)
      ROIs[[i]][[1]][,"Ms"] <- Mi * (K+1-kk)
      ROIs[[i]][[1]][,"Rs"] <- Ri
      ROIs[[i]][[2]] <- Mit
      
    }else if(IUCN_now[i] == 0){    # if this bird species has gone extinct 
      # set its Rs and cost = 0
      # store into ROIs
      ROIs[[i]][[1]][,"Is"] <- xi_all[[i]] * (K+1-kk)
      ROIs[[i]][[1]][,"Ms"] <- 0
      ROIs[[i]][[1]][,"Rs"] <- 0
      ROIs[[i]][[2]] <- matrix(0, nrow = n_xi, ncol = T_h)
    }
    
  } # end of i species loop
  
  return(ROIs)
}
