##################################################################
#
#  Optimize budget allocation among bird conservation projects
#    - Separable programming with piece-wise linear approximation of true ROIs
#    - Steps of constructing a separable programming model
#       - 1. Define alternative actions
#          - possible investment levels of each bird project;
#            each with associated reward (as parameter in objective function)
#            and cost (as parameter in total budget constraint)#       
#       - 2. Define auxilary points along ROIs: 
#            for continuous ROIs, the pre-defined investment levels (xi_all) are considered the auxilary 
#            decision variables
#       - 3. Define objective function
#       - 4. Identify constriants 
#            - Total budget constraint 
#            - All auxilary decision vars of the same project must sum to one
#            - SOS constraints on auxilary points 
#
#            - Data structure of each constraint: a list of three elements
#                - 1st element: left-hand side of the constraint equation (e.g. costs associated with each decision variables)
#                - 2nd element: eqaulity of inequlaity (e.g. <=)
#                - 3rd element: right-hand side of the constraint equation (e.g. available budgets)
#       
#       - 5. Using linear programming solver to find optimal solution (Gurobi)
#
#
#
#
#    - Inputs
#         - ROIs: current ROIs of each project
#         - B_avail_now: budget available for direct management 
#         - strategy: management strategy being assessed (S1, S2, S3, S4)
#         - pf_now: ongoing budget allocation portfolio 
#         - kk: current reallocation event id
#    - Outputs
#         - optpf: the optimal budget allocation portfolio for now till the end of time horizon


portfolio_opt <- function(ROIs, B_avail_now, strategy, pf_now, kk){
  
  # debug
  # ROIs <- ROIs_now
  # strategy <- set_strategy[pp]
  
  # optimal portoflio to return
  #    - 1st element of the list: allocated total monetary budget (from now to the end of time horizon)
  #    - 2nd element of the list: monetary expenses at each 8-yr period (unit time step)
  #                               (matrix, row is species, col is time step)
  #    - 3rd element of the list: expected total utility at end time under the current portfolio
  optpf <- list(rep(0, np_all),                        
                matrix(0, nrow = np_all, ncol = T_h),  
                0)                                     

  
  
  model <- list()  # The separable programming model object
    
    
    # - 1. Define alternative actions
    #     - possible investment levels of each bird project;
    #       each with associated reward (as parameter in objective function)
    #       and cost (as parameter in total budget constraint)
    
    # ----> already specified and asessed in ROIs
    
    
    # - 2. Define auxilary points along ROIs: 
    #      for continuous ROIs, the pre-defined investment levels (xi_all) are considered the auxilary 
    #      decision variables    
    #     - Set the auxilary points as binary decision variables, with SOS constraints
    
    # generate auxilary ROI points and incorporate them into augmented ROIs
    auxROIs <- mapply(function(L, i){
      df <- L[[1]]
      df2 <- df
      df2 <- subset(df2, Is >=0)  
      df2 <- df2[order(df2$Is),]
      row.names(df2) <- c(1:nrow(df2))
      df2$pid <- i
      df2$varname <- paste0("a",i,"n",row.names(df2)) 
      return(df2)
    }, ROIs, c(1:np_all), SIMPLIFY = F)
    
    auxROIs_all <- do.call("rbind",auxROIs)     # turn ROIs into long form data.frame
    
    candidate_pid <- unique(auxROIs_all$pid)      # bird projects eligible for portfolio optimization
    ncpid<- length(candidate_pid)                 # number of eligible bird projects
    nvar <- nrow(auxROIs_all)                     # number of decision variables
    
    
    model$varnames <- auxROIs_all$varname     # names of decision variables in the SP model
    model$vtype <- "C"                        # type of decision variables: continuous
    model$lb <- 0                             # lower bound of decision variables: zero
      
    # - 3. Define objective function
    model$obj <- auxROIs_all$Rs              # expected penality value associated with each decision variables
                                             #   - here it is the prob of ending up in each IUCN status x the penalty weights
    model$modelsense <- 'min'                # this is a minimization problem: 
                                             #   - minimize expected total penalty scores at end time
      
    # - 4. Identify constriants 
    #     - Total budget constraint 
    #     - All auxilary decision vars of the same project must sum to one
    #     - SOS constraints on auxilary points 
    #     - Species with same IUCN status should receive same amount of investment
    #
    #     - Data structure of each constraint: a list of three elements
    #         - 1st element: left-hand side of the constraint equation (e.g. costs associated with each decision variables)
    #         - 2nd element: eqaulity of inequlaity (e.g. <=)
    #         - 3rd element: right-hand side of the constraint equation (e.g. available budgets)
    
      # Total budget constraint
      tc <- list(auxROIs_all$Ms, "<=", B_avail_now) 
      
      # Sum-to-one constraints
      sto <- list(matrix(0,nrow=ncpid, ncol = nvar,byrow = T), rep("=",ncpid), rep(1,ncpid))
      for(i in 1:ncpid){
        sto[[1]][i,which(auxROIs_all$pid==candidate_pid[i])] <- 1   # identify all aux vars corresponding to project i
      }
      
      # Species with same IUCN status receive same amount of investment
      #  - note: the scripts below assume each IUCN status group always have >=2 species
      
      if(iss_constraint == TRUE){   # if force species of same IUCN status to have identical budget allocation
        p_iss <- mapply(function(i){ # record which species are in each IUCN status group (except EX)
          which(IUCN_now==i-1)
        }, i = c(2:n_status))  
        
        n_iss <- sapply(p_iss, function(x) length(x))  # number of species in each IUCN status group
        
        # we need n_iss - n_status number of contraints
        iss <- list(matrix(0,nrow=(sum(n_iss) - n_status+1), ncol = nvar,byrow = T), 
                    rep("=",(sum(n_iss) - n_status+1)), rep(0,(sum(n_iss) - n_status+1)))
        
        # Set all CR species to have equal investment level
        for(i in 1:(n_iss[1]-1)){
          iss[[1]][i, which(auxROIs_all$pid==p_iss[[1]][i])] <- auxROIs_all[which(auxROIs_all$pid==p_iss[[1]][i]),"Ms"]        # 1,2,3...., i-1th species must have equal investment to.. 
          iss[[1]][i, which(auxROIs_all$pid==p_iss[[1]][i+1])] <- -auxROIs_all[which(auxROIs_all$pid==p_iss[[1]][i]),"Ms"]   # 2,3,4,...., ith species
        }
        
        # Set all EN species to have equal investment level
        for(i in 1:(n_iss[2]-1)){
          iss[[1]][i+n_iss[1]-1, which(auxROIs_all$pid==p_iss[[2]][i])] <- auxROIs_all[which(auxROIs_all$pid==p_iss[[2]][i]),"Ms"]        # 1,2,3...., i-1th species must have equal investment to.. 
          iss[[1]][i+n_iss[1]-1, which(auxROIs_all$pid==p_iss[[2]][i+1])] <- -auxROIs_all[which(auxROIs_all$pid==p_iss[[2]][i]),"Ms"]   # 2,3,4,...., ith species
        }
        
        # Set all VU species to have equal investment level
        for(i in 1:(n_iss[3]-1)){
          iss[[1]][i+n_iss[1]+n_iss[2]-2, which(auxROIs_all$pid==p_iss[[3]][i])] <- auxROIs_all[which(auxROIs_all$pid==p_iss[[3]][i]),"Ms"]       # 1,2,3...., i-1th species must have equal investment to.. 
          iss[[1]][i+n_iss[1]+n_iss[2]-2, which(auxROIs_all$pid==p_iss[[3]][i+1])] <- -auxROIs_all[which(auxROIs_all$pid==p_iss[[3]][i]),"Ms"]   # 2,3,4,...., ith species
        }
        
        # Set all NTLC species to have equal investment level
        for(i in 1:(n_iss[4]-1)){
          iss[[1]][i+n_iss[1]+n_iss[2]+n_iss[3]-3, which(auxROIs_all$pid==p_iss[[4]][i])] <- auxROIs_all[which(auxROIs_all$pid==p_iss[[4]][i]),"Ms"]        # 1,2,3...., i-1th species must have equal investment to.. 
          iss[[1]][i+n_iss[1]+n_iss[2]+n_iss[3]-3, which(auxROIs_all$pid==p_iss[[4]][i+1])] <- -auxROIs_all[which(auxROIs_all$pid==p_iss[[4]][i]),"Ms"]   # 2,3,4,...., ith species
        }
        
        
      }else if(iss_constraint == FALSE){    # if allows species in a IUCN status group to differ in budget allocation,
                   # set iss to be empty matrices
        iss<-list()
        iss[[1]] <- matrix(0,nrow = 0,ncol=nvar)
        iss[[2]] <- character()
        iss[[3]] <- numeric()
        
      }
      
      
      
    
      model$A <- rbind(matrix(tc[[1]], nrow = 1, byrow = T), sto[[1]], iss[[1]])   # left-hand side of contraint equations
      model$sense <- c(tc[[2]], sto[[2]], iss[[2]])                                # (in)equaility signs of contraint equations
      model$rhs <- c(tc[[3]],sto[[3]], iss[[3]])                                   # right-hand side of contraint equations
      
      
      # Generate a set of SOS constraint for each species project
      sos_all <- list()
      for(i in 1:ncpid){
        sos_i <- list()
        sos_i$type <- 2
        sos_i$index <- which(auxROIs_all$pid==candidate_pid[i])
        sos_i$weight <- rep(1,length(sos_i$index))
        sos_all[[i]] <- sos_i 
      }
      
      model$sos <- sos_all    # assign SOS constraints to the model
       
      
      # debug
      
      
   # - 5. Using linear programming solver to find optimal solution (Gurobi)
    #
      params <- list(OutputFlag=0)      # prevent txt output in the console 
      result <- gurobi(model, params)   # solve the optimization problem
      
      # debug
      #result <- gurobi(model)
      
      
   # - 6. Retrieve all optimal decisions and store into optpf
   #    - the output from optimization model is for auxiliary + true investment level decision variables,
   #      need to infer the corresponding true budget allocation to each target project 
      
      auxROIs_all$opt <- result$x
      
      optout <- auxROIs_all %>%
        group_by(pid) %>%
        summarise(optIs = Ms%*%opt, optMs = Ms%*%opt)     
      
      optpf[[1]][candidate_pid] <- optout$optIs       # store the true budget allocation to each target project 
      optpf[[3]] <- result$objval                     # store the expected total penality scores under this optimal portfolio
      
      
    # Compute monetary expenses between each reallocaion events associated with the optimal portfolio 
    #  - directly infer from optpf[[1]] as it is already monetary
    mcst_future <- matrix(rep(optpf[[1]]/(K+1-kk),K+1-kk),  # montary cost structure in-between each future reallocation events
                          ncol=K+1-kk,byrow = F)
                                                            # = total budget / remaining number of reallocation event
    mcst_past <- matrix(0, nrow=np_all, ncol=kk-1)          # fill in zero to columns of past time steps
    mcst <- as.matrix(cbind(mcst_past, mcst_future))
    optpf[[2]] <- mcst


  return(optpf)
}  