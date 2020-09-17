##################################################################
#
#  Optimize budget allocation among weed eradication projects
#    - Separable programming with piece-wise linear approximation of true ROIs
#    - Steps of constructing a separable programming model
#       - 1. Define alternative actions
#          - possible investment levels of each weed project;
#            each with associated reward (as parameter in objective function)
#            and cost (as parameter in total budget constraint)#       
#       - 2. Define auxilary points along ROIs: 
#            for discrete ROIs, generate points that are
#              a. surrounding the possible surve/treatment levels, with distance = neighboring_distance  (e.g.3-0.0001, 3+0.0001 )
#              b. intercepting points half-way between the possible surve/treatment levels (e.g. 0.5, 1.5, 2.5)
#                  - these points will be forced to be zero, to prohibit selecting investment levels that are not integers
#            - Set the auxilary points as binary decision variables, with SOS constraints
#       - 3. Define objective function
#       - 4. Identify constriants 
#            - Total budget constraint 
#            - All auxilary decision vars of the same weed project must sum to one
#            - SOS constraints on auxilary points 
#            - Intercepting points sum = 0 ---- prohibit invalid investment levels (e.g. freq of survey = 1.5)
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
#    - Inputs
#         - ROIs: current ROIs of each project
#         - B: budget available for direct management 
#         - strategy: management strategy being assessed (S1, S2, S3)
#         - pf_now: ongoing budget allocation portfolio 
#         - success_now: whether a weed species project has succeeded and terminated 
#                        (-1: not emerged yet; 0: no, 1: yes)
#         - kk: current reallocation event id
#    - Outputs
#         - optpf: the optimal budget allocation portfolio for now till the end of time horizon


portfolio_opt <- function(ROIs, B, strategy, pf_now, success_now, kk){
  
  
  library(gurobi) 
  library(tidyverse)
  
  # debug
   # ROIs <- ROIs_now
   # B <- B_avail_now
   # strategy <- set_strategy[[pp]]

  # optimal portoflio to return
  #    - 1st element of the list: allocated survey/treatment frequency (from now to the end of time horizon)
  #                               (-1: the project has not emerged)
  #    - 2nd element of the list: monetary expenses at each time step associated with the survey/treatment plan 
  #                               (matrix, row is species, col is time step)(-1: the project has not emerged)
  #    - 3rd element of the list: expected total weed risk score under the portfolio
  optpf <- list(rep(0, np_all),                        
                matrix(0, nrow = np_all, ncol = T_h),  
                0)                                     

  
  
  model <- list()  # The separable programming model object
    
    
    # - 1. Define alternative actions
    #     - possible investment levels of each weed project;
    #       each with associated reward (as parameter in objective function)
    #       and cost (as parameter in total budget constraint)
    
    # ----> already specified and asessed in ROIs
    
    
    # - 2. Define auxilary points along ROIs: 
    #      for discrete ROIs, generate points that are
    #      a. surrounding the possible surve/treatment levels, with distance = neighboring_distance  (e.g.3-0.0001, 3+0.0001 )
    #      b. intercepting points half-way between the possible surve/treatment levels (e.g. 0.5, 1.5, 2.5)
    #           - these points will be forced to be zero, to prohibit selecting investment levels that are not integers
    #     - Set the auxilary points as binary decision variables, with SOS constraints
    
    # generate auxilary ROI points and incorporate them into augmented ROIs
    auxROIs <- mapply(function(L, i){
      # debug
      # df <- ROIs[[2]]
      #  i=2
      df <- L[[1]]
      aux_Is <- c(df$Is + neighbor_dist, df$Is - neighbor_dist, df$Is[-length(df$Is)]+diff(df$Is)/2)  # last term generate intercepting points
      aux_Ms <- c(df$Ms, df$Ms, rep(999,length(df$Ms)-1))                                             # set costs: 1) neighboring points have idnetical costs as the true possible investment level ; 2) arbitrary cost that is easy to identify (e.g. 999) to intercepting points
      aux_Rs <- c(df$Rs, df$Rs, rep(df$Rs[1],length(df$Ms)-1))                                        # set returns: 1) neighboring points have idnetical returns as the true possible investment level; 2) arbitrary returns to intercepting points
      df2 <- rbind(df, data.frame("Is"=aux_Is, "Ms"=aux_Ms, "Rs"=aux_Rs))
      df2 <- subset(df2, Is >=0)                                                                      # remove point with negative investment level
      df2 <- df2[order(df2$Is),]
      row.names(df2) <- c(1:nrow(df2))
      df2$pid <- i
      df2$varname <- paste0("a",i,"n",row.names(df2)) 
      return(df2)
    }, ROIs, c(1:np_all), SIMPLIFY = F)
    
    auxROIs_all <- do.call("rbind",auxROIs)     # turn ROIs into long form data.frame
    
    
    # remove all projects that have not emerged (Rs = -1)
    auxROIs_all <- subset(auxROIs_all, Rs >=0)    # remove these projects from optimization target list
    candidate_pid <- unique(auxROIs_all$pid)      # remaining weed projects eligible for portfolio optimization
    ncpid<- length(candidate_pid)                 # number of eligible weed projects
    nvar <- nrow(auxROIs_all)                     # number of decision variables: the auxiliary and true investment level points
    
    
    # If strategy = constrained reallocation & it is currently not at the 1st decision event,
    # ongoing projects must remain at their investment level by which they were first allocated efforts to.
    # So before optimization, we need to pre-allocate budget to these projects
    # and exclude them from optimization
    
    if(strategy == "S2" & kk > 1){
      cr_pid <- which((pf_now[[1]]!=-1) & success_now==0)   # weed project id that are to be excluded
      cr_budget <- mapply(function(df,xi) {                 # budget needed to continue the projects' ongoing budget allocation
        xi <- ifelse(xi == -1, 0, xi)    # for project that has not emerged (x=-1), costs are set to zero
        b<-df[[1]][which(df[[1]][,"Is"]==round(xi)),"Ms"]
        return(ifelse(b>=0, b, 0))
      }, df = ROIs, xi = pf_now[[1]])

      auxROIs_all <- subset(auxROIs_all, !pid %in% cr_pid)  # remove these constrained projects from optimization target list
      B<- B - sum(cr_budget)           # deduce such pre-allocation from available budget pool
      candidate_pid <- unique(auxROIs_all$pid)              # remaining weed projects eligible for portfolio optimization
      ncpid<- length(candidate_pid)                         # number of eligible weed projects
      nvar <- nrow(auxROIs_all)                             # number of decision variables 
    }
    # 
    # debug
    # m1 <- pf_now[[2]][,c(6:20)]
    # m2 <- t(cr_budget)[,6:20]
    # ifelse(m1==m2, 0, m1-m2)
    
    
    
    # If there are still projects to optimize (ncpid > 0),
    # proceed to separable programming (SP) optimization using Gurobi
    # Otherwise, skip optimization and return zero vector
    
    if(ncpid > 0){
      model$varnames <- auxROIs_all$varname     # names of decision variables in the SP model
      model$vtype <- "C"                        # type of decision variables: continuous
      model$lb <- 0                             # lower bound of decision variables: zero
      
    # - 3. Define objective function
    model$obj <- auxROIs_all$Rs              # utility value (i.e. reward) associated with each decision variables
                                             #   - reward is the prob of successful eradication by end time X weed risk value of the spcies
    model$modelsense <- 'max'                # this is a maximization problem: 
                                             #   - maximize expected total weed risk values reduction by end time
      
    # - 4. Identify constriants 
    #     - Total budget constraint 
    #     - All auxilary decision vars of the same weed project must sum to one
    #     - SOS constraints on auxilary points 
    #     - Intercepting points sum = 0 ---- prohibit invalid investment levels (e.g. freq of survey = 1.5)
    #
    #     - Data structure of each constraint: a list of three elements
    #         - 1st element: left-hand side of the constraint equation (e.g. costs associated with each decision variables)
    #         - 2nd element: eqaulity of inequlaity (e.g. <=)
    #         - 3rd element: right-hand side of the constraint equation (e.g. available budgets)
    
      # Total budget constraint
      tc <- list(auxROIs_all$Ms, "<=", B) 
      
      # Sum-to-one constraints
      sto <- list(matrix(0,nrow=ncpid, ncol = nvar,byrow = T), rep("=",ncpid), rep(1,ncpid))
      for(i in 1:ncpid){
        sto[[1]][i,which(auxROIs_all$pid==candidate_pid[i])] <- 1   # identify all aux vars corresponding to weed project i
      }
      
      # intercepting points sum to zero
      ipz <- list(ifelse(auxROIs_all$Ms==999,1,0), "=", 0)
      
      
      model$A <- rbind(matrix(tc[[1]], nrow = 1, byrow = T), sto[[1]], ipz[[1]])   # left-hand side of contraint equations
      model$sense <- c(tc[[2]], sto[[2]], ipz[[2]])                                # (in)equaility signs of contraint equations
      model$rhs <- c(tc[[3]],sto[[3]], ipz[[3]])                                   # right-hand side of contraint equations
      
      
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
        mutate(Is2 = round(Is,-log10(neighbor_dist)-1)) %>%
        summarise(optIs = Is2%*%opt, optMs = Ms%*%opt)     # round investment level to a digit left to neighboring_dist
      
      
      optpf[[1]][candidate_pid] <- optout$optIs       # store the true budget allocation to each target project 
      optpf[[3]] <- result$objval                     # store the expected total weed risk values reduction under this optimal portfolio
      
    }else{   # if no projects to optimize (i.e. all ongoing projects are constrained to their ongoing investment level)
      optpf[[1]] <- rep(0, np_all)
      optpf[[3]] <- 0
    }
    
    
    
    # If strategy = constrained reallocation and it is not the 1st-time allocation, 
    # ongoing projects did not participate in optimization.
    # We need to add their cotinuing efforts and expected utility back to optpf
    if(strategy == "S2" &  kk > 1){
      optpf[[1]][cr_pid] <- pf_now[[1]][which((pf_now[[1]]!=-1) & success_now==0)]  # get the final optimal portfolio over all projects
      optpf[[3]] <- sum(mapply(function(roi, pf){    # update expected return from these projects to portoflio utility
        x <- roi[[1]][which(roi[[1]][,"Is"]==round(pf)),"Rs"]
        x <- ifelse(x<0,0,x)
        return(x)
      }, roi = ROIs, pf = optpf[[1]]))
    }
    
    # Compute monetary expenses between each reallocaion events associated with the optimal portfolio 
    temp <- lapply(ROIs, `[[`, 2)
    # from ROIs, find the cost structure assocaited with the optimal investment level
    for(i in 1:np_all){
      lvlid <- which(ROIs[[i]][[1]][,"Is"] == round(optpf[[1]][i])) # which investment level is allocated
      mcst <- temp[[i]][lvlid,]  # montary cost structure at each time step
      mcst <- ifelse(mcst < 0, 0, mcst)
      optpf[[2]][i,] <- mcst
    }
    
    # for weed projects that have not emerged, 
    # set their allocated efforts to -1
    optpf[[1]][which(success_now==-1)] <- -1


  return(optpf)
}  