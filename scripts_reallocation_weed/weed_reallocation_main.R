#################################################
#
# Budget reallocation analysis: Invasive weeds in Victoria, Australia
#

# Chung-Huey Wu (chung.huey.wu@gmail.com)



# Set library path to folder with MDPtoolbox (and linprog)
.libPaths("/home/chunghueyw/CH2_reallocation/lib")


# Load packages
library(tidyverse)
require(tictoc)
require(doSNOW)
require(tidyverse)
require(gurobi)


# load commnad-line argument (in Spartan .slurm submission)
#  - Arg1: random seed number (1-500)
#  - Arg2: reallocation frequency (every 1,2,3,...10 yr)
#  - Arg3: expected number of new project per year (0.25,0.5,1,1.5,2)


command_args <- commandArgs(trailingOnly = TRUE)
freq <- as.numeric(command_args[1]) 
expected_nppyr <- as.numeric(command_args[2]) 


# DEBUG
# freq <- 1
# expected_nppyr <- 1



# Load functions
source("./R/main/weed/weed_IncursionScenario.R")
source("./R/main/weed/weed_StageMatrixMod.R")
source("./R/main/weed/weed_HazardRateMod.R")
source("./R/main/weed/weed_gen_ROI.R")
source("./R/main/weed/weed_probSuccessAdj.R")
source("./R/main/weed/weed_portfolio_opt.R")




#------------------------------------------------------------------------------
# 1) Import data of candidate weed projects -------------------------

# data of 50 candidate weed species
dat_weed <- read.csv("./data/AU_weed.csv")

# assign which weed species (species ID) are initially available at t = 0
p_t0 <- c(seq(1, 39, by = 2), 41:50)  # arbitrarily select 30 species
p_new <- which(!c(1:50) %in% p_t0)    # other 20 species that may emerge later

np_t0 <- length(p_t0)                 # number of initial species
np_new <- length(p_new)               # number of late-emerging species
np_all <- np_t0 + np_new              # total number of species




#------------------------------------------------------------------------------
# 2) Setup decision contexts, reassessement scheme, and other parameters ------------------

# Time horizon and reassessment-reallocation schedule
# - note: time start from t = 0
T_h <- 20                                  # time horizon (year) 
T_unit <- 1                                # time step unit (year)
n_T <- T_h/T_unit                          # number of time steps over the time horizon

T_Rk <- seq(0,T_h,by=freq)                 # which time steps to make budget reallocation decisions
T_Rk <- T_Rk[T_Rk<T_h]

T_Aj <- c(T_Rk, T_h)                       # which time steps to assess weed eradication progress


K <- length(T_Rk)                          # number of reallocation events
J <- length(T_Aj)                          # number of reassessment events


# Speed of new weed incursion
lambda <- expected_nppyr                  # average incursion per year, as parameter of a Poisson process


# Management budget, and costs of reassessment & reallocation
#   - throughout the analyses, all monetary values are discounted and converted to value at t = 0
B_unit <- 1000000        # budget unit (A$ 1 million)
B_all <- 29.181442       # budget levels used in our analysis; similar to 
# real world budget (see Dodd et al. 2017)

Ca <- 0                  # cost of assessing weed eradication progress and 
# computing ROI (for ongoing and new weed projects)(per species project)

Cr <- 0                  # cost of conducting budget reallocation analysis and decision-making (for entire portfolio)


# Total budget available for direct management of weed projects
B_avail <- B_all - Ca*np_t0*(J-1) - Cr*K      # B_all minus preallocated reassessment & reallocation budgets 

# Management strategies with varying flexibility in reallocation
set_strategy <- c("S1", "S2", "S2")
n_strategy <- length(set_strategy)


# Project Return on Investment (ROI):  
#  - Investment: management efforts allocated to a project till the end of time horizon
#                assumed to be consistent across time steps
#  - Return: probability of successsful eradication by end time X weed risk value of the species
min_xi <- 0         # minimum freq of survey/treatment allocated to a weed project
max_xi <- 3         # maximum freq of survey/treatment allocated to a weed project
bin_xi <- 1         # bin size of investment level
xi <- seq(min_xi, max_xi, by = bin_xi)   # all possible freq of survey/treatment
n_xi <- length(xi)                       # number of possible freq of survey/treatment
xi_all <- rep(list(xi), np_all)          # possible freq of survey/treatment for each weed species
# here assumed to be all the same
neighbor_dist <- 0.001                   # for separable programing, some auxilary points 
# surronding the true xi will be generated (e.g. 3-0.001,
# 3+0.001). neighbor_dist sets the distance of these 
# auxilary points from the true xi.


# Forward simulation
n_iter <- 50    # number of weed management simulations to run






#------------------------------------------------------------------------------
# 3) Prepare weed-specific system models ------------------

# - 3.1) Construct new weed incursion scenarios via simulation
#  - Using Poisson process to generate scenarios of new weed emergence at each time step over the time horizon
#  - To be used in forward simulation of weed management progress
#  - ----> see StageMatrixMod.R for script details        

pool <- IncursionScenario(p_new, np_new, lambda, T_h, n_iter)   # pool of n=n_iter scenearios of weed emergence      
scenario_pool <- pool[1]                                        # scnearios of which weed species emerged at each time period t 
# (i.e. between time steps t-1 and t)

# - 3.2) Stage-matrix model for estimating weed eradication progress and costs (Hester et.al., 2013; Dodd et.al., 2017)
#  - Model changes in total infestation area, area under active management, and area under monitoring (before declaring eradication)
#    based on characteristics of weed species and the freq of survey/treatment applied
#  - We use this model to estimate the monetary costs of eradicating the species with a given survey/treatment frequency from a starting 
#    time step till an end time. Such cost informatio is used for ROI estimation.
#  - The model is also used in forward simulation of weed project progress: for computing size of infestation area over time
#  - ----> see StageMatrixMod.R for script details        

# The parameters for the stage-matrix model (derived from Dodd et al., 2017)
#  - vectors of characteristics of all 50 weed species projects 
L <- dat_weed$L / B_unit    # cost of labour for searching weed ($/hr)
d <- dat_weed$d             # discount rate (0.06 for all species)
z <- dat_weed$z / B_unit    # control cost ($/ha)
za <- z                     # control cost in active area ($/ha)
zm <- z                     # control cost in monitoring area ($/ha)
AC <- dat_weed$AC / B_unit  # administrative cost ($/yr)
RC <- dat_weed$RC / B_unit  # commnunication cost ($/yr)
dto <- dat_weed$dto         # distance to nearest mgmt office (km)
cs <- dat_weed$cs           # climate suitability score
DP <- dat_weed$DP           # annual period of detectability prior to seed set (months)
tm <- dat_weed$tm           # tm: time to reproductive maturity (years)
se <- dat_weed$se           # se: whether the species has been eradicated elsewhere (1 when yes, 0 when no)
S <- dat_weed$S             # speed of search (meter/hr)
PL <- dat_weed$PL           # propagule longevity (yr)
ap <- dat_weed$ap           # additional monitoring period (yr)
JJ <- PL + ap               # monitoring time needed to confirm eradication
HL <- JJ * log(2)           # seed half-life
G <- -log(0.5)/HL           # seed germination rate
Wa <- dat_weed$Wa           # detection distance (meter) of the searcher in active area
dr <- dat_weed$dr           # detection ratio between Wa/Wm
Wm <- Wa/dr                 # detection distance (meter) of the searcher in monitoring area
cm <- (1/dr)*Wm             # search coverage in monitoring area
ca <- cm*dr                 # search coverage in active area; assume the same survey effort is implemented in all parcels of the landscape

P_KA <- dat_weed$P_KA       # proportion of detected plants that are killed after control is applied in active area
P_KM <- dat_weed$P_KM       # proportion of detected plants that are killed after control is applied in monitored area
He_threshold <- 1           # threshold of infested area below which eradication is declared


# - 3.3) Hazard rate model for estimating probability of successful eradication over time (Dodd et al., 2015)
#   - Regression model on past real-world data of how weed characteristics and management history affects eradication success
#   - We use the model to estiamte the probability of successsful eradication a weed species by an end time given how much efforts have 
#     been/will be invested in the species
#   - It is used in ROI estimation: for calculating the 'return' of the species'.
#   - ----> see HazardRateMod.R for script details

# regression coefficients sourced from Dodd et al. (2015)    
#   - to be used for predicting probability of eradication success
#   - for details about each regression variables, see HazardRateMod.R
HRMpars <- c(-17.11, -0.2299, 0.05092, 0.04962, 1.256, 0.5, -0.1823, 0.83, -2.07, 3.579, 0, 3.548)    




#------------------------------------------------------------------------------
# 4) Iterative reallocation analysis under forward simulation
#    - Implement our proposed reallocation approach
#    - Multiple management strategies with varying reallocation flexibility are parallelly examined in each simulation iteration: 
#       - S1-"Static": no budget reallocation
#       - S2-"Constrained reallocation": only allow budget released from terminated projects (eradication success) to 
#                                     be invested in newly emerged projects
#       - S3-"Full reallocation": budget can be reallocated among weed projects with no constriants
#       - The three strategies are simulated using the same set of weed incurison scenarios (see above), and using the same random seeds of
#         species eradication success (see below). So their conservation outcomes can be directly compared within each simulation iteration.
#
#    - Steps
#      - 4.1. Generate ROI for all exisiting weed projects
#           - use stage-matrix and hazard rate models to compute costs over time and expected total weed risk values at end time, respectively 
#      - 4.2. Find optimal budget allocation among budgets and switch to such new budget portflio
#           - use separable programming
#           - for S1-"Static" strategy: continue using the ongoing portfolio; 
#                                    remaining budget of terminated projects won't be available to other projects
#           - for S2-"Constrained reallocation": separable programming only optimize the allocation of released budgets 
#                                                to newly emerged weed projects  
#           - for S3-"Full reallocation": adopt the new optimal porfolio
#      - 4.3. Simulate weed eradication progress till the next reallocation event
#           - use hazard rate model to generate stochastic eradication success
#              - if a weed is still persist, use stage-matrix model to compute an updated infestation area size
#              - if a weed is successfully eradicated, infestation area are set to zero. The project no longer require investment.
#      - 4.4. Identify newly emerged weed species by the next reallocation event, add to the list of existing weed projects
#           - based on pre-simulated incursion scenarios
#
#      - Repeat 4.1-4.4 until the end of time horizon

# Setup parallel computing
cl <- makeCluster(4)
registerDoSNOW(cl)


# progress tracking in R
pb_P <- txtProgressBar(max = n_iter, style = 3)
progress_P <- function(n) setTxtProgressBar(pb_P, n)
opts_P <- list(progress = progress_P)


# Iterations over weed incursion scenarios
simu_weed_out <- foreach(s = 1:n_iter,
                         .combine = "rbind",
                         .errorhandling = "remove",
                         .options.snow = opts_P,
                         .packages = c("tidyverse", "gurobi"))%dopar%{
  
  # debug
  # s <- 1
  
  # Select new weed incursion scenario for this simulation iteration
  scenario <- scenario_pool[[1]][[s]]
  
  # Data frame to store weed projects status and management history over time in this scenario iteration
  #  - success: whether a species has not emerged (-1), is persistent (0), or successfully eradicated (1) at each assessment time step (t5-t20)
  #  - m_spent: cumulative monetary budget spent on each project by each assessment time step
  #  - efforts: frequency of survey/treatment implemented in a project at each time perid (e.g. value at t10 = freq used during t=5 to t=10) 
  simu_weed_iter <- data.frame(iter = rep(s, np_all*3*n_strategy),
                               spID = rep(c(1:np_all), 3*n_strategy),
                               metrics = rep(c("success", "m_spent", "efforts"), each = np_all),
                               strategy = rep(set_strategy, each = np_all*3))
  for(jj in 1:K){
    simu_weed_iter[,paste0("t",T_Aj[jj+1])] <- -1
  }
  
  
  
  # Generate random seeds of species' eradication success in this simulation iteration
  #  - a uniform random numbers between 0 and 1 for each species at each reassessment/reallcation time step. 
  #  - the seed will be shared by system trajectories under all management strategies in this simulation iteration
  #  - example: If the seed for Species 1 at t=5 is 0.5, and if under strategy S1/S2/S3 the prob of eradicating Species 1 
  #    at t=5 are estimated (using hazard rate model) to be 0.3, 0.6, 0.7 respectively, then only under S2 and S3 trajectories Species 1
  #    is successfully eradicated  (as 0.6/0.7 > 0.5).
  
  seed_success <- matrix(runif(np_all*K), nrow = np_all, ncol = K)
  
  
  
  # Parallel simulation of weed eradication progress under each mgmt strategy (S1:no reallocation, S2: constrained reallocation, S3: full reallocation)
  for(pp in 1:n_strategy){
    
    # debug
    # pp <- 1
    
    # select management strategy
    strategy <- set_strategy[pp]
    
    # - 4.0) Set up initial conditions of weed projects at t = 0, and initialize variables that tracks weed project progress
    
    # Variables that track management budget and current time step
    B_avail_now <- B_avail
    t_now <- 0
    
    # Variables that track project status and management efforts of each of the 50 weed species
    p_now <- p_t0                      # exsiting weed species projects 
    np_now <- length(p_now)            # number of existing weed species projects 
    H_now <- dat_weed[, "H1"]          # total infestation area
    A_now <- H_now                     # infested arae under active management;
    # projects always begin with whole area under active management
    Mj_now <- list()                   # area with no existing adult plants, under monitoring for j time steps (year)
    for(i in 1:np_all){                # - a list of vectors
      Mj_now[[i]] <-  rep(0, JJ[i])
    }
    
    success_now <- rep(-1, np_all)     # whether a weed species project has succeeded and terminated (-1: not emerged yet; 0: no, 1: yes)
    success_now[p_now] <- 0            # initial species projects are set to zero
    
    
    # Variable to track the current budget allocation
    #    - 1st element of the list: allocated survey/treatment frequency (from now to the end of time horizon)
    #                               (-1: the project has not emerged)
    #    - 2nd element of the list: monetary expenses at each time step associated with the survey/treatment plan 
    #                               (matrix, row is species, col is time step)(-1: the project has not emerged)
    #    - 3rd element of the list: expected total weed risk score under the portfolio
    pf_now <- list(rep(-1, np_all),     
                   matrix(-1, nrow = np_all, ncol = T_h),     
                   0)                  
    
    # Variable that tracks cumulative monetary expenditure in each weed project (time discounted to value at time = 0)
    spent_now <- rep(0, np_all)        
    
    
    # Iterate through reallocation events
    for(kk in 1:K){
      
      # debug
      # kk <- 1
      # kk <- kk +1
      
      # Identify which time step it is now
      t_now <- T_Rk[kk]
      rJ <- J-kk    # remaining number of reassessment events
      rK <- K-kk    # remaining number of reallocation events
      
      
      # - 4.1. Generate ROI for all exisiting weed projects
      #     - use stage-matrix and hazard rate models to compute costs over time and expected total weed risk values at end time, respectively 
      #     - if a weed has been successfully eradicated, rewards in are ROI set to zero
      
      ROIs_now <- gen_ROI(pp, kk)
      
      # check distribution of project costs and returns
      ROIs_all <- do.call("rbind", lapply(ROIs_now, "[[", 1))
      
      
      # - 4.2. Find optimal budget allocation among budgets and switch to such new budget portflio
      #     - use separable programming
      #     - for S1-"Static" strategy: continue using the ongoing portfolio; 
      #                                 remaining budget of terminated projects won't be available to other projects
      #     - for S2-"Constrained reallocation": separable programming only optimize the allocation of released budgets 
      #                                          to newly emerged weed projects  
      #     - for S3-"Full reallocation": adopt the new optimal porfolio
      
      if(set_strategy[pp] == "S1" & kk > 1){     # for S1-static strategy, keep using pf_now
        pf_next <- pf_now
      }else{    # proceed to separable programming optimization if it is first allocation event, or for stratgies S2 and S3 
        pf_next <- portfolio_opt(ROIs_now, B_avail_now, set_strategy[pp], pf_now, success_now, kk)  
      }
      
      # - 4.3. Simulate weed eradication progress till the next reallocation event
      #     - use hazard rate model to generate stochastic eradication success
      #        - if a weed is still persist, use stage-matrix model to compute an updated infestation area size
      #        - if a weed is successfully eradicated, infestation area are set to zero. The project no longer require investment.
      
      # Spend money: calculate the monetary costs of implementing the portfolio till the next reallocation event, deduce from budget pool
      # Find the time steps where costs are to be summed
      if(kk < length(T_Rk)){
        ts <- c((T_Rk[kk]+1):T_Rk[kk+1]) 
      }else{
        ts <- c((T_Rk[kk]+1):T_h)
      }     
      
      # Find the weed speies projects that are still ongoing (success = 0) from now to next reallocation event,      
      # for others (e.g. terminated projects), set costs = 0
      m_spent <- pf_next[[2]][,ts, drop = FALSE]
      m_spent[which(success_now!=0),] <- 0
      
      # Deduce sum of these costs from budget pool 
      B_avail_next <- B_avail_now - sum(m_spent)
      spent_next <- spent_now + rowSums(m_spent)  # record the expenditure
      
      # Use hazard rate model to generate stochastic eradication success by the next reallocation time step
      # set the target end time to be next realllocation event for hazard rate model
      if(kk < length(T_Rk)){
        t_end <- T_Rk[kk+1]
      }else{
        t_end <- T_h
      }
      
      # apply hazard rate model to compute probability of success eradication by the next realllocation event
      success_prob <- mapply(function(i, H, xi)probSuccessAdj(pp, i, H, round(xi), t_now, t_end),
                             i = c(1:np_all), H = H_now, xi = pf_next[[1]])
      
      # compare the probabilities against the pre-generated random seeds to determine eradication success
      # (-1: not emerged yet; 0: persistent; 1: eradication succeeded)
      success_next <- mapply(function(prob, seed, s_now){
        s_next <- ifelse(prob >= seed, 1, s_now)
        return(s_next)
      },prob = success_prob, seed = seed_success[,kk], s_now = success_now)
      
      
      # update size of infestation arae for each weed species under implemented survey/treatment till next reallocation event
      H_next <- H_now       # total infestation arae
      A_next <- A_now       # infested arae under active management
      Mj_next <- Mj_now     # area with no existing adult plants, under monitoring for j time steps (year)
      for(i in 1:np_all){
        #debug
        #i <- 1
        if(success_next[i] == 1){           # if weed successfully eradicated, set infestation area to zero
          H_next[i] <- 0                 
          A_next[i] <- 0                 
          Mj_now[[i]] <-  rep(0, JJ[i])
        }else{                              # if weed still persist, use stage-matrix model to estimate infestation progress
          new_states <- StageMatrixMod(i, pf_next[[1]][i], t_now, t_end, dat_weed[i,"H1"],
                                       H_now[i], A_now[i], Mj_now[[i]], neighbor_dist, He_threshold)
          H_next[i] <- new_states[[1]] 
          A_next[i] <- new_states[[2]]   
          Mj_next[[i]] <- new_states[[3]] 
        }
      } # end of i loop
      
      
      # - 4.4. Identify newly emerged weed species by the next reallocation event, add to the list of existing weed projects
      #    - based on pre-simulated incursion scenarios  
      p_to_add <- numeric()  # ID of newly emerged weed species
      
      if(kk < length(T_Rk)){  # before the last reallocation event
        for(tt in (T_Rk[kk]+1):T_Rk[kk+1]){
          p_to_add <- c(p_to_add, scenario[[tt]])
        }  
      }else{                  # at the last reallocation event
        for(tt in (T_Rk[kk]+1):T_h){
          p_to_add <- c(p_to_add, scenario[[tt]])
        }  
      }
      
      p_next <- c(p_now, p_to_add)   # updated list of existing weed species projects; now include newly emerged ones
      np_next <- length(p_next)   
      
      # update success_next: emerged projects change from -1 to 0
      success_next[p_to_add] <- 0
      
      # pre-allocate budget for reassessing the newly emerged projects till the end of time horizon
      B_avail_next <- B_avail_next - Ca*(length(p_to_add))*(rJ)   
      
      
      
      # Update state variables
      B_avail_now <- B_avail_next 
      spent_now <- spent_next 
      pf_now <- pf_next    
      p_now <- p_next
      np_now <- np_next         
      H_now <- H_next
      A_now <- A_next
      Mj_now <- Mj_next
      success_now <- success_next    
      
      # store results
      simu_weed_iter[c(1:np_all)+(pp-1)*np_all*3,4+kk] <- success_now
      simu_weed_iter[c((np_all+1):(2*np_all))+(pp-1)*np_all*3,4+kk] <- spent_now
      simu_weed_iter[c((2*np_all+1):(3*np_all))+(pp-1)*np_all*3,4+kk] <- pf_now[[1]]   
      
    } # end of kk iteration (reallocation events)
  } # end of pp strategy iteration
  
  # return
  simu_weed_iter
  
} # end of scneario iterations                                

# close parallel computing cluster
stopCluster(cl)



# 5) Summarize management performance under the three strategies --------------
# - 5.1) Total weed risk values remains at terminal time
# - 5.2) Relative performance (S2/S3 compared to S1; S3 compared to S2)
# - 5.3) Number of species successfully eradicated (for initial vs. late emerging species)

simu_weed_out_long <-  simu_weed_out
simu_weed_out_long$strategy <- rep(c("S1", "S2", "S3"), each = np_all*3)

# - 5.1) Total weed risk values remains at terminal time
risk_reduced_T <- simu_weed_out_long %>%
  rename(Sp_No=spID)%>%
  filter(metrics == "success") %>%
  left_join(dat_weed[,c("Sp_No", "Risk")]) %>%
  mutate(rrT = Risk * ifelse(t20==1,1,0)) %>%
  group_by(strategy, iter) %>%
  summarise(sumrrT = sum(rrT)) %>%
  ungroup() %>%
  group_by(strategy) %>%
  summarise(mean = mean(sumrrT), sd = sd(sumrrT))


# # plot
# datPlot1 <- simu_weed_out_long %>%
#   rename(Sp_No=spID)%>%
#   filter(metrics == "success") %>%
#   left_join(dat_weed[,c("Sp_No", "Risk")]) %>%
#   mutate(rrT = Risk * ifelse(t20==1,1,0)) %>%
#   group_by(strategy, iter) %>%
#   summarise(sumrrT = sum(rrT)/(10^6)) %>%   # total risk value reduction under each strategies, simmarized over iterations
#   ungroup() 
# 
# p <- ggplot(aes(x= strategy, y= sumrrT), data=datPlot1)
# p + geom_violin() + 
#   geom_boxplot(width = 0.1) +
#   theme_bw() + 
#   labs(x="Reallocation strategy", y = "Total weed risk reduction (million AWRA unit)")




# - 5.2) Relative performance (S2/S3 compared to S1; S3 compared to S2)
r_risk_reduced_T <- simu_weed_out_long %>%
  rename(Sp_No=spID)%>%
  filter(metrics == "success") %>%
  left_join(dat_weed[,c("Sp_No", "Risk")]) %>%
  mutate(rrT = Risk * ifelse(t20==1,1,0)) %>%
  group_by(strategy, iter) %>%
  summarise(sumrrT = sum(rrT)) %>%
  spread(key = strategy, value = sumrrT) %>%
  mutate(rS2 = 100*(S2-S1)/S1, rS3=100*(S3-S1)/S1,
         rS3_S2 = 100*(S3-S2)/S2) %>%
  ungroup() %>%
  summarise(mean_rS2 = mean(rS2), sd_rS2 = sd(rS2),    # relative performance (%) of S2 compared to S1, simmarized over iterations
            mean_rS3 = mean(rS3), sd_rS3 = sd(rS3),    # relative performance (%) of S3 compared to S1, simmarized over iterations
            mean_rS2_S3 = mean(rS3_S2), sd_rS3_S2 = sd(rS3_S2))  # relative performance (%) of S3 compared to S2, simmarized over iterations

# # plot
# datPlot2 <- simu_weed_out_long %>%
#   rename(Sp_No=spID)%>%
#   filter(metrics == "success") %>%
#   left_join(dat_weed[,c("Sp_No", "Risk")]) %>%
#   mutate(rrT = Risk * ifelse(t20==1,1,0)) %>%
#   group_by(strategy, iter) %>%
#   summarise(sumrrT = sum(rrT)) %>%
#   spread(key = strategy, value = sumrrT) %>%
#   mutate(rS2 = 100*(S2-S1)/S1, rS3=100*(S3-S1)/S1,
#          rS2_S3 = 100*(S3-S2)/S2) %>%
#   ungroup() %>%
#   gather(key = strategy, value = rrT, rS2, rS3)
# 
# p <- ggplot(aes(x= strategy, y= rrT), data=datPlot2)
# p + geom_violin() + 
#   geom_boxplot(width = 0.1) +
#   theme_bw() + 
#   labs(x="Reallocation strategy", y = "Relative performance in total weed risk reduction\ncompared to static S1 (%)")


# - 5.3) Number of species successfully eradicated (initial vs. late emerging species)
n_success_T <-simu_weed_out_long %>%
  filter(metrics == "success" & t20 == 1) %>%
  mutate(group = ifelse(spID %in% p_t0, "initial sp.", "emerging sp.")) %>%
  group_by(strategy, group, iter) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(strategy, group) %>%
  summarise(mean = mean(count), sd = sd(count))   # number of initial/late-emerging species eradicated under each strategies, 
# simmarized over iterations





# save entire workspace
filename <- paste0("Ch2_main_weed_", n_iter,
                   "_f", as.numeric(command_args[1]), 
                   "np", as.numeric(command_args[2]),
                   ".RData")
save.image(filename)
cat("Workspace successfully saved (after forward simulation.")

