########################################################################
#
# Budget reallocation analysis: threatended birds in Victoria, Australia
#

# Chung-Huey Wu (chung.huey.wu@gmail.com)


# Set library path to folder containing MDPtoolbox (and linprog)
.libPaths("/home/chunghueyw/CH2_reallocation/lib")


# Load packages
library(tidyverse)
library(tictoc)
library(doSNOW)
library(tidyverse)
library(gurobi)

# Load functions
source("./gen_ROI.R")
source("./portfolio_opt.R")
source("./IUCNtransition.R")



#------------------------------------------------------------------------------
# 1) Import data of bird projects 

# data of bird species
dat_bird <- read.csv("./data/AU_bird.csv")

# assign which bird species (species ID) are initially present at t = 0
#  - assume all species are initially present and no new emerging species
p_t0 <- dat_bird$SpID
np_t0 <- length(p_t0)        # number of initial species
np_all <- np_t0              # total number of species



#------------------------------------------------------------------------------
# 2) Setup decision contexts, reassessment scheme, and other parameters 

# Time horizon and reassessment-reallocation schedule
# - note: time start from t = 0
T_unit <- 8                                # time step unit (year)
T_h <- 16/T_unit                           # time horizon
n_T <- T_h/T_unit                          # number of time steps over the time horizon
T_Aj <- c(0,8, 16)/T_unit                  # which time steps to assess bird IUCN status
T_Rk <- c(0,8)/T_unit                      # which time steps to make budget reallocation decisions
J <- length(T_Aj)                          # number of reallocation events
K <- length(T_Rk)                          # number of reallocation events

# Utility values for a species being in each IUCN status
#  - in terms of penality to be minimized
setIUCNstatus <- c("EX", "CR", "EN", "VU", "NTLC")
n_status <- length(setIUCNstatus)
penalty_weights <- c(1.0, 0.5, 0.05, 0.005, 0)  # pentalty score for EX, CR, EN, VU, NTLC


# Management budget, and costs of reassessment & reallocation
#   - throughout the analyses, all monetary values are discounted and converted to value at t = 0
B_unit <- 100000         # budget unit (AU$ 0.1 million)
B_all <- 576/4           # budget level used in our analysis; similar to 
                         # real world budget (see McCarthy et al., 2008)

set_Ca <- c(0, 25400)/B_unit
Ca <- set_Ca[2]          # cost of reassessing IUCN status for all species 
Cr <- 0                  # cost of conducting budget reallocation analysis and decision-making (for entire portfolio)


# Total budget available for direct management of bird projects
#   - assume assessment costs at initial and end time steps are from
#     other sources, so Ca incurs only J-2 times over the time horizon

B_avail <- B_all - Ca*np_t0*(J-2) - Cr*K      # B_all minus preallocated reassessment & reallocation budgets 


# Management strategies with varying flexibility in reallocation
set_strategy <- c("S1","S2", "S3")
n_strategy <- length(set_strategy)
iss_constraint <- TRUE    # whether we force species in the same IUCN group to be allocated the same budget level


# Project Return on Investment (ROI):  
#  - Investment: monetary budget allocated to a bird project till the end of time horizon
#  - Return: expected penality score of the species by end time
min_xi <- 0                       # minimum per-time-step budget allocated to a bird project (unit is B_unit)
max_xi <- 10                      # maximum per-time-step budget allocated to a bird project (unit is B_unit)
bin_xi <- 0.25                    # bin size of per-time-step  investment level
xi <- seq(min_xi, max_xi, by = bin_xi)   # all possible per-time-step budget
n_xi <- length(xi)                       # number of possible per-time-step budget
xi_all <- rep(list(xi), np_all)          # possible per-time-step budget for each bird species
                                         # here assumed to be all the same
neighbor_dist <- 0.0001           # here assumed to be all the same

# Forward simulation
n_iter <- 5000       # number of bird management simulations to run


#------------------------------------------------------------------------------
# 3) Prepare bird-specific system models ------------------

# - 3.1) IUCN-transition model for computing the probability of transitioning between IUCN status based on investment level (McCarthy et al., 2008)
#  - ----> see IUCNtransition.R for script details        

# Regression coefficients of how monetary investment level correlates with past observed changes in bird IUCN status over 8-yr period (see McCarthy et al., 2008)
par_mean <- c(-2.944, -1.157, -2.641, -0.3843, -3.307, 0.03525, -4.346, -1.03)
par_std <- c(0.3209, 0.8135, 0.2895, 0.2724, 0.442, 0.06975, 0.8326, 1.519)



#------------------------------------------------------------------------------
# 4) Iterative reallocation analysis under forward simulation
#    - Implement our proposed reallocation approach
#    - Multiple management strategies with varying reallocation flexibility are parallelly examined in each simulation iteration: 
#       - S1-"Static": no budget reallocation
#       - S2-"Constrained reallocation": only allow budget released from terminated projects (bird species extinction) to be invested in other bird projects
#       - S3-"Full reallocation": budget can be reallocated among bird projects with no constriants
#       - The three strategies are simulated using the same set of IUCN-state transition random seeds (see below). So their conservation outcomes can be directly compared within each simulation iteration.
#
#    - Steps
#      - 4.1. Generate ROI for all exisiting bird projects
#           - use IUCN transition model to compute expected utility (based on terminal IUCN status) at end time 
#      - 4.2. Find optimal budget allocation among budgets and switch to such new budget portflio
#           - use separable programming
#           - Following the original study, we assume species with same IUCN status receive same amount of investment
#            - for S1-"Static-no-reassessment" strategy: continue using the ongoing portfolio; 
#              remaining budget of terminated projects won't be available to other projects
#            - for S2-"Static-with-reassessment" strategy: same as S1, but spend budget on reassessing all species
#            - for S2-"Constrained reallocation": budget of all ongoing projects cannot be modified. Since there are new bird projects emerging, there will be no reallocation. So S2 is equivalent to S2.
#             - for S3-"Full reallocation": adopt the new optimal porfolio    
#      - 4.3. Simulate bird conservation progress till the next reallocation event
#           - use IUCN transition model to generate probability of status changes, and then compare to random seeds to 
#             determine realized status transitions.
#
#      - Repeat 4.1-4.3 until the end of time horizon

# Setup parallel computing
cl <- makeCluster(10)
registerDoSNOW(cl)



# progress tracking in R
pb_P <- txtProgressBar(max = n_iter, style = 3)
progress_P <- function(n) setTxtProgressBar(pb_P, n)
opts_P <- list(progress = progress_P)


# Iterations over simulated bird management
simu_bird_out <- foreach(s = 1:n_iter,
                         .combine = "rbind",
                         .options.snow = opts_P,
                         .errorhandling = "remove",
                         .packages = c("tidyverse", "gurobi"))%dopar%{
  
  # for(s in 1:n_iter){  
  # debug
  # s <- 1
  
  # Data frame to store bird projects status and management history over time in this scenario iteration
  #  - IUCN status: at each reassessment timing (e.g. t=0,8,16)
  #        - id: EX=0, CR=1, EN=2, VU=3, NTLC=4
  #  - m_spent: cumulative monetary budget spent on each project by each assessment time step (unit: B_unit)
  
  simu_bird_iter <- data.frame(iter = rep(s, np_all*2*n_strategy),
                               spID = rep(c(1:np_all), 2*n_strategy),
                               metrics = rep(c("IUCN", "m_spent"), each = np_all),
                               strategy = rep(set_strategy, each = np_all*2),
                               t8 = -1,
                               t16 = -1)
  
  
  # Generate random seeds of species' IUCN status transition in this simulation iteration
  #  - Rationale: from a current IUCN status, there are a number of possible future status a species can transition into.
  #               The higher the monetary investment into a project the better chance the species will improve its IUCN status
  #               In a simulation iteration with multiple strategies examined, we need to ensure that the strategy with higher   #               investment level MUST have equal or better per-time-step IUCN transition outcome than strategies with lower   #               investment level. 
  #  - A uniform random numbers between 0 and 1 for each species at each reassessment/reallcation time step is generated. 
  #  - The seed will be shared by system trajectories under all management strategies in this simulation iteration
  #  - Example: If the seed for Species 1 being at t=8 is 0.5, and if under strategy S1 and S2, the prob of Species 1 transitioning into {EX, CR, EN, VU, NTLC}(worst ot best) is [0.3,0.3,0.2,0.1,0.1] and [0.1,0.2,0.3,0.3,0.1], the cumulative prob of transitioning prob are [0.3,0.6,0.8,0.9,1] and [0.1,0.3,0.6,0.9,1], then the seed of 0.5 falls into 'CR' under S1 and EN under S2. As a result, Speciea 1 will transition into 'CR' and 'EN' respectively.
  
  #set.seed(1000)
  seed_success <- matrix(runif(np_all*K), nrow = np_all, ncol = K)
  
  
  
  # Parallel simulation of bird conservation progress under each mgmt strategy (S1:no reallocation no reassessment, S2: reassessment but no reallocation, S2: constrained reallocation, S3: full reallocation)
  for(pp in 1:n_strategy){
    
    # debug
    # pp <-3
    
    # select management strategy
    strategy <- set_strategy[pp]
    
    # - 4.0) Set up initial conditions of bird projects at t = 0, and initialize variables that tracks bird project progress
    
    # Variables that track management budget and current time step
    B_avail_now <- B_avail
    t_now <- 0
    
    
    # if strategy is S1: no ressessment no reallocation, add the reassessment costs back to budget pool
    if(strategy == "S1"){
      B_avail_now <- B_avail_now + Ca*np_t0*(J-2)
    }
    
    # Variables that track project status and management efforts of each of the 270 bird species
    p_now <- p_t0                        # exsiting bird species projects 
    np_now <- length(p_now)           # number of existing bird species projects 
    IUCN_now <- dat_bird$ini_StatusID    # IUCN status at t = 0
    
    
    # Variable to track the current budget allocation
    #    - 1st element of the list: allocated total monetary budget (from now to the end of time horizon)
    #    - 2nd element of the list: monetary expenses at each 8-yr period (unit time step)
    #                               (matrix, row is species, col is time step)
    #    - 3rd element of the list: expected total utility at end time under the current portfolio
    pf_now <- list(rep(0, np_all),     
                   matrix(0, nrow = np_all, ncol = n_T),     
                   0)                  
    
    # Variable that tracks cumulative monetary expenditure in each bird project (no discounting, as in McCarthy et al., 2008)
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
      
      # - 4.1. Generate ROI for all exisiting bird projects
      #     - use IUCN transition model to compute expected penality (based on terminal IUCN status) at end time 
      ROIs_now <- gen_ROI(pp, kk, prob_scale)
      
      # - 4.2. Find optimal budget allocation among budgets and switch to such new budget portflio
      #   - use separable programming
      #   - for S1-"Static-no-reassessment" strategy: continue using the ongoing portfolio; 
      #                                             remaining budget of terminated projects won't be available to other projects
      #     - for S2-"Static-with-reassessment" strategy: same as S1, but spend budget on reassessing all species
      #     - for S2-"Constrained reallocation": budget of all ongoing projects cannot be modified. Since there are new bird projects emerging, there will be no reallocation. So S2 is equivalent to S2.
      #     - for S3-"Full reallocation": adopt the new optimal porfolio    
      if(set_strategy[pp] != "S3" & kk > 1){     # for S1,S2,S2 strategy, keep using the remaining 
        # budget allocation in pf_now
        pf_next[[1]] <- rowSums(as.data.frame(pf_now[[2]][,kk:K]))  # update the allocated total monetary budget (from now to the end of time horizon)
      }else{    # proceed to separable programming optimization if it is first allocation event, or for stratgies S2 and S3 
        pf_next <- portfolio_opt(ROIs_now, B_avail_now, set_strategy[pp], pf_now, kk)  
      }
      
      # - 4.3. Simulate bird conservation progress till the next reallocation event
      #   - use IUCN transition model to generate probability of status changes, and then compare to random seeds to 
      #     determine realized status transitions.
      
      # Spend money: calculate the monetary costs of implementing the portfolio till the next reallocation event, deduce from budget pool
      
      m_spent <- as.data.frame(pf_next[[2]][,kk])
      
      # Deduce sum of these costs from budget pool 
      B_avail_next <- B_avail_now - sum(m_spent)
      spent_next <- spent_now + rowSums(m_spent)  # record the expenditure
      
      # Use IUCN transition model to generate stochastic IUCN status change by the next reallocation time step
      # set the target end time to be next realllocation event for transition model
      t_end <- kk
      
      # apply transition model to compute probability of transition into each IUCN status by the next realllocation event
      IUCN_next_prob <- mapply(function(i, x){
        IUCNtransition(i, x, t_now, t_end, prob_scale)  # prob of ending up in [EX, CR, EN, VU, NTLC] status at end time
      }, i = c(1:np_all), x = pf_next[[1]])
      
      # scale to sum to one (for rounding errors)
      IUCN_next_prob <- apply(IUCN_next_prob,2,function(col)col/sum(col))
      
      # compare the probabilities against the pre-generated random seeds to determine IUCN status transition
      cumu_prob <- apply(IUCN_next_prob,2,cumsum)   # cumulative sum from prob of entering EX to NTLC 
      IUCN_next <- mapply(function(i, seed){
        min(which(cumu_prob[,i]>seed))-1
      },i = c(1:np_all), seed = seed_success[,kk])
      
      
      # Update state variables
      B_avail_now <- B_avail_next 
      spent_now <- spent_next 
      pf_now <- pf_next    
      IUCN_now <- IUCN_next
      
      # store results
      simu_bird_iter[c(1:np_all)+(pp-1)*np_all*2,4+kk] <- IUCN_now
      simu_bird_iter[c((np_all+1):(2*np_all))+(pp-1)*np_all*2,4+kk] <- spent_now
      
    } # end of kk iteration (reallocation events)
  } # end of pp strategy iteration
  
  gc()
  
  
  # return
  simu_bird_iter
  
  # debug
  # check if at any time species of same IUCN status has same budget allocation
  # IUCN_t8 <- simu_bird_iter %>%
  #   filter(metrics == "IUCN") %>%
  #   select(-t16)
  # expense_t16 <-simu_bird_iter %>%
  #   mutate(et16 = t16-t8) %>%
  #   filter(metrics == "m_spent") %>%
  #   select(-c(t8,t16)) %>% 
  #   left_join(IUCN_t8, by = c("iter", "spID", "strategy"))
  
} # end of scneario iterations                                

# close parallel computing cluster
stopCluster(cl)



# 5) Summaize management performance under the four strategies --------------
# - 5.1) Total pentality scores at end time
# - 5.2) Number of species in each IUCN status group at end time

simu_bird_out_long <-  simu_bird_out
simu_bird_out_long$strategy <- rep(c("S1", "S2", "S3"), each = np_all*2)

# - 5.1) Total penalty scores at terminal time
end_penalty <- simu_bird_out_long %>%
  filter(metrics == "IUCN") %>%
  mutate(penalty = penalty_weights[t16+1]) %>%
  group_by(strategy, iter) %>%
  summarise(sumP = sum(penalty)) %>%
  ungroup() %>%
  group_by(strategy) %>%
  summarise(mean = mean(sumP), sd = sd(sumP))


# plot
# datPlot1 <- simu_bird_out_long %>%
#   rename(Sp_No=spID)%>%
#   filter(metrics == "success") %>%
#   left_join(dat_bird[,c("Sp_No", "Risk")]) %>%
#   mutate(rrT = Risk * ifelse(t20==1,1,0)) %>%
#   group_by(strategy, iter) %>%
#   summarise(sumrrT = sum(rrT)/(10^6)) %>%   # total risk value reduction under each strategies, simmarized over iterations
#   ungroup() 
#   
# p <- ggplot(aes(x= strategy, y= sumrrT), data=datPlot1)
# p + geom_violin() + 
#   geom_boxplot(width = 0.1) +
#   theme_bw() + 
#   labs(x="Reallocation strategy", y = "Total bird risk reduction (million AWAR unit)")




# - 5.2) Number of species in each IUCN status group at end time
end_sp_status <- simu_bird_out_long %>%
  filter(metrics == "IUCN") %>%
  group_by(strategy, iter) %>%
  summarise(nEX = length(which(t16==0)),
            nCR = length(which(t16==1)),
            nEN = length(which(t16==2)),
            nVU = length(which(t16==3)),
            nNTLC = length(which(t16==4))) %>%
  ungroup()%>%
  group_by(strategy) %>%
  summarise(meannEX = mean(nEX), sdnEX = sd(nEX),
            meannCR = mean(nCR), sdnCR = sd(nCR),
            meannEN = mean(nEN), sdnEN = sd(nEN),
            meannVU = mean(nVU), sdnVU = sd(nVU),
            meannNTLC = mean(nNTLC), sdnNTLC = sd(nNTLC))


# # # plot
# datPlot2 <- simu_bird_out_long %>%
#   filter(metrics == "IUCN") %>%
#   group_by(strategy, iter) %>%
#   summarise(nEX = length(which(t16==0)),
#             nCR = length(which(t16==1)),
#             nEN = length(which(t16==2)),
#             nVU = length(which(t16==3)),
#             nNTLC = length(which(t16==4))) %>%
#   ungroup()
# 
# p <- ggplot(aes(x= strategy, y= nEX), data=datPlot2)
# p + geom_violin() +
#   geom_boxplot(width = 0.1) +
#   theme_bw() +
#   labs(x="Reallocation strategy", y = "Number of species extinct by the year 16")
# 




# save entire workspace
filename <- paste0("Ch2_main_bird_", n_iter,
                   "_ps", as.numeric(command_args[1]), 
                   "_rcost.RData")
save.image(filename)
cat("Workspace successfully saved (after forward simulation.")

