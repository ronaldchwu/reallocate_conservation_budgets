##########################################################################################
#
# SDP-MonteCarlo model to compare PPP, interative reallocation, and SDP optimal policy
#
#  - Main output: relative performance of SDP, IR, PPP over 5000 simulations
#                 (for one parameter set; need to examine across many parameter set (e.g. n=1000))

# ----------------------------------------------------------------------------------------
# BIRD CONTEXT
#   - 2-4 decision time steps
#   - 5 sp. projects (all initially present), each with state var
#       - 0: extinct
#       - 1: Critically endangered
#       - 2: Endangered
#       - 3: Threatened#
#   - Budget var: 0,1,2,....,240
#   - Action: 0,2,4,6,8
#   - Total budget: 20,30,40,50,60,70,80, 90, 100, 110, 120, 130, 140, 150, 160
#       - assessment cost for each sp. project = 1 per time step = 5 per time step for all species
#   - Utilities (end time)
#       - -0.005 if Threatneed; -0.05 if endagered; -0.5 if critically endangered; -1 if extinct
#   - Transition probability: assume change in state level at most in a time step
#       - Getting worse p10/p21/p32: decreasing function of investment level with decreasing marginal return (exponential decrease from cap)
#       - Improving p12/p23: increasing function of investment level (linear from 0.05 to cap)
#   - Three level of investment effectiveness
#       - low: prob cap at 0.1
#       - medium: prob cap at 0.4
#       - high: prob cap at 0.7
#   - Action contraints
#       - if budget cannot afford the expenses, transition into a special DEADEND state with -9999 utility


# Set library path to folder with MDPtoolbox (and linprog)
.libPaths("/home/chunghueyw/R/x86_64-pc-linux-gnu-library/3.5")

# load commnad-line argument (in Spartan .slurm submission)
#  - Arg1: random seed number (1-500)
#  - Arg2: prob of new project emergence (0.25, 0.5, 0.75)
#  - Arg3: number of initially present projects (0-5)

command_args <- commandArgs(trailingOnly = TRUE)
seed_num <- as.numeric(command_args[1]) 
probcap <- as.numeric(command_args[2])   # 0.1: low, 0.4:medium, 0.7: high
rcost_unit <- as.numeric(command_args[3]) # 0 or 1 


#DEBUG
# seed_num <- 1
# probcap <- 0.1
# rcost_unit <- 1


# load packages
library(tictoc)
library(doSNOW)
library(parallel)
library(plyr)
library(tidyverse)
library(MDPtoolbox)



# max number of decision time steps
nT <- 4

# number of species projects
N <- 5
ninitP <- 5  # number of projects initially present

# Define utility/penalty score
penalty_scores <- c(-1,-0.5,-0.05,-0.005)
names(penalty_scores) <- c(0,1,2,3)

# marginal effectiveness of investment level on state transition
#  - from a better to a poorer state: prob = exp(-rs * investment level) (exponential decay)
#  - from a poorer to a better state: prob = rs * investment level / (maximum investment level**up_slope_adj) (linear)
set.seed(seed_num)
rs <- runif(N, 0.75, 1.25)
up_slope_adj <- 1.5



# Define actions
A_one_sp <- c(0,2,4,6)                      # possible investment level for each species project
nA_one_sp <- length(A_one_sp)          
A_all <- expand.grid(rep(list(A_one_sp),N)) # possible investment portfolio for all projects
nA <- nrow(A_all)
A_all$A_id <- c(0:(nA-1))                   # investment portfolio id; start from 0

# Set budget scenario (initial budget, per time step income from t=0 on)
set_B <- c(20,30,40,50,60,70,80,
           90,100,110,120,130,140,150,160)


# Define system state vars
S_b <- seq(0,max(set_B))                                 # possible budget variable
S_one_sp <- c(0:3)                                       # possible state of a species project
nS_one_sp <- length(S_one_sp)
S_all_sp <- expand.grid(rep(list(S_one_sp),N))           # possible state of project portfolio
nS_all_sp <- nrow(S_all_sp)
S_all_sp$all_sp_id <- c(0:(nS_all_sp-1))                 # state portfolio if; start from 0
S_all <- expand.grid(S_all_sp$all_sp_id, S_b)
S_all <- rbind(S_all, data.frame(Var1=-999, Var2=-999))  # possible state of state portfolio + budget
nS <- nrow(S_all)
S_all$all_id <- c(0:(nS-1))                              # full state id; start from 0
deadendid <- max(S_all$all_id)                           # an additional state for directing invalid action into;
                                                         # a utility of -9999 will be associated with this state




#==================================================================================================
# COMPUTING transition matrix and reward matrix of the SDP problem --------------------------------


cat("Start computing P matrix....\n")
# list to store transition and reward matrix (P and R) across all actions
PR_all <- vector('list',nA)

# Iterate through all possible actions
PR_all <- foreach(a = 1:nA, .packages = c("tidyverse")) %do% {
  
  # debug
  #a <- 1009
  
  tic()
  cat(paste0("Start computing P and R under action ", a,".\n"))
  
  
  # Get action id and corresponding individual investment levels for eahc species project
  a_id <- A_all[a,"A_id"]
  actions <- unlist(A_all[a,c(1:N)])
  a_cost <- sum(actions)
  
  
  #-------------------------------------------------------------------------------------------------------------------  
  # Step 1: Compute transition matrix preP of species projects only (ignore budget)
  #         store results into the correct budget level block
  #    - Output: a matrix of size nS_all_sp x nS_all_sp
  
  # Sparse matrix to store the transition probability
  #   - transition from row state to col state
  preP <- Matrix::Matrix(0, nrow = nS_all_sp, ncol = nS_all_sp, sparse = TRUE)
  
  
  for(ss in 1:nrow(preP)){
  
    # DEBUG
    #ss <- 38
    
    # Get starting project state (row state)
    s1 <- S_all_sp[ss,c(1:N)]
    
    
    # Step 1-1. For each project, compute the prob of transitioning into a new state
    #    - output: a matrix (nS_one_sp x N)
    
    # matrix to store transition prob of individual project
    p_indiv_sp <- array(0, c(nS_one_sp,N))
    
    for(i in 1:N){
      if(s1[i]==0){         # if sp is extinct
        p_indiv_sp[0+1,i] <- 1
      }else if(s1[i]==1){   # if sp project is critically endangered
        p_indiv_sp[0+1,i] <- probcap*exp(-rs[i]*actions[i])                  # prob of decaying to poorer state is an exponential function of effort 
        p_indiv_sp[2+1,i] <- min(1-p_indiv_sp[0+1,i],probcap*rs[i]*actions[i]/(max(A_one_sp)*up_slope_adj))   # prob of improving to next state is a linear function of effort
        p_indiv_sp[1+1,i] <- 1-p_indiv_sp[2+1,i]-p_indiv_sp[0+1,i]
      }else if(s1[i]==2){   # if sp project is endagered
        p_indiv_sp[1+1,i] <- probcap*exp(-rs[i]*actions[i])                  # prob of decaying to poor state is an exponential function of effort 
        p_indiv_sp[3+1,i] <- min(1-p_indiv_sp[1+1,i],probcap*rs[i]*actions[i]/(max(A_one_sp)*up_slope_adj)) 
        p_indiv_sp[2+1,i] <- 1-p_indiv_sp[1+1,i]-p_indiv_sp[3+1,i]
      }else if(s1[i]==3){   # if sp project is threatended
        p_indiv_sp[2+1,i] <- probcap*exp(-rs[i]*actions[i])                  # prob of decaying to poor state is an exponential function of effort 
        p_indiv_sp[3+1,i] <- 1 - p_indiv_sp[2+1,i]
      }
    } # end of i loop
    
    # Step 1-2. Compute the joint transition prob of all projects, store into preP
    p_all_sp <- expand.grid(p_indiv_sp[,1], p_indiv_sp[,2], p_indiv_sp[,3],
                            p_indiv_sp[,4], p_indiv_sp[,5]) %>%
      mutate(prob = Var1*Var2*Var3*Var4*Var5, 
             ID = c(1:nS_all_sp)-1)
    
    preP[ss,] <- p_all_sp$prob/sum(p_all_sp$prob)  # adjust for rounding errors to ensure sum to 1
    
  } # end of ss loop
  
  
  #-------------------------------------------------------------------------------------------------------------------  
  # Step 2: Copy preP into the correct budget level block in the full transition matrix P
  #

  
  # Find row/column location of each budget blocks S_all
  #   - location 1,2,3.... = budget 0,1,2...
  block_start <- which(S_all$Var1==0)
  block_end <- which(S_all$Var1==nS_all_sp-1)
  

  # Identify budget level blocks at which the action portfolio is invalid (not affordable),
  # for these blocks, set prob = 1 for transitioning into deadend state
  if(a_cost==0){  # if no investment
    P_deadend <- Matrix::Matrix(0, nrow = 0, ncol = 0, sparse = TRUE)
  }else{          # in total investment level > 0
    deadend_rows <- c(block_start[1]:block_end[a_cost])    # all blocks with budget < a_cost will transition into deadend
    
    # create the sparse deadend matrix (to be later combined into a full transition matrix P)
    # - prob = 1
    P_deadend <- Matrix::Matrix(0, nrow = length(deadend_rows), ncol = nS, sparse = TRUE)
    P_deadend[, deadendid+1] <- 1
  }

  
  # For budget level blocks where actions are valid, copy paste preP into the correct column blocks
  
  # list to store transition matrix of each budget-level row blocks
  list_P <- list()
  
  # for the first row block, budget after expenditure will be zero; so paste preP matrix onto col block of b=0
  list_P[[a_cost+1]] <- cbind(preP,
                              Matrix::Matrix(0, nrow = nS_all_sp, ncol = nS-nS_all_sp, sparse = TRUE))
  
  # for subsequent row blocks, paste preP matrix onto col block of b-a_cost-1
  for(rb in (a_cost+2):(max(S_b))){
    
    list_P[[rb]] <- cbind(Matrix::Matrix(0, nrow = nS_all_sp, ncol = block_end[rb-1-a_cost], sparse = TRUE),
                         preP,
                         Matrix::Matrix(0, nrow = nS_all_sp, ncol = nS-nS_all_sp-block_end[rb-1-a_cost], sparse = TRUE))
  }
  
  # add the last row-block: budget is at maximum budget (max(S_b))
  rb <- max(S_b)+1   
  if(a_cost==0){
    
    # if a_cost == 0, budget remain the same; attach a single column of deadend state with prob = 0
    list_P[[rb]] <-cbind(Matrix::Matrix(0, nrow = nS_all_sp, ncol = block_end[rb-1], sparse = TRUE),
                        preP,
                        Matrix::Matrix(0, nrow = nS_all_sp, ncol = 1, sparse = TRUE))
  }else if(a_cost==max(S_b)){
    # if a_cost == max budget, budget after investment will go to zero
    list_P[[rb]] <-cbind(preP,
                        Matrix::Matrix(0, nrow = nS_all_sp, ncol = nS-nS_all_sp, sparse = TRUE))
  }else{
    list_P[[rb]] <-cbind(Matrix::Matrix(0, nrow = nS_all_sp, ncol = block_end[rb-1-a_cost], sparse = TRUE),
                        preP,
                        Matrix::Matrix(0, nrow = nS_all_sp, ncol = nS-nS_all_sp-block_end[rb-1-a_cost], sparse = TRUE))
  }
  

  # last row of P matrix: the deadend
  #  - once in deadend, prob=1 to stay in deadend
  P_last_r <- Matrix::Matrix(0, nrow = 1, ncol = nS, sparse = TRUE)
  P_last_r[1,nS] <- 1
  
  # combine into full transition matrix P
  if(a_cost==0){
    P <- rbind(do.call(rbind,list_P), P_last_r)
  }else{
    P <- rbind(P_deadend, do.call(rbind,list_P), P_last_r)
  }

  rm(list_P)
  
  #------------------------------------------------------------------------------------------------
  #  Step3: Compute reward matrix R
  #      - all zeros, as we only care about end-time utility, which will be specfied in terminal reward matrix (RT)
  
  R <- Matrix::Matrix(0, nrow = nS, ncol = nS, sparse = TRUE)
  
  
  # Return a list containing the P and R matrices under the given action
  PR <- list(P, R)
  cat(paste0("P and R under action ", a, " successfully generated.\n"))
  PR
  
  #debug
  #PR_all[[a]] <- PR
  
}  # end of action loop



# Organize results into full P, R matrices (over different actions)
P <- lapply(PR_all, '[[',1)
R <- lapply(PR_all, '[[',2)

cat("Full P and R successfully generated.\n")



# Define terminal reward matrix
#  - if a species project is in poor or improving state at terminal time, 
#    it incurs negative utility (penalty) equals to its risk score
#  - terminal reward is the sum of all such penalities
RT <- mapply(function(x){
  sum(penalty_scores[as.character(unlist(S_all_sp[x+1,1:N]))])
  }, x = S_all$Var1,
  SIMPLIFY = TRUE)
RT[nS] <- -9999     # assign -9999 to deadend state




#==============================================================
# Run stochastic dynamic programming
#  - for different time steps

results_2ts <- mdp_finite_horizon(P, R, discount = 0.95, N = 2, h = RT)  # note: discount rate does not matter here
                                                                         #       as we care about terminal reward only
 
opt_policy_2ts <- results_2ts$policy


results_3ts <- mdp_finite_horizon(P, R, discount = 0.95, N = 3, h = RT)  # note: discount rate does not matter here
#       as we care about terminal reward only

opt_policy_3ts <- results_3ts$policy


results_4ts <- mdp_finite_horizon(P, R, discount = 0.95, N = 4, h = RT)  # note: discount rate does not matter here
#       as we care about terminal reward only

opt_policy_4ts <- results_4ts$policy

optimal_policy_list <- list(opt_policy_2ts, opt_policy_3ts, opt_policy_4ts)
names(optimal_policy_list) <- c("2ts","3ts","4ts")


#==============================================================
# Forward simulation
#  - simulate mgmt performance under different action strategoes
#      - no action
#      - random
#      - one-off PPP: allocate budget at t=0 based on return-on-investment of existing project; do not allow further modification
#      - iterative reallocation: redo PPP at each decision time steps
#      - SDP optimal policy: follow optimal state-and time-dependent actions prescribed by SDP


# load functions
source("./gen_PPP.R")
source("./gen_actions.R")
source("./gen_stateprob.R")
source("./gen_randomA.R")
source("./gen_next_states.R")




# set per-time-step reassessment costs for all projects
rcost_1ts <- rcost_unit*N


# Decide # of iterations to simulate
n_iter <- 5000

# Defin mgmt policies
set_policy <- c("DN","random","PPP","IRc","IRuc","SDP")
n_policy <- length(set_policy)

# List to store results under each budget level
simu_out_tblist <- list()

for(tts in 2:4){
for(bb in 1:length(set_B)){
  
  # DEBUG
  # tts <- 4
  # bb <- length(set_B)
  # # 
  # Decide budget level
  b_init <- set_B[bb]                           
  
  # select the optimal policy for the number of time steps to simulate
  opt_policy <- optimal_policy_list[paste0(tts,"ts")][[1]]
  
  # doParallel setup
  cl <- makeCluster(8)
  registerDoSNOW(cl)
  
  
  # Generate transition matrix of site states (all sites only) under each action-----------
  # progress tracking
  pb <- txtProgressBar(max = n_iter, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  cat("Start forward simulation....")
  
  # ---- loop over many iterations (w/ same habitat quality change realization)
  simu_out <- foreach(iter = 1:n_iter, .combine = "rbind", .options.snow = opts, .packages = c("tidyverse")) %dopar% {
    
    
    # DEBUG
    #iter <- 1
    
    #-----------------------------------------
    # Step 1: Initiation
    
    
    
    # initial status of species projects and budget
    s_init <- c(1,1,2,2,3,b_init)                 # initial system full state
    allspID_init <- unlist(subset(S_all_sp,       # get ID of the project portfolio
                                  Var1 == s_init[1] &
                                    Var2 == s_init[2] &
                                    Var3 == s_init[3] &
                                    Var4 == s_init[4] &
                                    Var5 == s_init[5],
                                  select = "all_sp_id"))
    ID_init <- unname(unlist(subset(S_all,               # get ID of the full system state
                                    Var1 == allspID_init &
                                      Var2 == b_init,
                                    select = "all_id")))
    s_init <- c(s_init, ID_init)           
    s_now <- rep(list(s_init), n_policy)          # initiate full system state for each mgmt policy trajectory (as a list)
    
    
    # for IRc, IRuc, SDP policies, start with budget = TB - reassessment costs over the time horizon
    s_init_rc <- c(1,1,2,2,3,b_init-rcost_1ts*tts)                 # initial system full state
    allspID_init_rc <- unlist(subset(S_all_sp,       # get ID of the project portfolio
                                  Var1 == s_init_rc[1] &
                                    Var2 == s_init_rc[2] &
                                    Var3 == s_init_rc[3] &
                                    Var4 == s_init_rc[4] &
                                    Var5 == s_init_rc[5],
                                  select = "all_sp_id"))
    ID_init_rc <- unname(unlist(subset(S_all,               # get ID of the full system state
                                    Var1 == allspID_init_rc &
                                      Var2 == b_init-rcost_1ts*tts,
                                    select = "all_id")))
    s_init_rc <- c(s_init_rc, ID_init_rc)           
    
    s_now[[4]] <- s_init_rc
    s_now[[5]] <- s_init_rc
    s_now[[6]] <- s_init_rc
    
    
    # create data.frame to store state trajectories in this iteration
    simu_out_temp <- data.frame(iter=iter,
                                t=rep(c(0:(tts)), each = n_policy),
                                policy=rep(set_policy, tts+1),
                                sp1=s_init[1],
                                sp2=s_init[2],
                                sp3=s_init[3],
                                sp4=s_init[4],
                                sp5=s_init[5],
                                budget=s_init[6],
                                fullID=s_init[7],
                                a_sp1 = 0,         # actions taken (note that no actions are taken at t = tts + 1)
                                a_sp2 = 0,
                                a_sp3 = 0,
                                a_sp4 = 0,
                                a_sp5 = 0)
    
    # store initial system states
    for(p in 1:n_policy){  
      # debug
      #p <- 1
      rowloc <- p
      simu_out_temp[rowloc,c(2,4:10)] <- c(0,s_now[[p]])
    }
    
    
    #-----------------------------------------
    # Step 2: Loop over the tt time steps, simulate state transitions -----
    
    # action tracker
    a_last <- rep(list(rep(0,N)),n_policy)
    
    for(t in 0:(tts-1)){
      
      # debug
      # t<-0
      # t<-t+1
      # t
      # s_now
      # s_now_debug <- s_now
      # a_ids
      # #
      # s_now<-s_now_debug
      
      # Step 2-1. Select actions based on each mgmt policy
      a_ids <- gen_actions(s_now, a_last, t, opt_policy)   # list of actionids corresponding to each policy
      
      # Step 2-2. Simulate state transition (stochastic)
      s_next <- gen_next_states(s_now, a_ids)
      
      for(p in 1:n_policy){  # store system states after transition
        # debug
        #p <- 1
        rowloc <- p + (t+1)*(n_policy)
        simu_out_temp[rowloc,c(2,4:10)] <- c(t+1,s_next[[p]])
        simu_out_temp[rowloc-n_policy,c(11:15)] <- a_ids[[p]]  # store the actions allocated at each site
      }
      
      
      # update system states
      s_now_debug <- s_now
      s_now <- s_next
      a_last <- a_ids
      
    } # ---- end of time loop -----
    
    
    # return results
    simu_out_temp
    
    # debug
    #simu_out[[iter]] <- simu_out_temp
    
  } # ---- end of iteration loop; results stored in simu_out -----
  close(pb)
  stopCluster(cl)
  
  simu_out$TB <- b_init
  simu_out$ts <- tts
  simu_out_tblist[[bb+(tts-2)*length(set_B)]] <- simu_out
  
} # end of bb loop
} # end of tts loop

# ----- END of simulation -------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------------------------
# Compute relative performance (difference)
# 
# # combine simulation results at each budget levels
simu_out_all <- do.call(rbind, simu_out_tblist)

# compute end-time utility under each policy trajectory X budget levels
#   - row: budget level, col: DN, random, PPP, IR, SDP, cell: end-time utility
#   - replicate row blocks for each iteration
endU_temp <- simu_out_all %>%
  group_by(ts) %>%
  filter(t==max(t)) %>%
  ungroup() %>%
  mutate(EU = penalty_scores[as.character(sp1)] +
           penalty_scores[as.character(sp2)] +
           penalty_scores[as.character(sp3)] +
           penalty_scores[as.character(sp4)] +
           penalty_scores[as.character(sp5)]) %>%
  select(iter,TB,ts,policy,EU) %>%
  spread(key = policy, value = EU)
#
#
# # turn into long-form:
# #   - cols: iteration, budget-level, policy, relative performance
# endU_temp <- endU_temp %>%
#   select(iter,TB,DN,random,PPP,IR,SDP) %>%
#   gather(key = "policy", -iter, -TB, value = "value")
#
# # summarize over iterations
# endU <- endU_temp %>%
#   group_by(TB, policy) %>%
#   summarise(mean = mean(value),
#             sd = sd(value),
#             min = min(value),
#             Q1 = quantile(value)[2],
#             med = quantile(value)[3],
#             Q3 = quantile(value)[4],
#             max = max(value))
#
# # plot results (800x400)
# datPlot1<-endU_temp
# datPlot1$TB <- as.factor(datPlot1$TB)
# datPlot1$policy <- factor(datPlot1$policy, levels = c("DN","random","PPP","IR","SDP"))
# p <- ggplot(data = datPlot1, aes(x=TB,y=value, fill=policy))
# p + geom_boxplot(alpha = 0.7)+
#   theme_bw() +
#   theme(axis.title = element_text(size=15),
#         axis.text = element_text(size=12))+
#   scale_fill_viridis_d()+
#   labs(x="Total budget", y="End-time utility", title="Weed context (5 projects)")


# datPlot2<-endU
# datPlot2$TB <- as.factor(datPlot2$TB)
# datPlot2$policy <- factor(datPlot2$policy, levels = c("DN","random","PPP","IR","SDP"))
# p <- ggplot(data = datPlot2, aes(x=TB, y=mean, fill=policy))
# p + geom_bar(position = "dodge", stat="identity",alpha = 0.7)+
#   theme_bw()+
#   scale_fill_viridis_d()




#=---------------------------------------------------
# compute performance difference relative to SDP
#   - (policy-DN)/(SDP-DN) * 100%
#   - add cols: drandom, dPPP, dIR
endU_temp <- endU_temp %>%
  mutate(dDN=DN-SDP,
         drandom = random-SDP,
         dPPP = PPP-SDP,
         dIRc = IRc-SDP,
         dIRuc = IRuc-SDP)


# turn into long-form:
#   - cols: iteration, budget-level, policy, relative performance
endU_temp <- endU_temp %>%
  select(iter,TB,ts,dDN,drandom,dPPP,dIRc,dIRuc) %>%
  gather(key = "policy", -iter, -TB, -ts, value = "value")

# # plot results (1000x600)
datPlot2<-endU_temp
datPlot2$TB <- as.factor(datPlot2$TB)
datPlot2$ts <- as.factor(datPlot2$ts)
datPlot2$policy <- factor(datPlot2$policy, levels = c("dDN","drandom","dPPP","dIRc","dIRuc"))
p <- ggplot(data = datPlot2, aes(x=TB,y=value, fill=policy))
p + geom_boxplot(alpha = 0.7)+
  theme_bw() +
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=12))+
  scale_fill_viridis_d()+
  facet_grid(ts~.)+
  labs(x="Total budget", y="Difference in end-time utility compared to SDP",
       title=paste0("Bird context (5 projects), probcap=",probcap," , rcost_1ts=", rcost_unit*N))


# reduce workspace file size
rm(P)
rm(R)
rm(PR_all)
rm(RT)
rm(opts)
rm(P_deadend)
rm(preP)



# save entire workspace
filename <- paste0("CH2_SDP_MC_bird_", n_iter,
                   "_sn", as.numeric(command_args[1]), 
                   "pc", as.numeric(command_args[2]), 
                   "rcu", as.numeric(command_args[3]),
                   ".RData")
save.image(filename)
cat("Workspace successfully saved (after forward simulation.")



#---------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------
#
# Result synthesize
#


endU_temp_list <- list()
decision_similarity_list <- list()
count <- 0

for(sn in 1:5){
  for(pc in c(0.1,0.4)){
    for(rcu in c(0,1)){
      
      load(paste0("./CH2_SDP_MC_bird_5000_sn",sn,"pc",pc,"rcu",rcu,".RData"))
      count <- count + 1
      
      simu_out_all <- do.call(rbind, simu_out_tblist)
      simu_out_all$sn <- sn
      simu_out_all$pc <- pc
      simu_out_all$rcu <- rcu
      
      endU_temp <- simu_out_all %>%
        group_by(ts) %>%
        filter(t==max(t)) %>%
        ungroup() %>%
        mutate(EU = penalty_scores[as.character(sp1)] +
                 penalty_scores[as.character(sp2)] +
                 penalty_scores[as.character(sp3)] +
                 penalty_scores[as.character(sp4)] +
                 penalty_scores[as.character(sp5)]) %>%
        select(sn,iter,TB,ts,policy,EU) %>%
        spread(key = policy, value = EU)
      
      endU_temp_list[[count]] <- endU_temp
      
      # Decision similarity between iterative reallocation SDP
      decision_similarity_list[[count]]  <- simu_out_all %>%
        select(iter,t,policy,a_sp1,a_sp2,a_sp3,a_sp4,a_sp5,TB,ts,sn,pc,rcu) %>%
        filter(policy %in% c("IRuc","SDP")) %>%
        group_by(iter,t,ts, TB) %>%
        mutate(diff_a_sp1 = ifelse(a_sp1 == lag(a_sp1), 1, 0),
               diff_a_sp2 = ifelse(a_sp2 == lag(a_sp2), 1, 0),
               diff_a_sp3 = ifelse(a_sp3 == lag(a_sp3), 1, 0),
               diff_a_sp4 = ifelse(a_sp4 == lag(a_sp4), 1, 0),
               diff_a_sp5 = ifelse(a_sp5 == lag(a_sp5), 1, 0)) %>%
        mutate(diff_a_all = (diff_a_sp1 + 
                               diff_a_sp2 +
                               diff_a_sp3 +
                               diff_a_sp4 +
                               diff_a_sp5)/5) %>%
        ungroup() %>%
        group_by(t,ts,TB) %>%
        summarise(perc_similarity_mean = mean(diff_a_all, na.rm = TRUE),
                  perc_similarity_sd = sd(diff_a_all, na.rm = TRUE)) %>%
        mutate(sn = sn,
               pc = pc,
               rcu = rcu)
      
    }
  }
}



#-------------------------------------------------------------------------------------------------------------
# Compute relative performance (difference)

# combine simulation results at each budget levels
endU_temp <- do.call(rbind, endU_temp_list)
deci_simi_temp <- do.call(rbind, decision_similarity_list)



#=---------------------------------------------------
# compute performance difference relative to SDP
#   - (policy-DN)/(SDP-DN) * 100%
#   - add cols: drandom, dPPP, dIR
endU_temp <- endU_temp %>%
  mutate(dDN=DN-SDP,
         drandom = random-SDP,
         dPPP = PPP-SDP,
         dIRc = IRc-SDP,
         dIRuc = IRuc-SDP)


# turn into long-form:
#   - cols: iteration, budget-level, policy, relative performance
endU_temp <- endU_temp %>%
  select(sn,iter,TB,ts,dDN,drandom,dPPP,dIRc,dIRuc) %>%
  gather(key = "policy", -sn, -iter, -TB, -ts, value = "value")

# plot results (1280x800)
datPlot2<-endU_temp
datPlot2 <- subset(datPlot2, TB<=80)
datPlot2$TB <- as.factor(datPlot2$TB)
datPlot2$ts <- as.factor(datPlot2$ts)
datPlot2$policy <- factor(datPlot2$policy, levels = c("dDN","drandom","dPPP","dIRc","dIRuc"))
datPlot2 <- subset(datPlot2, sn %in% c(1:3))
p <- ggplot(data = datPlot2, aes(x=TB,y=value, fill=policy))
p + geom_boxplot(alpha = 0.7)+
  geom_abline(slope = 0, intercept = 0) +
  theme_bw() +
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=12))+
  scale_fill_viridis_d()+
  facet_grid(ts~sn)+
  labs(x="Total budget", y="Difference in end-time utility compared to SDP",
       title=paste0("Bird context, probcap=",probcap," , rcost_1ts=", rcost_unit*N))






# Plot showing decision similarities across budget level, # decisio time steps
#  - X: budget level
#  - Y: % decision similarity (at different time steps)
#  - facet: Y- # of decision time periods


# 800 x 600
datPlot3 <- deci_simi_temp %>% 
  ungroup() %>%
  group_by(t,ts,TB,pc,rcu) %>%
  filter(t!=ts) %>%  # no actions were taken at the end time step. Remove
  filter(pc==0.1 & rcu==1) %>%
  summarise(meanDS = 100*mean(perc_similarity_mean))
datPlot3$TB <- as.factor(datPlot3$TB)
datPlot3$t <- as.factor(datPlot3$t)


p <- ggplot(data = datPlot3, aes(x=TB, y=meanDS,
                                 group=t, linetype=t, shape=t, color=t))
p + geom_line() +
  geom_point() +
  theme_bw() +
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=12))+
  scale_colour_colorblind() +
  facet_grid(ts~.)+
  labs(x="Total budget", y="Decision similarity between\n iterative reallocation and SDP policies (%)",
       title=paste0("Bird context, probcap=",probcap," , rcost_1ts=", rcost_unit*N))




