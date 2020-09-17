##################################################################################
#
# Function to construct scenarios of new weed incursion, using Poisson process
# - inputs
#      - p_new: candidate new weed species that may emerge over time
#      - n_p_new: number of species in p_new set
#      - lambda: rate of new weed incursion; as the parameter of Poisson process
#      - n_iter: number of scnearios to generate
# - outputs
#     - scenario_pool: which weed spcecies will emerge at each time step from t=1 to end of time horizon
#     - p_incur_pool: list of weed species that will emerge in this scenario
#

IncursionScenario <- function(p_new, n_p_new, lambda, T_h, n_iter){
  
  # store n_iter number of weed incursion scenarios
  scenario_pool <- list()
  p_incur_pool <- list()
  
  for(s in 1:n_iter){
    
    if(n_p_new > 0){
      # variables to store weed incursion history over time
      incur_steps <- rep(0, T_h)  # num of weed project incurred in each time step
      incur_count <- 0    # track total num of weed project emerged over time
      scenario <- list()  # which weed project(s) emerged at each time step 
      
      # track remaining candidate species pool over time
      p_new_temp <- p_new
      n_p_new_temp <- n_p_new
      
      # simulate weed incursion over time based on Poisson process
      for(t in 1:T_h){
        
        # debug
        #t<-1
        #t<-t+1

        
        # randomly generate the num of weed incurred at this time step
        incur_steps[t] <- rpois(1, lambda)
        
        # if num of weed incursion drawn from Poisson process is larger than the num of remaining new projects, 
        # adjust the number drawn 
        if(incur_steps[t] > n_p_new_temp){
          incur_steps[t] <- max(incur_steps[t] - (incur_count - n_p_new_temp), 0)
        }
        
        # update total count of num incursion
        incur_count <- incur_count + incur_steps[t]
        
        
        # randomly assign which projects incurred in this time step
        if(n_p_new_temp>1 & incur_steps[t]>0){
          scenario[[t]] <- sample(p_new_temp, incur_steps[t])  # store which weed project incurred
        }else if(n_p_new_temp==1 & incur_steps[t]==1){
          scenario[[t]] <- p_new_temp
        }else if(incur_steps[t]==0){
          scenario[[t]] <- numeric()
        }
        
        # update remaining new weed that have not incurred
        p_new_temp <- p_new_temp[!p_new_temp %in% scenario[[t]]]
        n_p_new_temp <- length(p_new_temp)
        
        # debug
        # cat(paste0("#sp drawn: ",(incur_steps[t]),"\n"))        
        # cat(paste0("#sp incurred: ",incur_count,"\n"))
        # cat(paste0("#sp left: ",(n_p_new_temp),"\n"))        
        
      }
      
      # output variables to store the results
      scenario_iter <- scenario
      p_incur_iter <- unique(unlist(scenario))
      
    }else{    # if there is no new project at all
      scenario_iter <- list()
      p_incur_iter <- c()
    }
    
    # store resulting scenario in this iteration
    scenario_pool[[s]] <- scenario_iter
    p_incur_pool[[s]] <- p_incur_iter
  } # end of s loop
  
  return(list(scenario_pool, p_incur_pool))
}