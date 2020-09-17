###############################################################
#
# Hazard rate model: computing the probability of eradication success
#   - Based on real world data of how eradication successes correlate with weed project properties & management efforts (see Dodd et al. 2015)
#        - assumptions made: for all species, Climate Score = 3, Eradicated Elsewhere = yes, Time to Reproductive Maturity = 1yr.
#   - Modification from Dodd et al., 2017
#        - We deterministically use median value of regression coefficient under Monte Carlo sampling (provided in Supplmentary data of Dodd et al. 2015)
#        - We ignore potential species-specific random effects for the 50 species in our study.
#
#   - Inputs
#       - H: total infested area (ha)
#       - cs: climate score (1-4)
#       - DP: annual period of detectability prior to seed set (months)
#       - Wa: search width (m)
#       - dto: distance to nearest management office (km)
#       - PL: propagule longevity (years)
#       - avg_E: mean survey/treatment frequency (per year) over from t_start to t_end (see below)
#       - tm: time to reproductive maturity (years)
#       - se:  whether the species has been eradicated elsewhere (1 when yes, 0 when no)
#
#       - m_DP: 8.42 from Dodd et al. 2015
#       - m_PL: 8.12 from Dodd et al. 2015
#       - HRMpars: regression model parameters from real world data of how eradication successes correlate with weed project properties & management efforts 
#       - t_start: time step when management start (with no interuption till t_end)
#       - t_end: target time step to estimate prob. of eradication
#       - new_E: survey/treatment frequency (per year) allocated for subsequent management
#       - neighbor_dist: dummy investment level used in ROI/optimization 
#                        that in here actually means zero investment (survey freq = 0).
#
#   - Outputs
#       - pSuccess: probability of successfully eradication before the end of time horizon
#
###############################################################

HazardRateMod <- function(H, cs, DP, Wa, dto, PL, avg_E, tm, se,
                          HRMpars,
                          t_start, t_end, new_E, neighbor_dist){
  
  #debug

  
  # unpack HRM parameters (regression coefficients)
  #  - from Dodd et al. 2015
  a = HRMpars[1]     # intercept
  b = HRMpars[2:10]  # slopes for the predictor variables:
                     #    log(H), cs, DP - m_DP, log(Wa), log(dto), PL - m_PL, avg_E, tm, se
  c = HRMpars[11]    # random effect
  v = HRMpars[12]    # shape parameter of the Weibull distribution
  m_DP = 8.42    # mean detection period (months) of species used in Dodd et al. 2015
  m_PL = 8.12    # mean propagule longevity (years) of species used in Dodd et al. 2015
  
  
  
  # estimate the probability of eradication success 
  #  - using the regression model to predict the probability
  
  if(t_start != t_end){    # if there are still remaining time to do management
    X <- c(log(H), cs, DP - m_DP, log(Wa), log(dto), PL - m_PL, avg_E, tm, se)  # predicators of the regression model
    r <- exp(a + b%*%X + c)           # r: the reciprocal of the mean time to extirpation
    lambda <- c(gamma(1 + 1/v)*r)     # paramter in modelling mean time to extirpation using Weilbull distribution
    
    #debug
    # x <- c(1:20)
    # y <- 1 - exp(-lambda*x^v)# CDF x<10
    # plot(x,y, ylim = c(0,1))
 
    if(new_E==0|new_E==neighbor_dist){
      pSuccess <- 0          # if no survey/treatment is allocated to this project in the future, 
                             # the weed species will definitely persist (even there are past efforts)
    }else{
      pSuccess <-  1 - exp(-lambda*(t_end-t_start)^v)   # cumulative probability distribution of Weibull distribution
                                                 # --> probability successfully eradicate the species before the end of time horizon
    }
    if(H == 0){     # if the weed species has already been eradicated, prob set to zero
      pSuccess <- 0
    }
  
  }else if(t_start==t_end){  # if already at the terminal time, there is no chance of eradicating the species
    pSuccess <- 0
  }

  return(pSuccess)
}

