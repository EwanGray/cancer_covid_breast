#V2.2 For Breast Cancer, Screen-detected and Symptomatic comparison, PSA included

library(tidyr)
library(dplyr)
library(formattable)
library(flexsurv)  # log-logistic and gompertz distributions
library(ggplot2)   # plotting survival curves
library(gtools)

##
##  This code was adapted from code written to perform the analyses presented in the manuscript:
##
##  An inverse stage-shift model to estimate the excess mortality and health 
##  economic impact of delayed access to cancer services due to the COVID-19 pandemic
##
##  by Koen Degeling, Nancy N Baxter, Jon Emery, Fanny Franchini, Peter Gibbs,
##  G Bruce Mann, Grant McArthur, Benjamin J Solomon, Maarten J IJzerman
##
##  For more information or when using/adapting the code, please refer to/cite the
##  pre-print version of the manuscript available from MedRxiv:
##  https://www.medrxiv.org/content/10.1101/2020.05.30.20117630v1
##  Message from original authors:
##  This code was written by Dr Koen Degeling, Cancer Health Services Research,
##  Centre for Cancer Research & Centre for Health Policy, Faculty of Medicine,
##  Dentistry and Health Sciences, The University of Melbourne, Melbourne, Australia
##
##  The code was reviewed by Dr Fanny Franchini, Cancer Health Services Research,
##  Centre for Cancer Research & Centre for Health Policy, Faculty of Medicine,
##  Dentistry and Health Sciences, The University of Melbourne, Melbourne, Australia
##
##  For questions or other enquiries, please contact Koen Degeling by email:
##  koen.degeling@unimelb.edu.au
##

### 1. INITIALIZATION ----

# Clear environment
# rm(list = ls()); gc();

# Loading packages
#library(dplyr);     # handling data.frames
#library(flexsurv);  # log-logistic and gompertz distributions
#library(ggplot2);   # plotting survival curves

# Global parameters
days_in_year <- 365.25;

### 2. FUNCTIONS TO PERFORM THE ANALYSIS ----

# In this section custom functions are defined to perform the analysis. 
# The 'fun_fit_distrs' function fits parameteric distributions based on quantiles and 
# corresponding percentiles.
# The 'fun_plot_distrs' function plots the survival data and fitted distributions for 
# graphical inspection. 
# The 'fun_integrate_distrs' function determines the expected survival in life years 
# of a certain time horizon for a distribution fitted by the 'fun_fit_distrs' function 
# by integrating the corresponding survival curve.


fun_fit_distrs <- function(q, p) {
  
  # This function fits a series of parametric distributions to quantiles and corresponding 
  # percentiles. It does so through nonlinear least squares regression on the parameters on 
  # log-scale. It fits the exponential, Gamma, Gompertz, log-logistic, log-Normal and
  # Weibull distributions. It return the fitted distribution parameters and sum of the
  # squared residuals in a list.
  #
  # INPUTS:
  # - p           vector of percentiles for which q provided quantiles
  # - q           vector of quantiles corresponding to the percentiles provided to p
  #
  # USES:
  # - pgompertz   function for the Gompertz distribution from the flexsurv package
  # - pllogis     function for the log-logistic distribution from the flexsurv package
  #
  # OUTPUTS:
  # - list        a list with two named objects:
  #               - data      
  #               - distrs    list named according to the distribution names: exponential, 
  #                           gamma, gompertz, loglogistic, lognormal, and weibull. Each 
  #                           list item contains a list of three ojects:
  #                           - pars          the estimated parameters on log scale
  #                           - convergence   logical indicating convergence
  #                           - SSM           the sum of the squared residuals
  
  # Merge the two vectors 'p' and 'q' into a data.frame for use in the call of the 
  # nonlinear least squares regression. Also sort
  df_fit <- data.frame(p = p, q = q);
  
  # The start value for the exponential distribution is estimated based on the latest
  # survival data point. To easily extract that data point, we sort the data first.
  df_fit <- arrange(df_fit, q);
  
  # Estimating the mortality rate according to the standard rate-probability formula
  start_q <- tail(df_fit, 1)$q;
  start_p <- tail(df_fit, 1)$p;
  start_rate <- -(1/start_q)*log(start_p);
  
  # The estimated rate is used as start value for the regression. Note that parameters
  # are estimated on log scale.
  fit_exponential <- nls(formula = p ~ pexp(q = q, rate = exp(log_rate), lower.tail = FALSE),
                         start = list(log_rate = log(start_rate)),
                         data = df_fit); 
  
  # The Weibull distribution is fitted subsequently, because start values for its
  # parameters can be defined based on the rate of the exponential distribution and its
  # parameters themselves are useful to define start values for the other distributions.
  # The start values of the Weibull distribution used basically resemble the fitted
  # exponential distribution.
  fit_weibull <- nls(formula = p ~ pweibull(q = q, shape = exp(log_shape), scale = exp(log_scale), lower.tail = FALSE),
                     start = list(log_shape = log(1), log_scale = -coef(fit_exponential)["log_rate"]),
                     data = df_fit);
  
  # The rate of exponential distribution is used to define the start value for the rate
  # parameter of the Gamma distribution, and the start value for the shape parameters is
  # based on that of the Weibull distribution.
  fit_gamma <- nls(formula = p ~ pgamma(q = q, shape = exp(log_shape), rate = exp(log_rate), lower.tail = FALSE),
                   start = list(log_shape = coef(fit_weibull)["log_shape"], log_rate = coef(fit_exponential)["log_rate"]),
                   data = df_fit);
  
  # The Gompertz distribution is the most difficult one to fit. Because it is most similar
  # to the Gamma distribution, the fitted parameters of the Gamma distribution are used as
  # start values. If fitting the Gompertz is unsuccessful using built-in function, it is
  # performed in a custom way using optim directly.
  fit_gompertz <- tryCatch({nls(formula = p ~ pgompertz(q = q, shape = exp(log_shape), rate = exp(log_rate), lower.tail = FALSE),
                                start = list(log_shape = coef(fit_gamma)["log_shape"], log_rate = coef(fit_gamma)["log_rate"]),
                                data = df_fit)}, error = function(e) NULL);
  if(is.null(fit_gompertz)) {
    fit_gompertz <- optim(
      par = c(coef(fit_gamma)["log_shape"], coef(fit_gamma)["log_rate"]),
      fn = function(x, p, q) sum((p - pgompertz(q = q, shape = exp(x["log_shape"]), rate = exp(x["log_rate"]), lower.tail = FALSE))^2),
      p = df_fit$p,
      q = df_fit$q);
  }
  
  # For the log-logistic distribution, start values for the parameters are based on the 
  # parameters of the Weibull distribution.
  fit_loglogistic <- nls(formula = p ~ pllogis(q = q, shape = exp(log_shape), scale = exp(log_scale), lower.tail = FALSE),
                         start = list(log_shape = coef(fit_weibull)["log_shape"], log_scale = coef(fit_weibull)["log_scale"]),
                         data = df_fit);
  
  # For the log-Normal distribution, start values for both parameters are based on the scale
  # parameter of the Weibull distribution. Note that this is the only distribution that does
  # not require the parameters to be transformed back to the original scale.
  fit_lognormal <- nls(formula = p ~ plnorm(q = q, meanlog = log_mean, sdlog = log_sd, lower.tail = FALSE),
                       start = list(log_mean = coef(fit_weibull)["log_scale"], log_sd = coef(fit_weibull)["log_scale"]),
                       data = df_fit);
  
  
  
  # Define a list 'list_out' that contains the data used for fitting as argument 'data' 
  # and within argument 'distrs' for each distribution its parameters, convergence and 
  # sum of the squared residuals
  list_out <- list(
    
    # The data.frame used for fitting
    data = df_fit,
    
    # List including all distributions
    distrs = list(
      
      # Exponential
      exponential = list(
        pars = coefficients(fit_exponential),
        convergence = fit_exponential$convInfo$isConv,
        SSR = sum(residuals(fit_exponential)^2)
      ),
      
      # Gamma
      gamma = list(
        pars = coefficients(fit_gamma),
        convergence = fit_gamma$convInfo$isConv,
        SSR = sum(residuals(fit_gamma)^2)
      ),
      
      # Gompertz
      gompertz = if(class(fit_gompertz) == "nls") {
        list(
          pars = coefficients(fit_gompertz),
          convergence = fit_gompertz$convInfo$isConv,
          SSR = sum(residuals(fit_gompertz)^2)
        )
      } else {
        gompertz = list(
          pars = fit_gompertz$par,    
          convergence = (fit_gompertz$convergence == 0),
          SSR = fit_gompertz$value
        )
      },
      
      # Log-normal
      lognormal = list(
        pars = coefficients(fit_lognormal),
        convergence = fit_lognormal$convInfo$isConv,
        SSR = sum(residuals(fit_lognormal)^2)
      ),
      
      # Log-logistic
      loglogistic = list(
        pars = coefficients(fit_loglogistic),
        convergence = fit_loglogistic$convInfo$isConv,
        SSR = sum(residuals(fit_loglogistic)^2)
      ),
      
      # Weibull
      weibull = list(
        pars = coefficients(fit_weibull),
        convergence = fit_weibull$convInfo$isConv,
        SSR = sum(residuals(fit_weibull)^2)
      )
      
    )
  );
  
  # Return the list
  return(list_out);
  
};


fun_plot_distrs <- function(list_out, time_horizon) {
  
  # This function plots the survival data and fitted parametric distributions based on the
  # list_out object returned by the fun_fit_distrs function.
  #
  # INPUTS:
  # - list_out        list returned by the fun_fit_distrs function, containing the data and
  #                   fitted distributions
  # - time_horizon    the time horizon for which the parametric distributions are to be
  #                   plotted
  #
  # USES:
  # - ggplot          plotting functions from the ggplot2 package
  #
  # OUTPUTS:
  # - plot            the survival plot
  
  # Extracting the survival data and distributions from the list_out object.
  df_fit <- list_out$data;
  distrs <- list_out$distrs;
  
  # Setting some general parameters for the plot. For the type and width of the lines, there
  # will be a difference between the actual survival data (obs = observed) and the simulated
  # data according to the parametric distributions (sim = simulated).
  x_values <- seq(from = 0, to = time_horizon, by = 1/12);
  line_width <- c(obs = 1.5, sim = 1);
  line_type <- c(obs = 1, sim = 2);
  colour_blind_palette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7");
  
  # Plot the distributions in order from generally worst to best fitting distributions, so 
  # that the latter are plotted over the former
  ggplot() +
    
    # General plotting parameters
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual(values = colour_blind_palette) + 
    labs(title = "Observed vs. simulated survival",
         x = "Time in Years",
         y = "Survival Probability") +
    theme(legend.title = element_blank()) + 
    theme(axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.title = element_text(size = 14,face = "bold"),
          legend.text = element_text(size = 14),
          plot.title = element_text(size = 14, face = "bold")) + 
    
    # Observed data
    geom_line(aes(x = df_fit$q, 
                  y = df_fit$p, 
                  colour = "Data"),
              size = line_width["obs"],
              linetype = line_type["obs"]) + 
    
    # Exponential
    geom_line(aes(x = x_values, 
                  y = pexp(q = x_values, 
                           rate = exp(distrs$exponential$pars["log_rate"]),
                           lower.tail = FALSE),
                  colour = "Exponential"),
              size = line_width["sim"],
              linetype = line_type["sim"]) + 
    
    # Log-logistic
    geom_line(aes(x = x_values, 
                  y = pllogis(q = x_values, 
                              shape = exp(distrs$loglogistic$pars["log_shape"]),
                              scale = exp(distrs$loglogistic$pars["log_scale"]),
                              lower.tail = FALSE),
                  colour = "Log-logistic"),
              size = line_width["sim"],
              linetype = line_type["sim"]) + 
    
    # Log-normal
    geom_line(aes(x = x_values, 
                  y = plnorm(q = x_values, 
                             meanlog = distrs$lognormal$pars["log_mean"],
                             sdlog = distrs$lognormal$pars["log_sd"],
                             lower.tail = FALSE),
                  colour = "Log-normal"),
              size = line_width["sim"],
              linetype = line_type["sim"]) + 
    
    # Gamma
    geom_line(aes(x = x_values, 
                  y = pgamma(q = x_values, 
                             shape = exp(distrs$gamma$pars["log_shape"]),
                             rate = exp(distrs$gamma$pars["log_rate"]),
                             lower.tail = FALSE),
                  colour = "Gamma"),
              size = line_width["sim"],
              linetype = line_type["sim"]) +
    
    # Weibull
    geom_line(aes(x = x_values, 
                  y = pweibull(q = x_values, 
                               shape = exp(distrs$weibull$pars["log_shape"]),
                               scale = exp(distrs$weibull$pars["log_scale"]),
                               lower.tail = FALSE),
                  colour = "Weibull"),
              size = line_width["sim"],
              linetype = line_type["sim"]) +  
    
    # Gompertz
    geom_line(aes(x = x_values, 
                  y = pgompertz(q = x_values, 
                                shape = exp(distrs$gompertz$pars["log_shape"]),
                                rate = exp(distrs$gompertz$pars["log_rate"]),
                                lower.tail = FALSE),
                  colour = "Gompertz"),
              size = line_width["sim"],
              linetype = line_type["sim"]);
  
};


fun_integrate_dist <- function(dist, log_pars, time_horizon) { ## Koen, I have edited the function goal here
  
  # This function returns the expected survival in life years over a chosen time horizon 
  # for a distribution fitted using the 'fun_fit_distrs' function, by integrating the 
  # corresponding survival curve.
  #
  # //delete// This function plots the survival data and fitted parametric distributions based on the
  # list_out objected returned by the fun_fit_distrs function. //delete//
  #
  # INPUTS:
  # - dist            character identifying the distribution type
  # - log_pars        vector with the parameters of the distribution on log scale
  # - time_horizon    the time horizon for which the expected survival is to be calculated
  #
  # USES:
  # - pgompertz       function for the Gompertz distribution from the flexsurv package
  # - pllogis         function for the log-logistic distribution from the flexsurv package
  #
  # OUTPUTS:
  # - out             expected survival in life years over the time horizon 
  
  # Initialize the out object to NULL so it can be easily be checked whether the integration
  # was successful.
  out <- NULL;
  
  # Exponential
  if(dist == "exponential") {
    out <- integrate(
      f = pexp,
      lower = 0,
      upper = time_horizon,
      lower.tail = FALSE,
      rate = exp(log_pars["log_rate"])
    )$value;
  }
  
  # Gamma
  if(dist == "gamma") {
    out <- integrate(
      f = pgamma,
      lower = 0,
      upper = time_horizon,
      lower.tail = FALSE,
      shape = exp(log_pars["log_shape"]),
      rate = exp(log_pars["log_rate"])
    )$value;
  }
  
  # Gompertz
  if(dist == "gompertz") {
    out <- integrate(
      f = pgompertz,
      lower = 0,
      upper = time_horizon,
      lower.tail = FALSE,
      shape = exp(log_pars["log_shape"]),
      rate = exp(log_pars["log_rate"])
    )$value;
  }
  
  # Log-logistic
  if(dist == "loglogistic") {
    out <- integrate(
      f = pllogis,
      lower = 0,
      upper = time_horizon,
      lower.tail = FALSE,
      shape = exp(log_pars["log_shape"]),
      scale = exp(log_pars["log_scale"])
    )$value;
  }
  
  # Log-normal
  if(dist == "lognormal") {
    out <- integrate(
      f = plnorm,
      lower = 0,
      upper = time_horizon,
      lower.tail = FALSE,
      meanlog = log_pars["log_mean"],
      sdlog = log_pars["log_sd"]
    )$value;
  }
  
  # Weibull
  if(dist == "weibull") {
    out <- integrate(
      f = pweibull,
      lower = 0,
      upper = time_horizon,
      lower.tail = FALSE,
      shape = exp(log_pars["log_shape"]),
      scale = exp(log_pars["log_scale"])
    )$value;
  }
  
  # Check whether integration successful
  if(is.null(out)) stop("No integration has been performed. Perhaps the distribution is not supported?");
  
  # The time_horizon is included in the name of the object for convenience
  names(out) <- paste0("time_horizon_", time_horizon);
  
  # Return the output
  return(out);
  
};


####STAGE SHIFT MODELLING DEFINED AS A FUNCTION####

###BASE CASE PARAMETER VALUES###
#For parameters that are being varied
Base_HR_breast_delay_I <- 1.018
# The stage-specific survival data for breast cancer is taken from Scottish registry data using data from diagnoses 2010-2018:
Base_df_survival_breast <- data.frame(
  disease_stage = c("I", "II", "III", "IV", "unknown"),
  survival_0y = c(1.000, 1.000, 1.000, 1.000, 1.000),
  survival_1y = c(0.995, 0.989, 0.962, 0.929, 0.994),
  survival_2y = c(0.986, 0.964, 0.880, 0.789, 0.973),
  survival_3y = c(0.973, 0.935, 0.816, 0.739, 0.938),
  survival_4y = c(0.959, 0.906, 0.769, 0.676, 0.899),
  survival_5y = c(0.943, 0.879, 0.721, 0.65, 0.845)
);
# The predicted incidence for 2020 is also taken from the average of 2017/2018 in Scottish cancer registry:
Base_n_incidence_breast <- 8814/2;
# The distribution of the disease stage at treatment initiation is based on Scottish registry 2017-2018 sum:
Base_v_stage_numbers_breast <- c(I = 3561, II = 4042, III = 664, IV = 450, unknown = 94);


#############stage_shift function##########
##INPUTS: TTI hazard ratio (scalar), stage-specific survival data (data frame 5 years (columns) by 5 stages (rows)), total annual incidence (scalar), incidence for each stage (vector), whether to use within-stage TTI HR (1=use)

##OUTPUTS: List: [1] Table with excess mortality over 5 years, LY lost over 5 years, LY lost over 10 years, for both 3 and 6 month delays, [2] Predicted stage distribution, [3] predicted survival no delay, [4] predicted survival for 3 month delay, [5] predited survival for 3m and 6m delay w TTI HR applied

#Model function:

stage_shift <- function(HR_TTI,df_survival_stage,tot_incidence,stage_numbers, within_stage) {
  # #testing
  # HR_TTI <- Base_HR_breast_delay_I
  # df_survival_stage <- Base_df_survival_breast
  # tot_incidence <- Base_n_incidence_breast
  # stage_numbers <- Base_v_stage_numbers_breast
  
  ### 3. STAGE SHIFT MODELLING ----
  
  ## 3.1 BREAST CANCER ----
  
  ## DATA USED FOR MODELLING
  
  # The stage-specific survival data:
  df_survival_breast <- df_survival_stage;
  
  # The predicted incidence for 2020
  n_incidence_breast <- tot_incidence;
  
  # The distribution of the disease stage at treatment initiation:
  v_stage_numbers_breast <- stage_numbers
  v_stage_proportions_breast <- v_stage_numbers_breast / sum(v_stage_numbers_breast);
  
  # The costs of treating a specific stage of disease have been approximated from the
  # publication by Goldsbury et al. (2018):
  v_costs_breast <- c(I = 50699, II = 67069, III = 98632, IV = 105934, unknown = 100860);
  
  
  ## ESTIMATING TIME TO STAGE PROGRESSION (TTSP)
  
  # For breast cancer, TTSP is estimated by approximating how long a delay needs to be for the
  # 5-year survival of stage I patients to match that of stage II patients based on the hazard
  # ratio of a surgical delay on overall survival.
  
  # The required hazard ratio for a stage I patient to match the 5-year survival of a stage II 
  # patients is calculated as the ratio of the mortality rates. The 5-year survival of 100% 
  # for stage I patients is problematic, because it suggests no events, which does not allow 
  # a rate to be calculated. Therefore, we assume the 5-year survival probability to be 99%.
  #No longer 100% survival for stage I in data so correction no longer needed.
  #v2.1 converted to transitions between all stages, including transitioning up multiple stages, note that
  # transitions to stage 3/4 are handled so that stage 3/4 ratio remains the same. Assumes presentations at 3
  # are less afffected by disruption as more likely to be urgent/emergency.
  
  mortality_rate_breast_I <- -(1/5) * log(filter(df_survival_breast, disease_stage == "I")$survival_5y);
  mortality_rate_breast_II <- -(1/5) * log(filter(df_survival_breast, disease_stage == "II")$survival_5y);
  mortality_rate_breast_III <- -(1/5) * log(filter(df_survival_breast, disease_stage == "III")$survival_5y);
  mortality_rate_breast_IV <- -(1/5) * log(filter(df_survival_breast, disease_stage == "IV")$survival_5y);
  
  HR_breast_required_I_II <- mortality_rate_breast_II / mortality_rate_breast_I;
  HR_breast_required_I_III <- mortality_rate_breast_III / mortality_rate_breast_I;
  HR_breast_required_I_IV <- mortality_rate_breast_IV / mortality_rate_breast_I;
  HR_breast_required_II_III <- mortality_rate_breast_III / mortality_rate_breast_II;
  HR_breast_required_II_IV <- mortality_rate_breast_IV / mortality_rate_breast_II;
  HR_breast_required_III_IV <- mortality_rate_breast_IV / mortality_rate_breast_III;
  
  # Hazard ratio of surgical delays on overall survival for stage I patients:
  # - HR on overall survival is 1.018 for a 1-week delay (Khorana, 2019).
  # Can add different HR for each stage progression.
  HR_breast_delay_I <- HR_TTI;
  HR_breast_delay_II <- HR_TTI;
  HR_breast_delay_III <- HR_TTI;
  
  t_breast_delay_I <- 7 / days_in_year * 12;     # in months
  
  # The HRs are converted to regression coefficients by transferring to log scale to model
  # the coefficient as a linear function of time
  coef_breast_delay_I <- log(HR_breast_delay_I) / t_breast_delay_I;
  coef_breast_delay_II <- log(HR_breast_delay_II) / t_breast_delay_I;
  coef_breast_delay_III <- log(HR_breast_delay_III) / t_breast_delay_I;
  
  coef_breast_required_I_II <- log(HR_breast_required_I_II);
  coef_breast_required_I_III <- log(HR_breast_required_I_III);
  coef_breast_required_I_IV <- log(HR_breast_required_I_IV);
  coef_breast_required_II_III <- log(HR_breast_required_II_III);
  coef_breast_required_II_IV <- log(HR_breast_required_II_IV);
  coef_breast_required_III_IV <- log(HR_breast_required_III_IV);
  
  # TTSP (in months) can now be calculated
  TTSP_breast_I_II <- coef_breast_required_I_II/ coef_breast_delay_I;
  TTSP_breast_I_III <- coef_breast_required_I_III / coef_breast_delay_I;
  TTSP_breast_I_IV <- coef_breast_required_I_IV / coef_breast_delay_I;
  TTSP_breast_II_III <- coef_breast_required_II_III / coef_breast_delay_II;
  TTSP_breast_II_IV <- coef_breast_required_II_IV / coef_breast_delay_II;
  TTSP_breast_III_IV <- coef_breast_required_III_IV / coef_breast_delay_III;
  
  # proportion of stage 3 in 3/4 at baseline
  stage_III_IV_prop <- v_stage_numbers_breast[3]/(v_stage_numbers_breast[3]+v_stage_numbers_breast[4])
  
  ## CALCULATING THE STAGE SHIFT
  
  # Based on the TTSP, the proportion of patients who are expected to start treatment at a
  # more advanced stage can be calculated. Since only the expected TTSP is known, and not its
  # distribution, the only distribution that can be used is the exponential distribution for
  # which the rate parameter is defined as: 1 / expected value.
  
  #Calculate transition probabilities sequentially from smallest to largest
  # i.e. To progress to stage III need to also progress to stage II
  
  # Proportion of patients that will progress following a 3-month delay, adjusted for the fact
  # that the service disruption will only be for 3 months, i.e. 25% of 2020 patients:
  prop_progress_breast_I_IV_3mo <- pexp(q = 3, rate = 1 / TTSP_breast_I_IV)* (3/12);
  prop_progress_breast_I_III_3mo <- pexp(q = 3, rate = 1 / TTSP_breast_I_III)* (3/12) - prop_progress_breast_I_IV_3mo;
  prop_I_3to4 <- (prop_progress_breast_I_IV_3mo + prop_progress_breast_I_III_3mo);
  prop_progress_breast_I_III_3mo <- prop_I_3to4*stage_III_IV_prop
  prop_progress_breast_I_IV_3mo <- prop_I_3to4*(1 - stage_III_IV_prop)
  prop_progress_breast_I_II_3mo <-  pexp(q = 3, rate = 1 / TTSP_breast_I_II) * (3/12) - (prop_progress_breast_I_IV_3mo + prop_progress_breast_I_III_3mo);
  prop_progress_breast_II_IV_3mo <- pexp(q = 3, rate = 1 / TTSP_breast_II_IV)* (3/12);
  prop_progress_breast_II_III_3mo <- pexp(q = 3, rate = 1 / TTSP_breast_II_III)*(3/12) - prop_progress_breast_II_IV_3mo;
  prop_II_3to4 <- (prop_progress_breast_II_IV_3mo + prop_progress_breast_II_III_3mo);
  prop_progress_breast_II_III_3mo <- prop_II_3to4*stage_III_IV_prop
  prop_progress_breast_II_IV_3mo <- prop_II_3to4*(1 - stage_III_IV_prop)
  prop_progress_breast_III_IV_3mo <- 0 #pexp(q = 3, rate = 1 / TTSP_breast_III_IV)* (3/12); assume stage 3 still present
  
  # Proportion of patients that will progress following a 6-month delay, adjusted for the fact
  # that the service disruption will only be for 6 months, i.e. 50% of 2020 patients:
  prop_progress_breast_I_IV_6mo <- pexp(q = 6, rate = 1 / TTSP_breast_I_IV)* (6/12);
  prop_progress_breast_I_III_6mo <- pexp(q = 6, rate = 1 / TTSP_breast_I_III)*(6/12) - prop_progress_breast_I_IV_6mo;
  prop_I_3to4 <- (prop_progress_breast_I_IV_6mo + prop_progress_breast_I_III_6mo);
  prop_progress_breast_I_III_6mo <- prop_I_3to4*stage_III_IV_prop
  prop_progress_breast_I_IV_6mo <- prop_I_3to4*(1 - stage_III_IV_prop)
  prop_progress_breast_I_II_6mo <-  pexp(q = 6, rate = 1 / TTSP_breast_I_II) * (6/12) - (prop_progress_breast_I_IV_6mo + prop_progress_breast_I_III_6mo);
  prop_progress_breast_II_IV_6mo <- pexp(q = 6, rate = 1 / TTSP_breast_II_IV)* (6/12);
  prop_progress_breast_II_III_6mo <- pexp(q = 6, rate = 1 / TTSP_breast_II_III)*(6/12) - prop_progress_breast_II_IV_6mo;
  prop_II_3to4 <- (prop_progress_breast_II_IV_6mo + prop_progress_breast_II_III_6mo);
  prop_progress_breast_II_III_6mo <- prop_II_3to4*stage_III_IV_prop
  prop_progress_breast_II_IV_6mo <- prop_II_3to4*(1 - stage_III_IV_prop)
  prop_progress_breast_III_IV_6mo <- 0 #pexp(q = 6, rate = 1 / TTSP_breast_III_IV)* (6/12); assume stage 3 still present
  
  
  # Based on these proportions, the distribution of disease stage at treatment initiation
  # following a 3 month or 6 month service disruption and delay can be calculated:
  v_stage_proportions_breast_3mo <- c(
    I = v_stage_proportions_breast["I"] * (1 - (prop_progress_breast_I_II_3mo + prop_progress_breast_I_III_3mo + prop_progress_breast_I_IV_3mo)),
    II = v_stage_proportions_breast["II"] + (v_stage_proportions_breast["I"] * prop_progress_breast_I_II_3mo) - (v_stage_proportions_breast["II"] * prop_progress_breast_II_III_3mo) - (v_stage_proportions_breast["II"] * prop_progress_breast_II_IV_3mo),
    III = v_stage_proportions_breast["III"] + (v_stage_proportions_breast["I"] * prop_progress_breast_I_III_3mo) + (v_stage_proportions_breast["II"] * prop_progress_breast_II_III_3mo) - (v_stage_proportions_breast["III"] * prop_progress_breast_III_IV_3mo),
    IV = v_stage_proportions_breast["IV"] + (v_stage_proportions_breast["I"] * prop_progress_breast_I_IV_3mo) + (v_stage_proportions_breast["II"] * prop_progress_breast_II_IV_3mo) + (v_stage_proportions_breast["III"] * prop_progress_breast_III_IV_3mo),
    unknown = v_stage_proportions_breast["unknown"]
  );
  
  v_stage_proportions_breast_6mo <- c(
    I = v_stage_proportions_breast["I"] * (1 - (prop_progress_breast_I_II_6mo + prop_progress_breast_I_III_6mo + prop_progress_breast_I_IV_6mo)),
    II = v_stage_proportions_breast["II"] + (v_stage_proportions_breast["I"] * prop_progress_breast_I_II_6mo) - (v_stage_proportions_breast["II"] * prop_progress_breast_II_III_6mo) - (v_stage_proportions_breast["II"] * prop_progress_breast_II_IV_6mo),
    III = v_stage_proportions_breast["III"] + (v_stage_proportions_breast["I"] * prop_progress_breast_I_III_6mo) + (v_stage_proportions_breast["II"] * prop_progress_breast_II_III_6mo) - (v_stage_proportions_breast["III"] * prop_progress_breast_III_IV_6mo),
    IV = v_stage_proportions_breast["IV"] + (v_stage_proportions_breast["I"] * prop_progress_breast_I_IV_6mo) + (v_stage_proportions_breast["II"] * prop_progress_breast_II_IV_6mo) + (v_stage_proportions_breast["III"] * prop_progress_breast_III_IV_6mo),
    unknown = v_stage_proportions_breast["unknown"]
  );
  
  #Add additional delay hazard - assumes proportional hazards
  # - determine if this is to be used based on function input
  if(within_stage != 1){HR_breast_delay_I <- 1}
  #HR for 3 & 6 month delay periods (12 and 24 weeks)
  HR_breast_delay_3m <- exp(log(HR_breast_delay_I)*(3*4.345))
  HR_breast_delay_6m <- exp(log(HR_breast_delay_I)*(6*4.345))
  #Calculate the survival probabilities with this delay added
  #3 month
  df_survival_breast_delay_3m <- df_survival_breast %>%
    select(contains("survival")) # select columns
  df_survival_breast_delay_3m <- df_survival_breast_delay_3m ^ HR_breast_delay_3m #apply HR
  df_survival_breast_delay_3m <- bind_cols(select(df_survival_breast,contains("stage")),df_survival_breast_delay_3m) # remake matrix
  #6 month
  df_survival_breast_delay_6m <- df_survival_breast %>%
    select(contains("survival")) # select columns
  df_survival_breast_delay_6m <- df_survival_breast_delay_6m ^ HR_breast_delay_6m #apply HR
  df_survival_breast_delay_6m <- bind_cols(select(df_survival_breast,contains("stage")),df_survival_breast_delay_6m) # remake matrix
  
  # Checking whether the distribution are still appropriate, i.e. sum up to 1
  #sum(v_stage_proportions_breast_3mo);
  #sum(v_stage_proportions_breast_6mo);
  
  
  ## ESTIMATING THE IMPACT OF THE STAGE SHIFTS
  
  # According to each distribution of stages at treatment initiation, the combined survival
  # is calculated by weighting the year-specific probabilities
  
  # Extracting the survival data from the data.frame with stage-specific data
  m_survival_breast <- select(df_survival_breast, contains("survival"));
  m_survival_breast_3m <- select(df_survival_breast_delay_3m, contains("survival"));
  m_survival_breast_6m <- select(df_survival_breast_delay_6m, contains("survival"));
  
  # Calculating the weighted survival
  v_survival_breast_base <- apply(m_survival_breast, 2, function(col) sum(v_stage_proportions_breast * col));
  v_survival_breast_3mo <- apply(m_survival_breast_3m, 2, function(col) sum(v_stage_proportions_breast_3mo * col));
  v_survival_breast_6mo <- apply(m_survival_breast_6m, 2, function(col) sum(v_stage_proportions_breast_6mo * col));
  
  # To allow for extrapolation and integration to calculate expected survival in terms of 
  # life years, a series of distributions are fitted to the data: exponential, Gamma, Gompertz
  # log-Normal, log-logistic, and Weibull distributions. These are fitted based on the quantiles
  # and percentiles of the weighted survival using the custom fun_fit_distrs function. The most
  # appropriate distribution is selected by graphical inspection using the 'fun_plot_distrs'
  # function and based on the sum of the squared residuals of the non-linear least squares 
  # regression models used to estimate the parameters on log scale.
  
  # Define a factor with the quantiles corresponding to the weighted percentiles
  v_q <- 0:5;
  
  # Selecting the distribution: baseline
  # - Based on that beyond 5 years background mortality is more likely than cancer-specific
  #   mortality, and Gompertz distributions are known to represent background mortality well,
  #   the Gompertz is selected, providing the most conservative estimate of survival.
  l_distr_breast_base <- fun_fit_distrs(q = v_q, p = v_survival_breast_base); 
  fun_plot_distrs(list_out = l_distr_breast_base, time_horizon = 10);
  sapply(l_distr_breast_base$distrs, function(l) l$SSR);
  dist_breast_base <- "exponential";
  
  # Selecting the distribution: 3 month disruption and delay
  # - Given that the weighted survival for the 3 month disruption and delay scenario is almost
  #   identical to that of the baseline scenario, and no differences between survival due to
  #   different parametric distributions should be introduced, the Gompertz is selected.
  l_distr_breast_3mo <- fun_fit_distrs(q = v_q, p = v_survival_breast_3mo); 
  fun_plot_distrs(list_out = l_distr_breast_3mo, time_horizon = 10);
  sapply(l_distr_breast_3mo$distrs, function(l) l$SSR);
  dist_breast_3mo <- "exponential";
  
  # Selecting the distribution: 6 month disruption and delay
  # - Same as for the 3 month disruption and delay
  l_distr_breast_6mo <- fun_fit_distrs(q = v_q, p = v_survival_breast_6mo); 
  fun_plot_distrs(list_out = l_distr_breast_6mo, time_horizon = 10);
  sapply(l_distr_breast_6mo$distrs, function(l) l$SSR);
  dist_breast_6mo <- "exponential";
  
  # Obtaining the expected survival in life years for two time horizons: 5 and 10 years
  # - Because the built-in integrate function used within the fun_integrate_dist function
  #   is not vectorized by default, a loop approach is taken using sapply. Alternatively,
  #   the fun_integrate_dist function can be vectorized
  v_time_horizons <- c(5, 10);
  
  v_ly_breast_base <- sapply(v_time_horizons, function(time_horizon) {
    fun_integrate_dist(dist = dist_breast_base,
                       log_pars = l_distr_breast_base$distrs[[dist_breast_base]]$pars,
                       time_horizon = time_horizon)
  });
  
  v_ly_breast_3mo <- sapply(v_time_horizons, function(time_horizon) {
    fun_integrate_dist(dist = dist_breast_3mo,
                       log_pars = l_distr_breast_3mo$distrs[[dist_breast_3mo]]$pars,
                       time_horizon = time_horizon)
  });
  
  v_ly_breast_6mo <- sapply(v_time_horizons, function(time_horizon) {
    fun_integrate_dist(dist = dist_breast_6mo,
                       log_pars = l_distr_breast_6mo$distrs[[dist_breast_6mo]]$pars,
                       time_horizon = time_horizon)
  });
  
  # The expected costs of treating the 2020 population are estimated by weighting the
  # stage-specific healthcare costs by the distribution of disease stage at treatment
  # initatiation according to the different scenarios
  # - No differences in costs because we only consider stage I -> stage II shifts and
  #   the treatment costs are expected to be similar
  costs_breast_base <- sum(v_stage_proportions_breast * v_costs_breast);
  costs_breast_3mo <- sum(v_stage_proportions_breast_3mo * v_costs_breast);
  costs_breast_6mo <- sum(v_stage_proportions_breast_6mo * v_costs_breast);
  
  
  ## Below extracts results for original manuscript and estiamtes costs, not used in current version
  
  # # Time to stage progression (in years)
  # TTSP_breast_I / 12;
  # 
  # # Probability of stage I patient with a delay to progress to the subsequent stage
  # prob_progress_breast_I_3mo;
  # prob_progress_breast_I_6mo;
  # 
  # # Proportion of stage I patients from the 2020 population that will progress
  # prop_progress_breast_I_3mo;
  # prop_progress_breast_I_6mo;
  # 
  # # Excess mortality for the 2020 population is calculated as the difference in survival
  # # of a delay scenario compared to the baseline, multiplied by the predicted incidence
  # sum((v_survival_breast_base - v_survival_breast_3mo) * n_incidence_breast);
  # sum((v_survival_breast_base - v_survival_breast_6mo) * n_incidence_breast);
  # 
  # # Expected survival in life years PER PATIENT for the considered time horizons
  # v_ly_breast_base;
  # v_ly_breast_3mo;
  # v_ly_breast_6mo;
  # 
  # # Expected survival in life years FOR THE 2020 POPULATION for the considered time horizons
  # v_ly_breast_base * n_incidence_breast;
  # v_ly_breast_3mo * n_incidence_breast;
  # v_ly_breast_6mo * n_incidence_breast;
  # 
  # # Expected healthcare costs PER PATIENT
  # costs_breast_base;
  # costs_breast_3mo;
  # costs_breast_6mo;
  # 
  # # Expected healthcare costs FOR THE 2020 POPULATION (in millions)
  # costs_breast_base * n_incidence_breast / 10^6;
  # costs_breast_3mo * n_incidence_breast / 10^6;
  # costs_breast_6mo * n_incidence_breast / 10^6;
  
  #Breast only for now
  #(1 - v_survival_breast_base[6])*n_incidence_breast # Total deaths from annual cohort
  
  ### MAKE OUTPUT TABLE###
  #Excess mortality, LY lost over 5 year horizon & 10 year horizon, 3 month and 6 month delays from disruption
  Output_table <- data.frame(
    Duration = c("3 months","6 months"),
    excess_mortality = c((v_survival_breast_base[6] - v_survival_breast_3mo[6]) * n_incidence_breast,(v_survival_breast_base[6] - v_survival_breast_6mo[6]) * n_incidence_breast),
    life_years_lost_5 = c(((v_ly_breast_base[1] - v_ly_breast_3mo[1]) * n_incidence_breast),((v_ly_breast_base[1] - v_ly_breast_6mo[1]) * n_incidence_breast)),
    life_years_lost_10 = c(((v_ly_breast_base[2] - v_ly_breast_3mo[2]) * n_incidence_breast),((v_ly_breast_base[2] - v_ly_breast_6mo[2]) * n_incidence_breast))
  )
  
  #Table of stage distributions
  stage_dist_out <- as.data.frame(rbind(v_stage_proportions_breast, v_stage_proportions_breast_3mo,v_stage_proportions_breast_6mo))
  
  #List containing all outputs
  Output_ls <- list(Output_table,stage_dist_out,df_survival_breast,df_survival_breast_delay_3m,df_survival_breast_delay_6m)
  
  return(Output_ls)
}

Output_1 <- stage_shift(Base_HR_breast_delay_I,Base_df_survival_breast,Base_n_incidence_breast,Base_v_stage_numbers_breast,0)
stage_dist_table <- Output_1[[2]]*100
row.names(stage_dist_table) <- c("No delay","3 months","6 months")
formattable(stage_dist_table, align=c("c","c"),digits = 3)


#Screen detected survival data for stage I to III, all breast cancers for IV and unknown with Data from Scottish registry 2010-2017
screen_df_survival_breast <- data.frame(
  disease_stage = c("I", "II", "III", "IV", "unknown"),
  survival_0y = c(1.000, 1.000, 1.000, 1.000, 1.000),
  survival_1y = c(0.997, 0.994, 0.980, 0.929, 0.994),
  survival_2y = c(0.992, 0.984, 0.948, 0.789, 0.973),
  survival_3y = c(0.984, 0.970, 0.913, 0.739, 0.938),
  survival_4y = c(0.977, 0.956, 0.884, 0.676, 0.899),
  survival_5y = c(0.965, 0.933, 0.839, 0.65, 0.845)
)

#Symptomatic/clinical detected survival data for stage I to III, all breast cancers for IV and unknown with Data from Scottish registry 2010-2017
clinical_df_survival_breast <- data.frame(
  disease_stage = c("I", "II", "III", "IV", "unknown"),
  survival_0y = c(1.000, 1.000, 1.000, 1.000, 1.000),
  survival_1y = c(0.993, 0.986, 0.959, 0.929, 0.994),
  survival_2y = c(0.978, 0.957, 0.868, 0.789, 0.973),
  survival_3y = c(0.960, 0.929, 0.799, 0.739, 0.938),
  survival_4y = c(0.936, 0.888, 0.748, 0.676, 0.899),
  survival_5y = c(0.915, 0.859, 0.698, 0.65, 0.845)
)
#Incidence of screen detected & clinically detected breast cancers by stage, Scottish registry 2017
screen_detected_incidence <- 1384
clinical_detected_incidence <- 2016
screen_fraction <- screen_detected_incidence/(screen_detected_incidence +clinical_detected_incidence)

screen_v_stage_numbers_breast <- c(I = 911, II = 354, III = 36, IV = 3, unknown = 13);
clinical_v_stage_numbers_breast <- c(I = 667, II = 1016, III = 268, IV = 6, unknown = 33)

#Get output table for each of screen detected and clincal detected
Output_screen <- stage_shift(Base_HR_breast_delay_I,screen_df_survival_breast,Base_n_incidence_breast*screen_fraction,screen_v_stage_numbers_breast,0)
formattable(Output_screen[[1]],col.names = c("Duration of disruption","Excess mortality","LY lost over 5 years","LY lost over 10 years"), align=c("l","c"), digits = 3)


Output_clinical <- stage_shift(Base_HR_breast_delay_I,clinical_df_survival_breast,Base_n_incidence_breast*(1- screen_fraction),clinical_v_stage_numbers_breast,0)
formattable(Output_clinical[[1]],col.names = c("Duration of disruption","Excess mortality","LY lost over 5 years","LY lost over 10 years"), align=c("l","c"), digits = 3)


Output_combined <- (Output_screen[[1]][,2:4]) + (Output_clinical[[1]][,2:4])
Output_combined <- bind_cols(select(Output_screen[[1]],contains("Duration")),Output_combined)

formattable(Output_combined,col.names = c("Duration of disruption","Excess mortality","LY lost over 5 years","LY lost over 10 years"), align=c("l","c"), digits = 3)

Output_combined_2 <- (Output_screen[[1]][,2:4]*screen_fraction)
Output_combined_2 <- bind_cols(select(Output_screen[[1]],contains("Duration")),Output_combined_2)

formattable(Output_combined_2,col.names = c("Duration of disruption","Excess mortality","LY lost over 5 years","LY lost over 10 years"), align=c("l","c"), digits = 3)


###############PSA###############################
#Get 95% CI for excess mortality estimates accounting for sampling uncertainty in TTI HR and stage specific incidence

#Generate random draw of ITT HR & stage distribution
#TTI_HR - using standard errors from Khorana 2019
#1.018 (1.015-1.020) [HR for stage 1 - used for all transitions]

#RNG
set.seed(100) # seed for replication
n_rep <- 10000 #number of PSA iterations
TTI_HR_ran <-  exp(rnorm(n = n_rep,mean = log(1.018),sd = 0.001)) #values of TTI to use
v_stage_numbers_breast_ran <- rdirichlet(n_rep, Base_v_stage_numbers_breast) #values of stage distribution proportions
v_stage_numbers_breast_ran <- as.data.frame(round(v_stage_numbers_breast_ran*Base_n_incidence_breast,digits = 0)) %>%
  rename( "I" =V1, "II" = V2, "III" = V3, "IV" = V4, "unknown" = V5) #number not proportion

#Data frame to store results
PSA_output <- as.data.frame(matrix(ncol = 2,nrow = length(TTI_HR_ran),data = 0)) %>%
  rename( "3_months" =V1, "6_months" = V2)

#Screen-detected/clinical stratified analysis
#RNG incidence by strata
screen_v_stage_numbers_breast_ran <- rdirichlet(n_rep, screen_v_stage_numbers_breast) #values of stage distribution proportions screen-detected
clin_v_stage_numbers_breast_ran <- rdirichlet(n_rep, clinical_v_stage_numbers_breast) #values of stage distribution proportions clinical-detected

screen_v_stage_numbers_breast_ran  <- as.data.frame(round(screen_v_stage_numbers_breast_ran *Base_n_incidence_breast*screen_fraction,digits = 0)) %>%
  rename( "I" =V1, "II" = V2, "III" = V3, "IV" = V4, "unknown" = V5) #number not proportion
clin_v_stage_numbers_breast_ran  <- as.data.frame(round(clin_v_stage_numbers_breast_ran *Base_n_incidence_breast*(1-screen_fraction),digits = 0)) %>%
  rename( "I" =V1, "II" = V2, "III" = V3, "IV" = V4, "unknown" = V5) #number not proportion

#Data frame to store results
PSA_output_screen <- as.data.frame(matrix(ncol = 6,nrow = length(TTI_HR_ran),data = 0)) %>%
  rename( "3_months_screen" =V1, "6_months_screen" = V2, "3_months_clinical" = V3, "6_months_clinical" = V4, "3_months_comb" = V5, "6_months_comb" = V6)

#Loop to run model at random values and keep results
j <- 1 #counter
for (x in 1:length(TTI_HR_ran)){
  
  iter_v_stage_numbers_breast_screen <- c(I = screen_v_stage_numbers_breast_ran[x,1], II = screen_v_stage_numbers_breast_ran[x,2], III = screen_v_stage_numbers_breast_ran[x,3], IV = screen_v_stage_numbers_breast_ran[x,4], unknown = screen_v_stage_numbers_breast_ran[x,5])
  
  iter_v_stage_numbers_breast_clin <- c(I = clin_v_stage_numbers_breast_ran[x,1], II = clin_v_stage_numbers_breast_ran[x,2], III = clin_v_stage_numbers_breast_ran[x,3], IV = clin_v_stage_numbers_breast_ran[x,4], unknown = clin_v_stage_numbers_breast_ran[x,5])
  
  Output_screen <- stage_shift(TTI_HR_ran[x],screen_df_survival_breast,Base_n_incidence_breast * screen_fraction,iter_v_stage_numbers_breast_screen,0)
  Mortality_tab <- Output_screen[[1]]
  PSA_output_screen[j,1] <- c(Mortality_tab$excess_mortality[1]) 
  PSA_output_screen[j,2] <- c(Mortality_tab$excess_mortality[2])
  
  Output_clinical <- stage_shift(TTI_HR_ran[x],clinical_df_survival_breast,Base_n_incidence_breast * (1-screen_fraction),iter_v_stage_numbers_breast_clin,0)
  Mortality_tab <- Output_clinical[[1]]
  PSA_output_screen[j,3] <- c(Mortality_tab$excess_mortality[1]) 
  PSA_output_screen[j,4] <- c(Mortality_tab$excess_mortality[2])

Output_combined <- (Output_screen[[1]][,2:4]) + (Output_clinical[[1]][,2:4])
PSA_output_screen[j,5] <- c(Output_combined$excess_mortality[1]) 
PSA_output_screen[j,6] <- c(Output_combined$excess_mortality[2])
  
j <- j + 1  
  }

#Get 95% interval for excess deaths over 5 years for 3 & 6 month delays
quantile(PSA_output_screen[,1],probs = c(0.025,0.975))
quantile(PSA_output_screen[,2],probs = c(0.025,0.975))
quantile(PSA_output_screen[,3],probs = c(0.025,0.975))
quantile(PSA_output_screen[,4],probs = c(0.025,0.975))
quantile(PSA_output_screen[,5],probs = c(0.025,0.975))
quantile(PSA_output_screen[,6],probs = c(0.025,0.975))
