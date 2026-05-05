#################################################################################
#--------------------------------------------------------------------------------
# While the parameters are known when generating the data, allowing us to 
# determine the skewness of the distribution, in real-world situations these 
# parameters, and therefore the skewness, are unknown.
# Consequently, the skewness must be estimated using a nonparametric approach, 
# such as the following:
#--------------------------------------------------------------------------------
#################################################################################

Calc_Skew <- function(times, cens, plot = TRUE, quantiles = c(0.25, 0.5, 0.75)){
  
  library(GofCens)
  library(ggplot2)
  library(survival)
  
  ind <- c()
  km_cens <- data.frame()
  min_surv <- numeric()
  sk_est <- numeric()
  
  datcens_fit <- survfit(Surv(times, cens) ~ 1)
  
  km_cens <- data.frame("times" = datcens_fit$time,
                        "surv" = datcens_fit$surv)
  
  if(plot){
    p1 <- ggplot(km_cens, aes(x = times, y = surv)) +
      geom_step(linewidth = 1, colour = "blue") +
      labs(title = "Kaplan-Meier Survival Curve",
           x = "Time",
           y = "Survival probability") +
      theme_minimal() + ylim(c(0,1))
    print(p1)
  }
  
  min_surv <- min(datcens_fit$surv)
  
  if(min_surv > quantiles[1]){
    q_ind <- sapply(c(round(1 - min_surv,2), 0.5, round(min_surv,2)), function(p) {
      ind <- which(abs(datcens_fit$surv - p) < 1e-2)
      if(length(ind) == 0) NA else ind[1]
    })
      }else{
    q_ind <- sapply(quantiles, function(p) {
      ind <- which(abs(datcens_fit$surv - p) < 1e-2)
      if(length(ind) == 0) NA else ind[1]
    })
  }
  
  q <- c(ifelse(length(q_ind[1])==0, NA, datcens_fit$time[q_ind[1]]),
         ifelse(length(q_ind[2])==0, NA, datcens_fit$time[q_ind[2]]),
         ifelse(length(q_ind[3])==0, NA, datcens_fit$time[q_ind[3]]))
  
  sk_est <- as.numeric((q[1] - 2 * q[2] + q[3]) / (q[3] - q[1]))
  
  return(sk_est)
  
}



#-------------------------------------------------------------------------------
# Symmetric Random Censoring Example
#-------------------------------------------------------------------------------

set.seed(28657)

dat <- rweibull(243, shape = 3.5, scale = 20)
cens <- rexp(243, 1/18.5)
tobs <- pmin(dat, cens)
delta <- as.numeric(dat <= cens)
datcens <- data.frame(times = tobs, cens = delta)
datnocens <- data.frame(times = dat, cens = rep(1, n=243))


Calc_Skew(times = datcens$times, cens = datcens$cens, plot = TRUE)


datcens_fit <- survfit(Surv(datcens$times, datcens$cens) ~ 1)

km_cens <- data.frame("times" = datcens_fit$time,
                      "surv" = datcens_fit$surv)

(q <- quantile(datnocens$times, c(0.25,0.5,0.75)))

as.numeric((q[1] - 2 * q[2] + q[3]) / (q[3] - q[1]))

ind1 <- which(round(datcens_fit$surv,2) == 0.75)
ind2 <- which(round(datcens_fit$surv,2) == 0.5)
ind3 <- which(round(datcens_fit$surv,2) == 0.25)
time_q1_obs <- ifelse(length(ind1)==0, NA, datcens_fit$time[ind1])
time_q2_obs <- ifelse(length(ind2)==0, NA, datcens_fit$time[ind2])
time_q3_obs <- ifelse(length(ind3)==0, NA, datcens_fit$time[ind3])

(q <- c(time_q1_obs, time_q2_obs, time_q3_obs))

as.numeric((q[1] - 2 * q[2] + q[3]) / (q[3] - q[1]))

Calc_Skew(times = datcens$times, cens = datcens$cens)




#-------------------------------------------------------------------------------
# Symmetric Administrative Censoring Example
#-------------------------------------------------------------------------------

set.seed(28657)

dat <- rweibull(243, shape = 3.5, scale = 20)
cens <- c()
quant <- quantile(dat, 0.7)
cens <- ifelse(dat < quant, 1, 0)
datnocens <- data.frame(times = dat, cens = rep(1, 243))
dat[dat > quant] <- quant
datcens <- data.frame(times = dat, cens = cens)


datcens_fit <- survfit(Surv(datcens$times, datcens$cens) ~ 1)

km_cens <- data.frame("times" = datcens_fit$time,
                      "surv" = datcens_fit$surv)

ggplot(km_cens, aes(x = times, y = surv)) +
  geom_step(linewidth = 1, colour = "blue") +
  labs(title = "Kaplan-Meier Survival Curve",
       x = "Time",
       y = "Survival probability") +
  theme_minimal() + ylim(c(0,1))


(q <- quantile(datnocens$times, c(0.25,0.5,0.75)))

as.numeric((q[1] - 2 * q[2] + q[3]) / (q[3] - q[1]))

ind <- which(round(datcens_fit$surv,2) == c(0.75,0.5,0.25))
ind1 <- which(round(datcens_fit$surv,2) == 0.75)
ind2 <- which(round(datcens_fit$surv,2) == 0.5)
ind3 <- which(round(datcens_fit$surv,2) == 0.25)
time_q1_obs <- ifelse(length(ind[1])==0, NA, datcens_fit$time[ind[1]])
time_q2_obs <- ifelse(length(ind[2])==0, NA, datcens_fit$time[ind[2]])
time_q3_obs <- ifelse(length(ind[3])==0, NA, datcens_fit$time[ind[3]])

(q <- c(time_q1_obs, time_q2_obs, time_q3_obs))

as.numeric((q[1] - 2 * q[2] + q[3]) / (q[3] - q[1]))

min_surv <- min(datcens_fit$surv)

ind <- which(round(datcens_fit$surv,2) == c(round(1 - min_surv,2), 0.5, round(min_surv,2)))
ind1 <- which(round(datcens_fit$surv,2) == round(1 - min_surv,2))
ind2 <- which(round(datcens_fit$surv,2) == 0.5)
ind3 <- which(round(datcens_fit$surv,2) == round(min_surv,2))
time_q1_obs <- ifelse(length(ind[1])==0, NA, datcens_fit$time[ind[1]])
time_q2_obs <- ifelse(length(ind[2])==0, NA, datcens_fit$time[ind[2]])
time_q3_obs <- ifelse(length(ind[3])==0, NA, datcens_fit$time[ind[3]])

(q <- c(time_q1_obs, time_q2_obs, time_q3_obs))

as.numeric((q[1] - 2 * q[2] + q[3]) / (q[3] - q[1]))

Calc_Skew(times = datcens$times, cens = datcens$cens)

















































