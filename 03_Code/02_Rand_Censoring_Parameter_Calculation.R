################################################################################
#-------------------------------------------------------------------------------
# To generate the radom censored data we need to know the value of the parameter
# lambda that generates the exact percentage of censoring we need. 
# Now we are searching for this value by simulating:
#-------------------------------------------------------------------------------
################################################################################


#############################################################
# 60% of censoring
#############################################################


#----------------------------------
# Weibull Symmetric, shape = 3.5
#----------------------------------

sec <- seq(0.01, 0.09, 0.001)
cens_percentage <- numeric()
best_lambda <- numeric()

for(j in 1:100){
  for(i in sec){
  
    dat <- rweibull(10000, 3.5, 20)
    censt <- rexp(10000, i)
    tobs <- pmin(dat, censt)
    delta <- as.numeric(dat <= censt)
    datos <- data.frame(times = tobs, cens = delta)
    cens_percentage <- c(cens_percentage, sum(datos$cens==0)/length(datos$cens))
  
  }
  best_lambda <- c(best_lambda, mean(sec[which(round(cens_percentage, 1) == 0.6)]))
  best_lambda
  cens_percentage <- numeric()
}

(best_lambda_weibull_symmetric <- mean(best_lambda))


#----------------------------------
# Weibull Asymmetric, shape = 1.7
#----------------------------------

cens_percentage <- numeric()
best_lambda <- numeric()

for(j in 1:100){
  for(i in sec){
    
    dat <- rweibull(10000, 1.7, 20)
    censt <- rexp(10000, i)
    tobs <- pmin(dat, censt)
    delta <- as.numeric(dat <= censt)
    datos <- data.frame(times = tobs, cens = delta)
    cens_percentage <- c(cens_percentage, sum(datos$cens==0)/length(datos$cens))
    
  }
  best_lambda <- c(best_lambda, mean(sec[which(round(cens_percentage, 1) == 0.6)]))
  best_lambda
  cens_percentage <- numeric()
}

(best_lambda_weibull_asymmetric <- mean(best_lambda))


#----------------------------------
# Weibull Extremely Asymmetric, shape = 1
#----------------------------------

cens_percentage <- numeric()
best_lambda <- numeric()

for(j in 1:100){
  for(i in sec){
    
    dat <- rweibull(10000, 1, 20)
    censt <- rexp(10000, i)
    tobs <- pmin(dat, censt)
    delta <- as.numeric(dat <= censt)
    datos <- data.frame(times = tobs, cens = delta)
    cens_percentage <- c(cens_percentage, sum(datos$cens==0)/length(datos$cens))
    
  }
  best_lambda <- c(best_lambda, mean(sec[which(round(cens_percentage, 1) == 0.6)]))
  best_lambda
  cens_percentage <- numeric()
}

(best_lambda_weibull_extra_asymmetric <- mean(best_lambda))


table_60_W <- data.frame(
  Skewness = c("Symmetric", "Asymmetric", "Extremely Asymmetric"),
  Lambda = c(best_lambda_weibull_symmetric, best_lambda_weibull_asymmetric, best_lambda_weibull_extra_asymmetric),
  "Lamba frac" = c("1/18.5", "1/15.95", "1/13.2"),
  Calc = c(1/18.5, 1/15.95, 1/13.2))




#----------------------------------
# Lognormal Symmetric, sdlog = 0.01
#----------------------------------

sec <- seq(0.05, 0.2, 0.001)
cens_percentage <- numeric()
best_lambda <- numeric()

for(j in 1:100){
  for(i in sec){
    
    dat <- rlnorm(10000, sdlog = 0.01, meanlog = 2)
    censt <- rexp(10000, i)
    tobs <- pmin(dat, censt)
    delta <- as.numeric(dat <= censt)
    datos <- data.frame(times = tobs, cens = delta)
    cens_percentage <- c(cens_percentage, sum(datos$cens==0)/length(datos$cens))
    
  }
  best_lambda <- c(best_lambda, mean(sec[which(round(cens_percentage, 1) == 0.6)]))
  best_lambda
  cens_percentage <- numeric()
}

(best_lambda_lognormal_symmetric <- mean(best_lambda))


#----------------------------------
# Lognormal Asymmetric, sdlog = 0.4
#----------------------------------

cens_percentage <- numeric()
best_lambda <- numeric()

for(j in 1:100){
  for(i in sec){
    
    dat <- rlnorm(10000, sdlog = 0.4, meanlog = 2)
    censt <- rexp(10000, i)
    tobs <- pmin(dat, censt)
    delta <- as.numeric(dat <= censt)
    datos <- data.frame(times = tobs, cens = delta)
    cens_percentage <- c(cens_percentage, sum(datos$cens==0)/length(datos$cens))
    
  }
  best_lambda <- c(best_lambda, mean(sec[which(round(cens_percentage, 1) == 0.6)]))
  best_lambda
  cens_percentage <- numeric()
}

(best_lambda_lognormal_asymmetric <- mean(best_lambda))


#----------------------------------
# Lognormal Extremely Asymmetric, sdlog = 1.5
#----------------------------------

cens_percentage <- numeric()
best_lambda <- numeric()

for(j in 1:100){
  for(i in sec){
    
    dat <- rlnorm(10000, sdlog = 1.5, meanlog = 2)
    censt <- rexp(10000, i)
    tobs <- pmin(dat, censt)
    delta <- as.numeric(dat <= censt)
    datos <- data.frame(times = tobs, cens = delta)
    cens_percentage <- c(cens_percentage, sum(datos$cens==0)/length(datos$cens))
    
  }
  best_lambda <- c(best_lambda, mean(sec[which(round(cens_percentage, 1) == 0.6)]))
  best_lambda
  cens_percentage <- numeric()
}

(best_lambda_lognormal_extra_asymmetric <- mean(best_lambda))


table_60_L <- data.frame(
  Skewness = c("Symmetric", "Asymmetric", "Extremely Asymmetric"),
  Lambda = c(best_lambda_lognormal_symmetric, best_lambda_lognormal_asymmetric, best_lambda_lognormal_extra_asymmetric),
  "Lamba frac" = c("1/7.99", "1/8.01",  "1/7.25"),
  Calc = c(1/7.99, 1/8.01, 1/7.25))






#############################################################
# 30% of censoring
#############################################################


#----------------------------------
# Weibull Symmetric, shape = 3.5
#----------------------------------

sec <- seq(0.01, 0.09, 0.001)
cens_percentage <- numeric()
best_lambda <- numeric()

for(j in 1:100){
  for(i in sec){
    
    dat <- rweibull(10000, 3.5, 20)
    censt <- rexp(10000, i)
    tobs <- pmin(dat, censt)
    delta <- as.numeric(dat <= censt)
    datos <- data.frame(times = tobs, cens = delta)
    cens_percentage <- c(cens_percentage, sum(datos$cens==0)/length(datos$cens))
    
  }
  best_lambda <- c(best_lambda, mean(sec[which(round(cens_percentage, 1) == 0.3)]))
  best_lambda
  cens_percentage <- numeric()
}

(best_lambda_weibull_symmetric <- mean(best_lambda))


#----------------------------------
# Weibull Asymmetric, shape = 1.7
#----------------------------------

cens_percentage <- numeric()
best_lambda <- numeric()

for(j in 1:100){
  for(i in sec){
    
    dat <- rweibull(10000, 1.7, 20)
    censt <- rexp(10000, i)
    tobs <- pmin(dat, censt)
    delta <- as.numeric(dat <= censt)
    datos <- data.frame(times = tobs, cens = delta)
    cens_percentage <- c(cens_percentage, sum(datos$cens==0)/length(datos$cens))
    
  }
  best_lambda <- c(best_lambda, mean(sec[which(round(cens_percentage, 1) == 0.3)]))
  best_lambda
  cens_percentage <- numeric()
}

(best_lambda_weibull_asymmetric <- mean(best_lambda))


#----------------------------------
# Weibull Extremely Asymmetric, shape = 1
#----------------------------------

cens_percentage <- numeric()
best_lambda <- numeric()

for(j in 1:100){
  for(i in sec){
    
    dat <- rweibull(10000, 1, 20)
    censt <- rexp(10000, i)
    tobs <- pmin(dat, censt)
    delta <- as.numeric(dat <= censt)
    datos <- data.frame(times = tobs, cens = delta)
    cens_percentage <- c(cens_percentage, sum(datos$cens==0)/length(datos$cens))
    
  }
  best_lambda <- c(best_lambda, mean(sec[which(round(cens_percentage, 1) == 0.3)]))
  best_lambda
  cens_percentage <- numeric()
}

(best_lambda_weibull_extra_asymmetric <- mean(best_lambda))


table_30_W <- data.frame(
  Skewness = c("symmetric", "Asymmetric", "Extremely Asymmetric"),
  Lambda = c(best_lambda_weibull_symmetric, best_lambda_weibull_asymmetric,
             best_lambda_weibull_extra_asymmetric),
  "Lamba frac" = c("1/48.9", "1/46.25", "1/45.9"),
  Calc = c(1/48.9, 1/46.25, 1/45.9))





#----------------------------------
# Lognormal Symmetric, sdlog = 0.01
#----------------------------------

sec <- seq(0.01, 0.15, 0.001)
cens_percentage <- numeric()
best_lambda <- numeric()

for(j in 1:100){
  for(i in sec){
    
    dat <- rlnorm(10000, sdlog = 0.025, meanlog = 2)
    censt <- rexp(10000, i)
    tobs <- pmin(dat, censt)
    delta <- as.numeric(dat <= censt)
    datos <- data.frame(times = tobs, cens = delta)
    cens_percentage <- c(cens_percentage, sum(datos$cens==0)/length(datos$cens))
    
  }
  best_lambda <- c(best_lambda, mean(sec[which(round(cens_percentage, 1) == 0.3)]))
  best_lambda
  cens_percentage <- numeric()
}

(best_lambda_lognormal_symmetric <- mean(best_lambda))


#----------------------------------
# Lognormal Asymmetric, sdlog = 0.4
#----------------------------------

cens_percentage <- numeric()
best_lambda <- numeric()

for(j in 1:100){
  for(i in sec){
    
    dat <- rlnorm(10000, sdlog = 0.3, meanlog = 2)
    censt <- rexp(10000, i)
    tobs <- pmin(dat, censt)
    delta <- as.numeric(dat <= censt)
    datos <- data.frame(times = tobs, cens = delta)
    cens_percentage <- c(cens_percentage, sum(datos$cens==0)/length(datos$cens))
    
  }
  best_lambda <- c(best_lambda, mean(sec[which(round(cens_percentage, 1) == 0.3)]))
  best_lambda
  cens_percentage <- numeric()
}

(best_lambda_lognormal_asymmetric <- mean(best_lambda))


#----------------------------------
# Lognormal Extremely Asymmetric, sdlog = 1.5
#----------------------------------

cens_percentage <- numeric()
best_lambda <- numeric()

for(j in 1:100){
  for(i in sec){
    
    dat <- rlnorm(10000, sdlog = 1.5, meanlog = 2)
    censt <- rexp(10000, i)
    tobs <- pmin(dat, censt)
    delta <- as.numeric(dat <= censt)
    datos <- data.frame(times = tobs, cens = delta)
    cens_percentage <- c(cens_percentage, sum(datos$cens==0)/length(datos$cens))
    
  }
  best_lambda <- c(best_lambda, mean(sec[which(round(cens_percentage, 1) == 0.3)]))
  best_lambda
  cens_percentage <- numeric()
}

(best_lambda_lognormal_extra_asymmetric <- mean(best_lambda))


table_30_L <- data.frame(
  Skewness = c("Symmetric", "Asymmetric", "Extremely Asymmetric"),
  Lambda = c(best_lambda_lognormal_symmetric, best_lambda_lognormal_asymmetric, best_lambda_lognormal_extra_asymmetric),
  "Lamba frac" = c("1/20.5", "1/21.55",  "1/32.5"),
  Calc = c(1/20.5, 1/21.55, 1/32.5))






################################################################################
#-------------------------------------------------------------------------------
# Results:
#-------------------------------------------------------------------------------
################################################################################


table_60_W
table_60_L

table_30_W
table_30_L







