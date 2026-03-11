################################################################################
#-------------------------------------------------------------------------------
# We need a dataset without censoring, so using the parameters that fit our 
# skewness requirements, we generate uncensored data as follows:
#-------------------------------------------------------------------------------
################################################################################



set.seed(123)
gendata <- function(distr, shape, scale, meanlog, sdlog, cens, n){
  if (distr == "weibull"){
    
    dat <- rweibull(n, shape, scale)
    
    #generate the cesoring 
    if(cens == 1){ #no censoring
      
      delta <- rep(1, n)
      tobs <- dat
      
    } else if(cens == 0.4){ # 60% of cens
      
      if(shape == 3.5){
        censt <- rexp(n, 1/18.5)
        tobs <- pmin(dat, censt)
        delta <- as.numeric(dat <= censt)
      } else if(shape==1.7){
        censt <- rexp(n, 1/15.95)
        tobs <- pmin(dat, censt)
        delta <- as.numeric(dat <= censt)
      } else if(shape == 1){
        censt <- rexp(n, 1/13.2)
        tobs <- pmin(dat, censt)
        delta <- as.numeric(dat <= censt)
      }
      
    } else if(cens == 0.7){ # 30% of cens
      
      if(shape == 3.5){
        censt <- rexp(n, 1/48.9)
        tobs <- pmin(dat, censt)
        delta <- as.numeric(dat <= censt)
      } else if(shape == 1.7){
        censt <- rexp(n, 1/46.25)
        tobs <- pmin(dat, censt)
        delta <- as.numeric(dat <= censt)
      } else if(shape == 1){
        censt <- rexp(n, 1/45.9)
        tobs <- pmin(dat, censt)
        delta <- as.numeric(dat <= censt)
      }
      
    }
    
  }else if(distr == "lognormal"){
    
    dat <- rlnorm(n, meanlog, sdlog)
    
    #generate the cesoring 
    if(cens == 1){ #no censoring
      
      delta <- rep(1, n)
      tobs <- dat
      
    } else if(cens == 0.4){ # 60% of cens
      
      if(sdlog == 0.01){
        censt <- rexp(n, 1/7.99)
        tobs <- pmin(dat, censt)
        delta <- as.numeric(dat <= censt)
      } else if(sdlog == 0.4){
        censt <- rexp(n, 1/8.01)
        tobs <- pmin(dat, censt)
        delta <- as.numeric(dat <= censt)
      } else if(sdlog == 1.5){
        censt <- rexp(n, 1/7.25)
        tobs <- pmin(dat, censt)
        delta <- as.numeric(dat <= censt)
      }
      
    } else if(cens == 0.7){ # 30% of cens
      
      if(sdlog == 0.01){
        censt <- rexp(n, 1/20.5)
        tobs <- pmin(dat, censt)
        delta <- as.numeric(dat <= censt)
      } else if(sdlog == 0.4){
        censt <- rexp(n, 1/21.55)
        tobs <- pmin(dat, censt)
        delta <- as.numeric(dat <= censt)
      } else if(sdlog == 1.5){
        censt <- rexp(n, 1/32.5)
        tobs <- pmin(dat, censt)
        delta <- as.numeric(dat <= censt)
      }
      
    }
  }
  
  return(data.frame(times = tobs, cens = delta))
  
}


set.seed(28657)

sets <- list(distr = c("weibull", "lognormal"),
             shape = c(3.5, 1.7, 1),
             scale = 20,
             meanlog = 2,
             sdlog = c(0.01, 0.4, 1.5),
             cens=c(1,0.4,0.7),
             n = c(75,150,300))

leng <- sapply(sets, length)

for (a in 1:leng[1]) {
  
  distr <- sets$distr[a]
  
  if (distr == "weibull") {
    # Iterate over relevant parameters for Weibull
    for (b in 1:leng[2]) {
      for (c in 1:leng[3]) {
        for (f in 1:leng[6]) {
          for (g in 1:leng[7]) {
            data <- vector("list", 1000) #generate 1000 Data sets for each condition
            for (i in 1:1000) {
              data[[i]] <- gendata(
                distr = distr,
                shape = sets$shape[b],
                scale = sets$scale[c],
                cens = sets$cens[f],
                n = sets$n[g]
              )
            }
            comment(data) <- paste(
              "Distribution", distr,
              "Shape", sets$shape[b],
              "Scale", sets$scale[c],
              "Censoring", sets$cens[f],
              "n", sets$n[g]
            )
            file_name <- paste(
              "C:/Users/arnau.gomez/Desktop/GofCensSimulatios-Study/02_Data/Weibull_RandomCensoring/",
              "Distribution", distr,
              "Censoring", sets$cens[f],
              "n", sets$n[g],
              "Shape", sets$shape[b],
              ".RData", sep = '_'
            )
            save(data, file = file_name)
            rm(data, i)
          }
        }
      }
    }
  } else if(distr == "lognormal"){
    # Iterate over relevant parameters for Lognormal
    for (d in 1:leng[4]) {
      for (e in 1:leng[5]) {
        for (f in 1:leng[6]) {
          for (g in 1:leng[7]) {
            data <- vector("list", 1000)
            for (i in 1:1000) {
              data[[i]] <- gendata(
                distr = distr,
                meanlog = sets$meanlog[d],
                sdlog = sets$sdlog[e],
                cens = sets$cens[f],
                n = sets$n[g]
              )
            }
            comment(data) <- paste(
              "Distribution", distr,
              "Meanlog", sets$meanlog[d],
              "Sdlog", sets$sdlog[e],
              "Censoring", sets$cens[f],
              "n", sets$n[g]
            )
            file_name <- paste(
              "C:/Users/arnau.gomez/Desktop/GofCensSimulatios-Study/02_Data/Lognormal_RandomCensoring/",
              "Distribution", distr,
              "Censoring", sets$cens[f],
              "n", sets$n[g],
              "Sdlog", sets$sdlog[e],
              ".RData", sep = '_'
            )
            save(data, file = file_name)
            rm(data, i)
            
          }
        }
      }
    }
  }
}





















