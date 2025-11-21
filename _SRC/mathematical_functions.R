log_to_normal<-function(log_mean, log_sd){
  
  ### function to transform a mean and standard deviation on the log scale back to its original mean and standard deviation on the normal scale.
  
  sigma<-sqrt(log_sd)
  mu<-exp(log_mean + 0.5 * sigma^2)
  stdev<-sqrt( (exp(sigma^2) - 1) * exp(2 * log_mean + sigma^2) )
  return(list(normal_mean=mu, normal_sd=stdev))
}

range01 <- function(x){
  ### scale range of x so that it falls into 0 and 1
  idx<-which(!is.na(x))
  return ((x-min(x[idx]))/(max(x[idx])-min(x[idx])))
}

ma <- function(x, n = 5){
  # moving average along x
  stats::filter(x, rep(1 / n, n), sides = 2)
}

Dunn = function(c1){
  n = nrow(c1$membership)
  k = length(unique(c1$cluster))
  nu = sum(c1$membership^2)/n
  c(nu,(k*nu - 1)/(k-1))
}

col_ratio<-function(x){
  return(x/sum(x,na.rm=T)*100)
}

yamartino_avg <- function(obs, drop.na=F){
  if (drop.na){
    obs<-obs[!is.na(obs)]
  }
  if(length(obs) == 0 ){
    return(NA_integer_)
  }
  
  # to radians
  obs <- obs * pi / 180
  s_a <- mean(sin(obs), na.rm = TRUE)
  c_a <- mean(cos(obs), na.rm = TRUE)
  
  avg_dir <- atan2(s_a, c_a)
  
  
  eps <- sqrt(1-(s_a**2 + c_a **2))
  std <- asin(eps) * (1+(2/sqrt(3) - 1) * eps**3)
  
  
  # to degrees
  avg_dir <- avg_dir * 180 / pi
  std <- std * 180/pi
  if(avg_dir < 0){
    avg_dir <- avg_dir + 360
  }
  
  return(avg_dir)
}

yamartino_std <- function(obs, drop.na=F){
  if (drop.na){
    obs<-obs[!is.na(obs)]
  }
  if(length(obs) == 0 ){
    return(NA_integer_)
  }
  
  # to radians
  obs <- obs * pi / 180
  s_a <- mean(sin(obs), na.rm = TRUE)
  c_a <- mean(cos(obs), na.rm = TRUE)
  
  eps <- sqrt(1-(s_a**2 + c_a **2))
  std <- asin(eps) * (1+(2/sqrt(3) - 1) * eps**3)
  
  
  # to degrees
  std <- std * 180/pi
  
  
  return(std)
}