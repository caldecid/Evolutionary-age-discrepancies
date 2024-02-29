p0t <- function(lambda, mu, rho, t) {
  1 - (rho * (lambda - mu) / (rho * lambda + (lambda * (1 - rho) - mu) * exp((mu - lambda) * t)))
}

# the definite integral, Lambda(t1, t2) 
integratedRate <- function(lambda, mu, rho, t1, t2) {
  2 * (mu * t2 - log(lambda * rho - (lambda * (rho - 1) + mu) * exp((mu - lambda) * t2))) -
    2 * (mu * t1 - log(lambda * rho - (lambda * (rho - 1) + mu) * exp((mu - lambda) * t1))) 
}


ageDensityFunction <- function(lambda, mu, rho, t) {
  # get the integrated rate
  integrated_rate <- integratedRate(lambda, mu, rho, 0, t)
  # compute the density of the age
  age_density <- 2 * lambda * p0t(lambda, mu, rho, t) * dpois(0, integrated_rate)
  return(age_density) 
}



ageFunction <- function(lambda, mu, rho, t, v) {
  # get the integrated rate
  integrated_rate <- integratedRate(lambda, mu, rho, 0, t)
  # compute the density of the age 
  if ( t == v ) {
    age_density <- dpois(0, integrated_rate) 
  } else {
    age_density <- 2 * lambda * p0t(lambda, mu, rho, t) * dpois(0, integrated_rate) }
  return(age_density)
}


meanAgeFunction <- function(lambda, mu, rho, v) {
  # get the integrated rate
  integrated_rate <- integratedRate(lambda, mu, rho, 0, v)
  # compute the probability of zero events 
  p_no_events <- dpois(0, integrated_rate)
  # compute the mean, given that the age is not v 
  mean_with_events <- integrate(function(t) {
    t * ageDensityFunction(lambda, mu, rho, t) }, lower = 0, upper = v)$value
  # combine the probabilities
  mean_age <- p_no_events * v + mean_with_events
  return(mean_age) 
}


medianAgeFunction <- function(lambda, mu, rho, v) {
  # get the integrated rate
  integrated_rate <- integratedRate(lambda, mu, rho, 0, v)
  # compute the probability of zero events 
  p_no_events <- dpois(0, integrated_rate)
  # compute the median by optimization 
  median_age <- optim(par = v / 2, fn = function(t) {
    # compute the probability older than t 
    p <- integrate( function(s) {
      ageDensityFunction(lambda, mu, rho, s) 
    }, lower = t, upper = v)$value + p_no_events
    # compute the distance from the median 
    return(abs(p - 0.5))
  }, lower = 0, upper = v, method = "Brent")$par
  return(median_age) 
}


# SIM
forwardSimulate <- function(lambda, mu, rho, time) {
  # repeat until we get one extant lineage 
  repeat {
    # initialize simulation
    species <- data.frame(start = time, end = NA, status = "alive")
    # simulate forward 
    current_time <- time 
    repeat {
      # get the living species
      current_species <- which(species$status == "alive")
      # terminate if there are no living species
      if (length(current_species) == 0) { 
        break
      }
      # compute the rate of events
      event_rate <- (lambda + mu) * length(current_species)
      # generate a waiting time
      waiting_time <- rexp(1, event_rate) 
      current_time <- current_time - waiting_time
      # terminate if we go beyond the present 
      if (current_time < 0) {
        break 
      }
      # otherwise, perform an event
      if (runif(1) < lambda / (lambda + mu)) {
        # speciation event
        # choose affected species
        this_species <- current_species[sample.int(length(current_species), size = 1)]
        # terminate the species 
        species$status[this_species] <- "dead" 
        species$end[this_species] <- current_time
        # add new species
        new_species <- data.frame(start = c(current_time, current_time), end = NA, status = "alive")
        species <- rbind(species, new_species) 
      } else {
        # extinction event
        # choose affected species
        this_species <- current_species[sample.int(length(current_species), size = 1)]
        # terminate the species 
        species$status[this_species] <- "dead" 
        species$end[this_species] <- current_time
      }
    }
    # set missing times to 0 species$end[is.na(species$end)] <- 0
    # simulate incomplete sampling
    current_species <- which(species$status == "alive")
    species$status[current_species] <- ifelse(rbinom(length(current_species), size = 1, prob = rho) == 1, "alive", "dead")
    # check that there is exactly one extant species 
    if ( sum(species$status == "alive") == 1 ) {
      break 
    }
  }
  # done with simulation 
  return(species$start[species$status == "alive"])
}