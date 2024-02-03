### Yang and Rannala
p0t <- function(rho, l, m, t){
  # prob of 1 or more descendents
  p = rho * (l - m) / ( rho * l + (l * (1 - rho) - m) * exp((m - l) * t) )
  return(p)
}

p1t <- function(rho, l, m, t){
  p = 1 / rho * p0t(rho, l, m, t) ^2 * exp((m - l) * t)
  return(p)
}

get_sp_age_prob_new <- function(lam, mu, node_age,
                            rho = 1,
                            # discretization of species age
                            n_bins_per_myr = 100,
                            # discretization of the BD process (bin size)
                            dt = 0.00001){
  tot_t = seq(0,node_age, length.out=n_bins_per_myr * node_age)
  # probability of the node age with no descendants in the absence of extinction
  # i.e. prob of no speciation
  p_node_age_PB = p1(lam, 0, dt) ^ (tot_t / dt)
  # probability of no events along a lineage
  p_no_events = p1(lam, mu, dt) ^ (tot_t / dt)
  # conditional prob of no events times the probability of sampling the one descendant
  p = p_no_events / p_node_age_PB * rho
  # prob of speciation age x prob sister lineage with >=1 descendendents
  p_at_node = p[length(tot_t)] * p0t(rho, lam, mu, tot_t[length(tot_t)])
  # prob of species age younger than node age x prob of no descendents on the other side
  p_before_node = p_no_events[-length(tot_t)] * (1 - p0t(rho, lam, mu, tot_t[-length(tot_t)]))
  if (max(p_before_node) > 0){
    # normalize by the prob of dpeciation at the node
    p_before_node = p_before_node / sum(p_before_node) * (1 - p_at_node)
  }
  # prob vector along branch
  p_vec = c(p_before_node, p_at_node)
  # p_vec = p_vec / sum(p_vec)
  res = NULL
  res$time = tot_t
  res$prob = p_vec
  return(res)
}