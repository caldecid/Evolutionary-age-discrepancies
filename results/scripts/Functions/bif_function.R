
# Function for correcting bifurcating speciation --------------------------

### Kendall 1946
p1 <- function(lam, mu, t){
  p1t <- ((lam - mu)^2 * exp(-(lam-mu)* t)) / (lam - mu * exp(-(lam - mu)* t)) ^2 
  return(p1t)
}


get_sp_age_prob <- function(lam, mu, node_age,
                            # discretization of species age
                            n_bins_per_myr = 100,
                            # discretization of the BD process (bin size)
                            dt = 0.000001){
  tot_t = seq(0,node_age, length.out=n_bins_per_myr * node_age)
  # probability of the node age with no descendants in the absence of extinction
  # i.e. prob of no speciation
  p_node_age_PB = p1(lam, 0, dt) ^ (tot_t / dt)
  # probability of no events along a lineage
  p_no_events = p1(lam, mu, dt) ^ (tot_t / dt)
  # conditional prob of no events
  p = p_no_events / p_node_age_PB
  # prob of specieme age
  p_at_node = p[length(tot_t)]
  # prob of species age younger than node age (old version)
  # p_before_node = p[-length(tot_t)] / sum((p[-length(tot_t)])) * (1 - p_at_node)
  # new version:
  p_before_node = p_no_events[-length(tot_t)] / sum((p_no_events[-length(tot_t)])) * (1 - p_at_node)
  # prob vector along branch
  p_vec = c(p_before_node, p_at_node)
  res = NULL
  res$time = tot_t
  res$prob = p_vec
  return(res)
}

