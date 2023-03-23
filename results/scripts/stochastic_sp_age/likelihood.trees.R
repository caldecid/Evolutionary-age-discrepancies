##function
bd_likelihood <- function(s, e, l, m){
  # s, e: vectors with speciation and extinction times
  # l, m: speciation and extinction rates
  lik_events <- log(l) * length(s) + log(m) * length(e[e > 0])
  lik_waiting_time <-  - (l + m) * sum(s - e) 
  return(lik_events + lik_waiting_time)
}




Tree.lam.01 <- sim.bd.taxa(n = 20, lambda = 0.1, mu = 0.05,
                              complete = TRUE, numbsim = 1)[[1]]




###simulating species taxonomy; only budding speciation

Taxa.lam.01 <- sim.taxonomy(Tree.lam.01, beta = 0, lambda.a = 0)

Ages <- aggregate(Taxa.lam.01$start, list(Taxa.lam.01$sp), FUN = 'max')
Ages$te <- aggregate(Taxa.lam.01$end, list(Taxa.lam.01$sp), FUN = 'min')[,2]
colnames(Ages)[1:2] <- c('species', 'ts')
Ages

x.1 <- bd_likelihood(s = Ages$ts, e = Ages$te, l = 0.1, m = 0.05)


##################shuffling 

Tree.2 <- ape::ladderize(Tree.lam.01, right = TRUE)

Taxa.2 <- sim.taxonomy(Tree.2, beta = 0, lambda.a = 0)

Ages.2 <- aggregate(Taxa.2$start, list(Taxa.2$sp), FUN = 'max')
Ages.2$te <- aggregate(Taxa.2$end, list(Taxa.2$sp), FUN = 'min')[,2]
colnames(Ages.2)[1:2] <- c('species', 'ts')
Ages.2

x.2<- bd_likelihood(s = Ages.2$ts, e = Ages.2$te, l = 0.1, m = 0.05)

#################randomizing Nodes#####################



