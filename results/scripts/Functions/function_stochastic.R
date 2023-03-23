##sourcing the libraries and the directories
source(file.path(getwd(), "/source.R"))


# function stochastic species ages ----------------------------------------

##function for rotating randomly 80% of the nodes of a phylogeny and repeating
## this procedure 1000 times. The aim is to generate a distribution of probable
## ages for each tip. 

stochastic.sp.age <- function(tree){
  
  ##tree = phylogenetic tree
  
  ##rotate trees function
  trees.1000 <- rotate_nodes(tree)

  ##simulating taxonomy (only budding speciation)
  taxa.1000 <- lapply(trees.1000, FossilSim::sim.taxonomy, beta = 0)
  
  ##calculating ages
  taxa.age <- Map(getAgesSpecies, taxa.1000, trees.1000)
  
  ##binding rows
  clado.df <- do.call("rbind", taxa.age)
  
###get dataframe nested by label (species)
  getNestedDataframe <- function(x) {
    
    #x = dataframe with the species'ages of the 1000 trees
    
    ##nesting the dataframe
    clado.df.nest <- x %>% select(label, Age) %>% 
      group_by(label) %>% nest()
    
    ##confidence intervalS
    CI <- vector("list", length = nrow(clado.df.nest))
    
    for(i in seq_along(clado.df.nest[[2]])){
      CI[[i]] <- hdi(clado.df.nest[[2]][[i]], 0.95)
    }
    
    for(i in seq_along(clado.df.nest[[2]])){
      clado.df.nest[[2]][[i]] <- as.data.frame(table(clado.df.nest[[2]][[i]])) 
      clado.df.nest[[2]][[i]]$Probability <- clado.df.nest[[2]][[i]]$Freq/1000
      clado.df.nest[[2]][[i]] <- clado.df.nest[[2]][[i]][,-2] ##eliminating freq
      clado.df.nest$CI_low[i] <- CI[[i]][1] ##assigning low CI
      clado.df.nest$CI_up[i] <- CI[[i]][2] ##assigning high CI
    }
    colnames(clado.df.nest) <- c("species", "dist_ages", "CI_low", "CI_up")
    return(clado.df.nest)
  }
  
  suppressWarnings(getNestedDataframe(clado.df))
  ##output a dataframe with species names and a nested dataframe with 
  #distribution of ages, and the Confidence interval of that distribution
  
}








