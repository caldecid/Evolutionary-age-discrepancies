
# Rotating nodes of the phylogeny -----------------------------------------

# rotating nodes APE ------------------------------------------------------

##function for rotating 100% of the nodes of a phylogeny (1000 times) using
##ape library

rotateNodes <- function(Tree) {
  TreeRotated <- Tree
  NumTips <- Ntip(TreeRotated)
  Nodes <- (NumTips + 1):(2 * NumTips - 1)
  ShuffledNodes <- sample(Nodes, size = length(Nodes))
  for(j in ShuffledNodes) {
    TreeRotated <- ape::rotate(TreeRotated, node = j)
    TreeString <- write.tree(TreeRotated)
    TreeRotated <- read.tree(text = TreeString)
  }
  return(TreeRotated)
}

