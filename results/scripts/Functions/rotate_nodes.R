
# Rotating nodes of the phylogeny -----------------------------------------


# rotated nodes phytools ------------------------------------------------------

##function for rotating 80% of the nodes of the phylogeny (1000 times)
##using phytools

rotate_nodes <- function(tree){
  ##list to store the rotated trees
  trees.1000 <- vector("list", length = 1000)
  
  ###rotating the nodes of the tree
  Done <- FALSE
  success_rotate <- 0 
  while (!Done) {
    for(i in 1:length(trees.1000)){
      
      trees.1000[[i]] <- tree
      ###only works properly until 80% of the nodes
      rotated_tree <- tryCatch(phytools::rotateNodes(trees.1000[[i]],
                           nodes = sample(unique(trees.1000[[i]][["edge"]][,1]),
                                            size = 0.8*Nnode(trees.1000[[i]]))),
                               error = function(e) NA, warning = function(w) NA)
      if (!all(is.na(rotated_tree))) {
        trees.1000[[i]] <- rotated_tree
        success_rotate <- success_rotate + 1
      }
      if (success_rotate == length(trees.1000)) {
        Done <- TRUE
      }
    }
  }
  return(trees.1000)
}



# rotating nodes APE ------------------------------------------------------

##function for rotating 100% of the nodes of a phylogeny (1000 times) using
##ape. Torsten's implementation

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

