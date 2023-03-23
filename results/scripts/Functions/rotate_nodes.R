
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
  ##list to store the rotated trees
  trees.1000 <- vector("list", length = 1000)
  
  NumTips <- Ntip(Tree)
  Nodes <- (NumTips+1):(2 * NumTips - 1)
  Done <- FALSE
  while (!Done) {
    x <- sample(Nodes, size = round(1.0 * length(Nodes)))
    TreeRotated <- Tree
    success_rotate <- 0
    for(j in 1:length(x)){
      TreeTmp <- tryCatch(ape::rotate(TreeRotated, node = x[j]),
                          error = function(e) NA,
                          warning = function(w) NA)
      if (!all(is.na(TreeTmp))) {
        # Why does rotate() messing up the tree structure?
        # Saver to convert to string and back to tree
        TreeString <- write.tree(TreeTmp)
        TreeTmp <- read.tree(text = TreeString)
        TreeRotated <- TreeTmp
        success_rotate <- success_rotate + 1
      }
    }
    if (success_rotate == length(x)) {
      png_null_device(4, 4)
      plot_ok <- tryCatch(plot.phylo(TreeRotated, plot = FALSE),
                          error = function(e) NA,
                          warning = function(w) NA)
      dev.off()
      if (all(!is.na(plot_ok))) {
        # All attempts to rotate where successfull and phylogeny can be plotted
        Done <- TRUE
      }
    }
  }
  return(TreeRotated)
}

