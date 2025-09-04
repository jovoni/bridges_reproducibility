
get_medalt_tree = function(Xs, tree_func) {
  chrs = names(Xs)
  alleles = names(Xs[[1]])
  
  Ds = lapply(chrs, function(chr) {
    Ds = lapply(alleles, function(a) {
      input_X = Xs[[chr]][[a]]
      N = nrow(input_X)
      D = matrix(0, ncol = N, nrow = N)
      for (i in 1:(N)) {
        for (j in i:N) {
          #new_d = d_root(input_X[i,], input_X[j,])
          new_d = medalt_dist(list(input_X[i,]), list(input_X[j,]))
          D[j,i] = D[i,j] = new_d
        }
      }
      rownames(D) = colnames(D) = rownames(input_X)    
      D
    })
    
    Reduce("+", Ds)
  })
  D = Reduce("+", Ds)
  
  list(D=D, tree=tree_func(D))
}

get_root_tree = function(Xs, tree_func) {
  chrs = names(Xs)
  alleles = names(Xs[[1]])
  
  Ds = lapply(chrs, function(chr) {
    Ds = lapply(alleles, function(a) {
      input_X = Xs[[chr]][[a]]
      N = nrow(input_X)
      D = matrix(0, ncol = N, nrow = N)
      for (i in 1:(N)) {
        for (j in i:N) {
          #new_d = d_root(input_X[i,], input_X[j,])
          new_d = d_root(input_X[i,], input_X[j,])
          D[j,i] = D[i,j] = new_d
        }
      }
      rownames(D) = colnames(D) = rownames(input_X)    
      D
    })
    
    Reduce("+", Ds)
  })
  D = Reduce("+", Ds)
  
  list(D=D, tree=tree_func(D))
}

medalt_dist <- function(node1, node2) {
  d <- 0
  for (i in seq_along(node1)) {
    d <- d + disthelper(node1[[i]], node2[[i]])
  }
  return(d)
}

disthelper <- function(node1, node2) {
  if (0 %in% node1 || 0 %in% node2) {
    return(zerodisthelper(node1, node2))
  }
  return(distcalc(node1, node2))
}

distcalc <- function(node1, node2) {
  stopifnot(length(node1) == length(node2))
  if (length(node1) == 1) {
    return(abs(node1[1] - node2[1]))
  } else {
    d <- 0
    newlist <- node1 - node2
    while (length(newlist) > 0) {
      if (newlist[1] == 0) {
        newlist <- newlist[-1]
      } else if (newlist[1] > 0) {
        k <- 0
        for (i in seq_along(newlist)) {
          if (newlist[i] > 0) {
            k <- i
          } else {
            break
          }
        }
        for (i in 1:k) {
          newlist[i] <- newlist[i] - 1
        }
        d <- d + 1
      } else if (newlist[1] < 0) {
        k <- 0
        for (i in seq_along(newlist)) {
          if (newlist[i] < 0) {
            k <- i
          } else {
            break
          }
        }
        for (i in 1:k) {
          newlist[i] <- newlist[i] + 1
        }
        d <- d + 1
      }
    }
    return(abs(d))
  }
}

zerodisthelper <- function(node1, node2) {
  n1 <- rev(node1)
  n2 <- rev(node2)
  temp1 <- c()
  temp2 <- c()
  
  for (i in seq_along(n1)) {
    x1 <- n1[i]
    x2 <- n2[i]
    if (x1 == 0) {
      if (x2 == 0) {
        temp1 <- c(temp1, x1)
        temp2 <- c(temp2, x2)
      } else {
        return(1000000)
      }
    } else {
      temp1 <- c(temp1, x1)
      temp2 <- c(temp2, x2)
    }
  }
  
  return(distcalc(temp1, temp2))
}


d_root = function(x,y) {
  sum(sqrt(abs(x-y)))
}
