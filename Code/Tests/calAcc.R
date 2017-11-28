library(clue)

calAcc <- function (label, result) {
  result = bestMap(label, result)
  return(sum(label == result)/length(label))
}

bestMap = function(L1, L2) {
  # bestmap: permute labels of L2 match L1 as good as possible
  # [newL2] = bestMap(L1,L2);
   
  #===========    
  L1 = as.vector(L1)
  L2 = as.vector(L2)
  if (length(L1) != length(L2)) {
    stop('size(L1) must == size(L2)')
  }
  
  Label1 = unique(L1)
  nClass1 = length(Label1)
  Label2 = unique(L2)
  nClass2 = length(Label2)
  
  nClass = max(nClass1,nClass2)
  G = matrix(0, nClass, nClass)
  for (i in 1:nClass1) {
    for (j in 1:nClass2) {
      G[i,j] = sum(L1 == Label1[i] & L2 == Label2[j])
    }
  }
  # hungResult = hungarian(-G)
  # c = hungResult[[1]]
  # t = hungResult[[2]]
  c = solve_LSAP(t(G), maximum=TRUE)
  newL2 = rep.int(0, length(L2))
  for (i in 1:nClass2) {
    newL2[L2 == Label2[i]] = Label1[c[i]]
  }
  
  
  return(newL2)
  
  #=======backup old===========
  
  L1 = L1 - min(L1) + 1      #   min (L1) <- 1;
  L2 = L2 - min(L2) + 1      #   min (L2) <- 1;
  #===========    make bipartition graph  ============
  nClass = max(max(L1), max(L2))
  G = matrix(0, nrow(nClass), ncol(nClass))
  for (i in 1:nClass) {
    for (j in 1:nClass) {
      G[i,j] = sum(L1 == i & L2 == j)
    }
  }
  #===========    assign with hungarian method    ======
  hungResult = hungarian(-G)
  c = hungResult[[1]]
  t = hungResult[[2]]
  newL2 = matrix(0, nClass, 1)
  for (i in 1:nClass) {
    newL2[L2 == i] = c[i]
  }
  
  return(newL2)
}