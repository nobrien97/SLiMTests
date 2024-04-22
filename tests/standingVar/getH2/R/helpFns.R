decompressRel <- function(relPos, relVals, n) {
  # converts compressed relatedness matrix to a regular form
  result <- matrix(double(n*n), nrow = n)
  diag(result) <- 1
  
  # Convert relPos/relVals to list form
  listLength <- n * (n -1) / 2
  relList <- vector(mode = "numeric", length = listLength)
  relList[relPos+1] <- relVals # correct for zero based index with + 1
  result[upper.tri(result)] <- relList
  result <- result + t(result) - diag(diag(result)) # make symmetric
  
  colnames(result) <- 1:n
  rownames(result) <- 1:n
  return(result)
}


decompressHap <- function(compHaplos, n, m, ploidy = 2L) {
  # Converts a compressed haplotype matrix to a sparse form
  result <- matrix(integer(n*m), nrow = n, ncol = m)
  result[compHaplos+1] <- 1
  return(result)
}


hapToGen <- function(haplos, ploidy = 2L) {
  # Converts a matrix of haplotypes to a genotype matrix

    # First make sure we have a proper ploidy for our haplotypes
  stopifnot(nrow(haplos) %% ploidy == 0)
  
  res <- matrix(rep(0, nrow(haplos)%/%ploidy))
  haplo_seq <- seq(from = 1, to = nrow(haplos), by = ploidy)
  res <- haplos[haplo_seq,]
  if (ploidy == 1) {
    return(res)
  }
  for (i in 1:(ploidy-1)) {
    res <- res + haplos[haplo_seq + i,]
  }
  row.names(res) <- NULL
  return(res)
}

colVar <- function(m) {
  # Calculates variances of columns in a matrix/dataframe
  result <- double(ncol(m))
  for (i in seq_len(ncol(m))) {
    result[i] <- var(m[,i]) 
  }
  return(result)
}

addFixedHaplos <- function(haplos, fixedPos) {
  haplos[,fixedPos+1] <- 1
  return(haplos)
}
