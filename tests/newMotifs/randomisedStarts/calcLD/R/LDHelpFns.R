# Help functions for calculating LD

# Relabel genotypes according to fitness

RelabelGenotypeFrequencies <- function(d_rankings, l_parentalFreqs) {
  ranks <- d_rankings %>% select(starts_with("wpar"))
  # Get the labels
  genotype_names <- paste0("p", substr(colnames(ranks), 5, 6))
  
  # order rankings by fitness
  order_indices <- t(apply(ranks, 1, order))

  ab_id <- genotype_names[order_indices[,1]]
  AB_id <- genotype_names[order_indices[,4]]
  
  # Intermediates can be randomly assigned
  Ab_id <- genotype_names[order_indices[,2]]
  aB_id <- genotype_names[order_indices[,3]]
  
  # Create output
  out <- list(pAB = numeric(length(AB_id)),
              pAb = numeric(length(AB_id)),
              paB = numeric(length(AB_id)),
              pab = numeric(length(AB_id)))
  
  # Fill output
  # Convert l_parentalFreqs to a matrix for vectorized indexing
  n <- length(AB_id)
  parentalMat <- t(do.call(cbind, l_parentalFreqs))
  
  # Vectorized filling of the output lists
  out$pAB <- mapply(function(geno, i) parentalMat[geno, i], AB_id, seq_len(n))
  out$pAb <- mapply(function(geno, i) parentalMat[geno, i], Ab_id, seq_len(n))
  out$paB <- mapply(function(geno, i) parentalMat[geno, i], aB_id, seq_len(n))
  out$pab <- mapply(function(geno, i) parentalMat[geno, i], ab_id, seq_len(n))  
  return(out)
}

CalcLD <- function(l_freqs, metric = "D") {
  names(l_freqs$pAB) <- NULL
  names(l_freqs$pAb) <- NULL
  names(l_freqs$paB) <- NULL
  names(l_freqs$pab) <- NULL
  
  switch (metric,
    "D" = return( ( l_freqs$pAB * l_freqs$pab ) - ( l_freqs$pAb * l_freqs$paB ))
  )
  
}
