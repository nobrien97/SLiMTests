# Help functions for calculating LD

# Relabel genotypes according to fitness

RelabelGenotypeFrequencies <- function(d_rankings, l_parentalFreqs) {
  # Get the labels
  genotype_names <- paste0("p", substr(colnames(d_rankings)[6:9], 5, 6))
  ab_id <- genotype_names[apply(d_rankings[, 6:9], 1, which.min)]
  AB_id <- genotype_names[apply(d_rankings[, 6:9], 1, which.max)]
  
  # Intermediates can be randomly assigned
  Ab_id <- genotype_names[apply(d_rankings[, 6:9], 1, function(x) {(1:4)[-c(which.max(x), which.min(x))][1]})]
  aB_id <- genotype_names[apply(d_rankings[, 6:9], 1, function(x) {(1:4)[-c(which.max(x), which.min(x))][2]})]
  
  # Create output
  out <- list(pAB = numeric(length(AB_id)),
              pAb = numeric(length(AB_id)),
              paB = numeric(length(AB_id)),
              pab = numeric(length(AB_id)))
  
  # Fill output
  for (i in seq_along(AB_id[1:10000])) {
    out[["pAB"]][[i]] <- l_parentalFreqs[[AB_id[i]]][[i]]
    out[["pAb"]][[i]] <- l_parentalFreqs[[Ab_id[i]]][[i]]
    out[["paB"]][[i]] <- l_parentalFreqs[[aB_id[i]]][[i]]
    out[["pab"]][[i]] <- l_parentalFreqs[[ab_id[i]]][[i]]
  }
  
  return(out)
}
