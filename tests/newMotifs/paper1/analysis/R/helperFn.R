ModelFromIndex <- function(id) {
  motifs <- c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH")
  return(motifs[id])
}

ModelFromIndexWithR <- function(id) {
  motifs <- rep(c("NAR", "PAR", "FFLC1", "FFLI1", "FFBH"), times = 3)
  return(motifs[id])
}

RFromIndex <- function(id) {
  r <- rep(c(1e-10, 1e-5, 1e-1), each = 5)
  return(r[id])
}


se <- function(x, na.rm = F) {
  if (na.rm)
    x <- x[!is.na(x)]
  
  return(sd(x)/sqrt(length(x)))
}


CI <- function(x, quantile = 0.975, na.rm = F) {
  return(qnorm(quantile) * se(x, na.rm))
}

rad2deg <- function(rad) {(rad * 180) / (pi)}
deg2rad <- function(deg) {(deg * pi / 180)}
