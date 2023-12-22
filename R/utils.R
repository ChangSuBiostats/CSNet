# helper functions

# convert covariance matrices to correlation matrices
get_cor_from_cov <- function(cov_mat){
  sd_vec <- sqrt(diag(cov_mat))
  cor_mat <- cov_mat / outer(sd_vec, sd_vec)
  # pre-process NA due to 0 or negative variance estimates
  na_genes <- which(is.na(diag(cor_mat)))
  cor_mat[na_genes,] <- 0
  cor_mat[,na_genes] <- 0
  diag(cor_mat) <- 1
  # pre-process correlation that is out of [-1,1]
  cor_mat[cor_mat > 1] <- 0
  cor_mat[cor_mat < -1] <- 0
  return(cor_mat)
}

# -
# evaluate matrix norms
# -
F_norm <- function(Sigma, Sigma_star){
  error <- Sigma - Sigma_star
  return(sqrt(sum(error**2)))
}

op_norm <- function(Sigma, Sigma_star){
  #require(svd)
  error <- Sigma - Sigma_star
  s <- svd(error)$d[1]
  #s <- propack.svd(error,3)$d[1]
  return(s)
}

eval_TPR <- function(Sigma, Sigma_star, off_diagonal=F){
  if(off_diagonal){
    # Evaluate TPR with off-diagonal entries only
    Sigma <- Sigma[upper.tri(Sigma)]
    Sigma_star <- Sigma_star[upper.tri(Sigma_star)]
  }
  sum((Sigma*Sigma_star)!=0) / sum(Sigma_star!=0)
}

eval_FPR <- function(Sigma, Sigma_star, off_diagonal=F){
  if(off_diagonal){
    # Evaluate FPR with off-diagonal entries only
    Sigma <- Sigma[upper.tri(Sigma)]
    Sigma_star <- Sigma_star[upper.tri(Sigma_star)]
  }
  sum((Sigma!=0)*(Sigma_star==0)) / sum(Sigma_star==0)
}

# -
# thresholding operators
# -
generalized_th <- function(mat, th, type = 'hard', th_diag=T){
  diag_vals <- diag(mat)
  if(type == 'hard'){
    th_mat <- hard_thresholding(mat, th)
  }else if(type == 'soft'){
    th_mat <- soft_thresholding(mat, th)
  }else if(type == 'SCAD'){
    th_mat <- SCAD_thresholding(mat, th)
  }
  if(!th_diag) diag(th_mat) <- diag_vals
  return(th_mat)
}

hard_thresholding <- function(mat, th){
  mat[abs(mat) <= th] <- 0
  return(mat)
}

soft_thresholding <- function(mat, th){
  ind <- (abs(mat) > th)# & valid_ind
  val <- mat[ind]
  new_val <- sign(val) * (abs(val) - th)
  # set values smaller than th to 0
  mat[(abs(mat) <= th)] <- 0
  # set values greater than th to sign(val)(val-th)
  mat[ind] <- new_val
  return(mat)
}

# SCAD thresholding function is a linear interpolation between
# soft thresholding up to 2\lambda and hard thresholding after a\lambda
# value a = 3.7 was recommended by Fan and Li (2001),
SCAD_func <- function(x, th, a=3.7){
  if(abs(x) <= th){
    return(0)
  }else if(abs(x) <= 2*th){
    return(sign(x) * (abs(x) - th))
  }else if(abs(x) <= a*th){
    return(sign(x) * ((a-1)/(a-2)*(abs(x) - 2*th) + th))
  }else{
    return(x)
  }
}

SCAD_thresholding <- function(mat, th, a=3.7){
  ind_hard <- abs(mat) <= th
  mat[ind_hard] <- 0
  ind_soft <- !ind_hard & (abs(mat) <= 2*th)
  mat[ind_soft] <- sign(mat[ind_soft]) * (abs(mat[ind_soft]) - th)
  ind_inter <- (abs(mat) > 2*th) & (abs(mat) <= a*th)
  mat[ind_inter] <- sign(mat[ind_inter]) * ((abs(mat[ind_inter]) - 2*th) * (a-1) / (a-2) + th)
  return(mat)
}