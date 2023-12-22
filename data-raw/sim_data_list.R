# Simulate data for demonstrating the use of CSNet
#
# The following codes are a simplifed version of the simulation codes
# used in the original paper.
# For more details, please see the original version at
# https://github.com/ChangSuBiostats/CSNet_analysis

n <- 600
p <- 100
log_var <- 8.0
K <- 2
rhos <- c(0.8, 0.8)
betas <- c(2,1)
NB_exper_pars <- list(Tsize=200, S=6e+07)
cor_model <- 'AR_10'
seed <- 2345

# simulate two cell-type-specific cluster
cor_p <- round(p / (K + 1))
sub_cl <- lapply(1:K, function(k) (((k-1)*cor_p+1) :(k*cor_p)))


# simulate and fix cell type proportions
set.seed(1)
pi_ct1 <- rbeta(n, betas[1], betas[2])
pi_m <- matrix(c(pi_ct1, 1-pi_ct1), ncol = 2)
colnames(pi_m) <- paste('cell type', c(1,2))

# simulate covariance
sigma_sq_1 <- rep(exp(log_var), p)
sigma_sq_2 <- rep(exp(log_var), p)
sigma_sq_2[(2*cor_p+1): p] <- rep(exp(log_var+1), p-2*cor_p) 
sigma_sq_star_list <- list(sigma_sq_1, sigma_sq_2)

# generate correlation matrices
R_star_list <- list()
NB_par_list <- Sigma_star_list <- R_star_list <- list()

# truncated AR(1) as in Table 1
cor_m <- matrix(0, nrow = cor_p, ncol = cor_p)
cors <- rhos[1]^(0:(cor_p-1))
m <- 10
for(j in 1:cor_p){
  max_ind <- min(cor_p, j+m)
  cor_m[j, j:max_ind] <- cors[1:(max_ind-j+1)]
}
cor_m[lower.tri(cor_m)] <- t(cor_m)[lower.tri(cor_m)]

for(k in 1:K){
  cor_ind <- sub_cl[[k]]
  Sigma_star <- diag(sigma_sq_star_list[[k]])
  Sigma_star[cor_ind, cor_ind] <- cor_m * sigma_sq_star_list[[k]][cor_ind]
  Sigma_star_list[[k]] <- Sigma_star
  R_star_list[[k]] <- Sigma_star / outer(sqrt(sigma_sq_star_list[[k]]),
                                         sqrt(sigma_sq_star_list[[k]]))
}

# generate cell-type-specific mean expressions

set_phi <- function(sigma_sq, Tsize=200, S=6e+07, use_approx=F){
  phi <- 1/S * sqrt(sigma_sq / Tsize + 1/4) - 1/(2*S)
  return(phi)
}

mu_star_list <- list()
for(k in 1:K){
  gene_mu_vec <- numeric(p)
  for(j in 1:p){
    phi <- set_phi(sigma_sq_star_list[[k]][j],
                   Tsize=NB_exper_pars$Tsize, S=NB_exper_pars$S)
    gene_mu <- NB_exper_pars$S * NB_exper_pars$Tsize * phi
    gene_mu_vec[j] <- gene_mu
  }
  mu_star_list[[k]] <- gene_mu_vec
}

# simulate negative binomial samples 

gen_cor_NB_copula <- function(n, cor_mat, sigma_sq, mu){
  p <- nrow(cor_mat)
  cor_mat_up <- chol(cor_mat)
  copula <- matrix(rnorm(p*n),nrow = p)
  copula <- pnorm(t(cor_mat_up) %*% copula)
  
  exp_matrix <- matrix(NA, nrow=p, ncol=n)
  rownames(exp_matrix) <- paste0("Gene", 1:p)
  
  size_vec <- mu^2 / (sigma_sq - mu)
  for (i in 1:p){
    exp_matrix[i,] <- qnbinom(copula[i,], mu = mu[i], size = size_vec[i])
  }
  return(t(exp_matrix))
}


gene_exp_list <- list()
for(k in 1:K){
  cor_ind <- sub_cl[[k]]
  # make sure the random seed for two cell types are different
  ct_specific_seed <- seed + (k^2+43)
  gene_exp_mat <- matrix(0, n, p)
  
  # generate correlated expressions via copula
  tmp_cor <- gen_cor_NB_copula(n, R_star_list[[k]][cor_ind, cor_ind],
                               sigma_sq_star_list[[k]][cor_ind], 
                               mu_star_list[[k]][cor_ind])
  gene_exp_mat[, sub_cl[[k]]] <- tmp_cor
  
  # generate independent NB samples
  tmp_indpt <- matrix(NA, n, p-cor_p)
  # extract mean and variance statistics for independent genes
  indpt_ind <- which(!(1:p %in% sub_cl[[k]]))
  gene_mu_sub <- mu_star_list[[k]][indpt_ind]
  gene_var_sub <- sigma_sq_star_list[[k]][indpt_ind]
  set.seed(ct_specific_seed)
  for(j in 1:(p-cor_p)){
    gene_mu <- gene_mu_sub[j]
    gene_var <- gene_var_sub[j]
    gene_p <- gene_mu / gene_var
    tmp_indpt[,j] <- stats::rnbinom(n, 
                                    size = gene_mu * gene_p / (1-gene_p),
                                    mu = gene_mu)
  }
  gene_exp_mat[, indpt_ind] <- tmp_indpt
  
  gene_exp_list[[k]] <- gene_exp_mat
}


# obtain bulk expressions
bulk_exp <- matrix(0, nrow = n, ncol = p)
for(k in 1:K){
  bulk_exp <- bulk_exp + gene_exp_list[[k]] * pi_m[,k]
}

sim_data_list <- list(P = pi_m,
                      ct_X = gene_exp_list)

usethis::use_data(sim_data_list)
