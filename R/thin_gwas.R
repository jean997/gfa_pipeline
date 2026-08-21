thin_gwas <- function(beta_hat, se_beta_hat, R,
                       eps = c(0.5, 0.5), R_LD = ld_mat_list){
  n <- nrow(beta_hat)
  p <- ncol(beta_hat)
  K <- length(eps) # number of folds
  stopifnot(all.equal(dim(se_beta_hat), c(n, p)), dim(R) == c(p, p),
            all(eps > 0), isTRUE(all.equal(sum(eps), 1)))

  zeros <- matrix(0, nrow = n, ncol = p)

  # gen_bhat_from_b generates random variables with the correct covariance.
  ld_size <- sum(sapply(R_LD, nrow))
  af <- rep(0.5, ld_size) # af doesn't matter here
  N <- matrix(1, nrow = p, ncol = p)
  W <- lapply(seq_len(K), function(k){
    (GWASBrewer:::gen_bhat_from_b(b_joint = zeros, N = N,
                                  trait_corr = R,
                                  R_LD = R_LD, af = af, output_geno_scale = "sd")$beta_hat)*se_beta_hat
  })

  S <- Reduce(`+`, Map(function(w, e) sqrt(e) * w, W, eps))  # sum_j sqrt(e_j) W_j

  folds <- Map(function(w, e) {
    b <- e * beta_hat + sqrt(e) * w - e * S
    list(beta_hat = b / e, se_beta_hat = se_beta_hat / sqrt(e))
  }, W, eps)

  names(folds) <- paste0("fold", seq_len(K))
  folds
}
