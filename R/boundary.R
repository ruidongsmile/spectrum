boundary_opr <- function(sp_cmplx, order){
  splx <- 1:(order + 1)
  simpx_2 <- sp_cmplx[sapply(sp_cmplx, length) == order + 1]

  if(length(simpx_2) != 0){
    if (order > 0){
      simpx_1 <- sp_cmplx[sapply(sp_cmplx, length) == order]

      # boundary <- matrix(0, ncol = length(simpx_2), nrow = length(simpx_1))

      mtx_list <- simpx_2 |>
        lapply(
          FUN = function(s){
            splx |>
              lapply(
                FUN = function(idx){
                  s[splx != idx]}) |>
              lapply(FUN = function(s){
                1* sapply(simpx_1, FUN = function(x){
                  all(x == s)
                })}) |>
              unlist() |>
              matrix(ncol = length(simpx_1), byrow = TRUE)
          })

      bd_list <- mtx_list |>
        lapply(
          FUN = function(x){
            diag(sapply(splx, FUN = function(idx){(-1)^(idx+1)}), order + 1) %*% x |>
              colSums()
          }
        )

      bd_mtx <- matrix(unlist(bd_list), nrow = length(simpx_2), byrow = T) |>
        t()
    } else if(order == 0){
      bd_mtx <- matrix(rep(0, length(simpx_2)), nrow = 1)
    }
  } else{
    simpx_1 <- sp_cmplx[sapply(sp_cmplx, length) == order]
    bd_mtx <- matrix(rep(0, length(simpx_1)), ncol = 1)
  }


  return(bd_mtx)
}


pers_lap <- function(cmplx_K, cmplx_L, order = 1, W_q_minus_1_K=NULL, W_q_K=NULL, W_q_plus_1_L=NULL){
  # order is the order of simplices
  bd_opr_q_K <- cmplx_K |>
    boundary_opr(order = order)

  bd_opr_q_plus_1_L <- cmplx_L |>
    boundary_opr(order = order + 1)

  n_q_L <- length(cmplx_L[sapply(cmplx_L, length) == order + 1])
  n_q_K <- length(cmplx_K[sapply(cmplx_K, length) == order + 1])

  n_q_minus_1_K <- length(cmplx_K[sapply(cmplx_K, length) == order])
  n_q_plus_1_L <- length(cmplx_L[sapply(cmplx_L, length) == order + 2])

  if(is.null(W_q_minus_1_K)){
    W_q_minus_1_K <- if(n_q_minus_1_K != 0){rep(1, n_q_minus_1_K)} else{1}
  }
  W_q_minus_1_K_mat <- diag(W_q_minus_1_K)

  if(is.null(W_q_K)){
    W_q_K <- if(n_q_K != 0){rep(1, n_q_K)} else{1}
  }
  W_q_K_mat <- diag(W_q_K)


  if(is.null(W_q_plus_1_L)){
    W_q_plus_1_L <- if(n_q_plus_1_L !=0){rep(1, n_q_plus_1_L)}else{1}
  }
  W_q_plus_1_L_mat <- diag(W_q_plus_1_L)


  if(order == 0){
    lap_q_down_K <- 0
  }
  # else if(all(bd_opr_q_plus_1_L == 0)){
  #  lap <- W_q_K_mat %*% t(bd_opr_q_K) %*% diag(1 / W_q_minus_1_K) %*% bd_opr_q_K
  #}
  else{
    lap_q_down_K <- W_q_K_mat %*% t(bd_opr_q_K) %*% diag(1 / W_q_minus_1_K) %*% bd_opr_q_K
  }

  if(n_q_L == n_q_K){
    lap_q_up_L <- bd_opr_q_plus_1_L %*% W_q_plus_1_L_mat %*% t(bd_opr_q_plus_1_L) %*% diag(1 / W_q_K)

    lap <- lap_q_down_K + lap_q_up_L
    # return(lap_q_down_K + lap_q_up_L)
  } else{
    D_q_plus_1_L <- as.matrix(bd_opr_q_plus_1_L[seq(from=n_q_K + 1, to =n_q_L), ])

    gauss_red <- matlib::gaussianElimination(t(D_q_plus_1_L), diag(rep(1, ncol(D_q_plus_1_L))))

    R_q_plus_1_L <- t(gauss_red[, 1:nrow(D_q_plus_1_L), drop=FALSE])
    Y_mat <- t(gauss_red[, (nrow(D_q_plus_1_L) + 1):ncol(gauss_red), drop=FALSE])

    zero_idx <- which(colSums(abs(R_q_plus_1_L) > 1e-10) == 0)

    if (identical(zero_idx, integer(0))){
      lap <- lap_q_down_K
    } else{
      Z_mat <- Y_mat[, zero_idx, drop=FALSE]
      bd_opr_q_plus_1_L_K <- (bd_opr_q_plus_1_L %*% Y_mat)[1:n_q_K, zero_idx, drop=FALSE]

      lap_q_up_K_L <- bd_opr_q_plus_1_L_K %*% solve(t(Z_mat) %*% diag(1 / W_q_plus_1_L) %*% Z_mat) %*% t(bd_opr_q_plus_1_L_K) %*% diag(1 / W_q_K)
      lap <- lap_q_up_K_L + lap_q_down_K
    }
  }
  return(lap)
}








