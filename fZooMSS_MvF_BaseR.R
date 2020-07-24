fZooMSS_MvF_BaseR <- function(ngrps, curr_min_size, curr_max_size, A_iter, C_iter, Nb_iter, S_iter, A, B, C, Nb, S){

  for(i in 1:ngrps){

    idx_curr <- (curr_min_size[i]+1):curr_max_size[i] ## Set size range index for current group

    for(j in idx_curr){ ## Find the abundance at the next size class with standard MvF
      Nb_iter[i,j] <- (S_iter[i,j] + A_iter[i,j]*Nb[i,j-1])/(C_iter[i,j])

      if(j >= (idx_curr[1]+1)){ ## Find abundance with MvF with diffusion
        k <- j - 1
        Nb[i,k] <- (S[i,k] + A[i,k] * Nb[i,k-1] + B[i,k] * Nb_iter[i,k+1]) / C[i,k]
      }

      # MvF without diffusion for last size class
      if(j == idx_curr[length(idx_curr)]){
        Nb[i,j] <- 0
      }
    }
  }
  return(Nb)
}