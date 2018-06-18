givens <- function(L,k){
  p <- dim(L)[1]
  if( k>p ) return ("Wrong input of k!")
  Lk <- L[-k,]
  mk <- k
  while( mk < p ){
    mx <- Lk[mk,mk:(mk+1)]
    lmx <- sqrt(sum(mx*mx))
    Lk[mk,mk:(mk+1)] <- c(lmx,0)
    if( mk < p-1 ) Lk[(mk+1):(p-1), mk:(mk+1)] <- Lk[(mk+1):(p-1), mk:(mk+1)] %*% givens(mx, lmx)
    mk <- mk + 1
  }
  Lk[,-p]
}
