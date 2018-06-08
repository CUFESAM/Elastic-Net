msig <- function(p,rho){
  if(p==1){
    return ( matrix(1,p,p) )
  }
  else{
    return ( rbind( cbind( msig(p-1,rho), rho^((p-1):1)), rho^((p-1):0) ) )
  }
}
