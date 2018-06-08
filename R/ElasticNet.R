ElasticNet <- function(XTX, XTY, lam1=NULL,lam2){
  p <- nrow(XTX)
  w2 <- diag(XTX)
  alpha1 <- 1/(1+lam2)
  alpha2 <- 1-alpha1
  H <- alpha1*(XTX)+alpha2*diag(w2)
  w <- sqrt(w2)
  w1 <- 1/w
  relamb <- NULL
  reb <- NULL
  A <- rep(F,p)
  VA <- NULL
  j <- which.max(abs(XTY)*w1)
  A[j] <- !A[j]
  VA <- c(VA,j)
  L <- matrix(w[j],1)  #L是活跃集内XTX的乔布斯基分解
  lamb <- w1[j]*abs(XTY[j])
  b <- rep(0,p)
  relamb <- c(relamb, lamb)
  reb <- rbind(reb ,b)
  while(TRUE){
    CC <- w1*(XTY - drop(H %*% b))
    SCC <- sign(CC)
    SCCA <- SCC[VA]
    td <- forwardsolve(L,w[VA]*SCCA) #w[VA]表示在原路径上系数扩大方差的倍数（步长缩小方差的倍数）
    d <- backsolve(t(L),td)
    a <- w1[-VA]*drop(matrix(alpha1*XTX[-VA,VA],ncol = sum(A))%*%matrix(d,ncol=1))
    gam <- rep(0,p)
    ww <- -b[VA]/d
    gam[VA] <- ifelse( ww>0 & ww<lamb, ww, lamb)
    #mm <- max(gam[VA])+0.00001
    if(sum(A)<p) {
      gam[-VA] <-ifelse(a*lamb<=CC[-VA], (lamb-CC[-VA])/(1-a),(lamb+CC[-VA])/(1+a))##get it
      #gam[-VA] <- ifelse(gam[-VA]>0,gam[-VA],mm)
    }
    ###########################################################
    #马景义修改
    j <- which.min(gam)
    gammin <- gam[j]
    if((!is.null(lam1))&&lamb-gammin<=lam1){
      b[VA] <- b[VA] + (lamb-lam1)*d
      return(b)
    }
    b[VA] <- b[VA] + gammin*d
    lamb <- lamb - gammin
    relamb <- c(relamb, lamb)
    reb <- rbind(reb, b)

    #######################################
    if(lamb==0) break
    jj <- which(as.vector(VA)==j)
    if(length(jj)==0){
      XTXAJ <- alpha1*XTX[VA,j,drop=T]
      XTXJJ <- w2[j]
      L <- forupdate(L,XTXAJ,XTXJJ)
      VA <- c(VA,j)
    }else{
      L <- mgives(L,jj)
      VA <- VA[-jj]
    }
    A[j] <- !A[j]
  }
  list(relamb=relamb, reb=reb)
}
