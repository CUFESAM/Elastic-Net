ElasticNetCV = function(x, y){         
  np <- dim(x)                 
  n <- np[1]                            
  p <- np[2] 
  if(is.data.frame(x)){
    xn <- names(x)
  }else{
    xn <- paste("x",1:p,sep="")
  }
  
  x <- as.matrix(x, nrow=n, ncol=p)
  xm <- apply(x,2,mean)
  ym <- mean(y)
  xtx <- t(x)%*%x-n*outer(xm,xm)
  xty <- drop(t(x)%*%y)-n*xm*ym
  
  b <- NULL            
  lambdav <- NULL                          
  CV_mse <- NULL
  com <- NULL
  
  list(
    Elasticnet_ = function(lambda1= -1,lambda2){ 
      b_ <- elasticnet(xtx,xty,lam1=lambda1,lam2=lambda2) 
      if( lambda1 > 0 ) {
        b0 <- ym-sum(b_$b*xm)
        b <- c(b0,b_$b)
        names(b) <- c("inter",xn)
        return(b)
      }else{
        reb <- b_$reb
        relamb <- b_$relamb
        b0 <- ym - drop(reb%*%xm)
        reb <- cbind(relamb, b0, reb)
        rownames(reb) <- NULL
        colnames(reb) <- c("lambda1","inter",xn)
        return(reb)
      }
    },
    
    cv.choosemodel = function(lambda2=NULL, k){
      index <- sample(c(rep(1:k, n%/%k), sample(1:k, n%%k)),n) 
      if( is.null(lambda2) ) lambda2 <- seq(0,1,0.1)
      lam2_seq <- lambda2 
      nlam <- length(lam2_seq)
      mse_sum <- list()
      lamsearch_list <-list()
      reb_list <- list()
      lamk_list <-c()
      for(j in 1:nlam){
        res <- elasticnet(xtx,xty,lam2 = lam2_seq[j])  
        lamsearch_list[[j]] <- res$relamb
        lamk_list[j] <- length(res$relamb)
        reb_list[[j]] <- res$reb
        mse_sum[[j]] <- rep(0,lamk_list[j])
      }
      for(i in 1:k){
        ni <- sum(index==i)
        nii <- n-ni 
        xi <- x[index==i,,drop=FALSE]
        xmi <- apply( xi , 2, mean)
        yi <- y[index==i]
        ymi <- mean(yi) 
        xmi_ <- (n*xm-ni*xmi)/nii
        ymi_ <- (n*ym-ni*ymi)/nii 
        xixi <- sweep(xi,2,xmi,"-")
        yiyi <- yi - ymi
        xitxi <- xtx - t(xixi)%*%xixi - ni*n/nii*outer(xmi-xm,xmi-xm)
        xityi <- xty - drop(t(xixi)%*%yiyi) - (ni*n/nii)*(ymi-ym)*(xmi-xm)
        for(j in 1:nlam){
          lam_b <- elasticnet(xitxi,xityi,lam2 = lam2_seq[j])
          reb <- lam_b$reb 
          relamb <- lam_b$relamb  
          reb0 <- ymi_ - drop(reb%*%xmi_)
          lam_search <- lamsearch_list[[j]]
          lamk <- lamk_list[j]
          if(relamb[1] < lam_search[1]){
            relamb <- c(lam_search[1],relamb)
            reb <- rbind(0,reb)
          }
          
          lamik <- length(relamb)
          xout <- apply(outer(relamb,lam_search[-lamk],">="), 2, sum)
          dd <- (relamb[xout]-lam_search[-lamk])/(relamb[xout]-relamb[xout+1])
          bapp <- rbind((reb[xout+1,]-reb[xout,])*dd+reb[xout,],reb[lamik,])
          tbapp <- t(bapp)
          b0lam <- ymi_ - apply(tbapp*xmi_,2,sum) 
          mse_sum[[j]] <- mse_sum[[j]]+apply((sweep(xi%*%tbapp,2,b0lam,"+")-yi)^2,2,sum)
        }
      }
      mmse=lapply(mse_sum,'/',n)
      
      mmsea <- unlist(lapply(mmse,which.min))
      mmsev <- unlist(lapply(mmse,min))
      lambda1 <- numeric(nlam)
      bb <- matrix(nrow=nlam, ncol=p)
      
      for (i in 1:nlam){
        lambda1[i] <- lamsearch_list[[i]][mmsea[i]]
        bb[i,] <- reb_list[[i]][mmsea[i],]
      }
      bb <- cbind(ym-drop(bb%*%xm),bb)
      comb <- cbind(lam2_seq,lambda1,mmsev,bb)
      rownames(comb) <- NULL
      colnames(comb) <- c("lambda2","lambda1","CVMSE","inter",xn)
      com <<- comb
      
      minmsea <- which.min(mmsev)
      minmsev <- mmsev[minmsea]
      minlambda1 <- lambda1[minmsea]
      minlambda2 <- lambda2[minmsea]
      
      minb <- bb[minmsea,]
      b <<- minb
      CV_mse <<- minmsev
      lambdav <<- c(minlambda1,minlambda2)
    },
    output = function(){
      cat('mse:',CV_mse,"\n")
      cat('b:',b,"\n")
      cat("lambda1",lambdav[1],"\n")
      cat("lambda2",lambdav[2],"\n") 
    },
    predict = function(data){
      rows <- dim(data)[1]
      data_new <- as.matrix(cbind(rep(1,rows),data))
      out <- data_new %*% b
      return(out)
    }
  )
}
