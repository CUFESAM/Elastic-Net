ElasticNetCV = function(x, y){         
  np <- dim(x)                 
  n <- np[1]                                #lassolarsCV的私有属性 
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
  #马景义修改
  #k = 10
  
  ###########
  #lassolarsCV的私有属性，功能丰富时有待继续增加
  b <- NULL            
  #mse <- NULL            
  lambdav <- NULL                          
  CV_mse <- NULL
  com <- NULL
  ############
  list(
    #lassolarsCV的方法（可从外部调佣）
    Elasticnet = function(lambda1= -1,lambda2){ 
      #####在给定lambda的场合，直接用lambda计算估计系数，不进行交叉验证
      b_ <- elasticnetcpp(xtx,xty,lam1=lambda1,lam2=lambda2) 
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
      if( is.null(lambda2) ) lambda2 <- seq(0,0.01,0.001)
      lam2_seq <- lambda2 
      nlam <- length(lam2_seq)
      mse_sum <- list()
      lamsearch_list <-list()
      reb_list <- list()
      lamk_list <-c()
      for(j in 1:nlam){
        res <- elasticnetcpp(xtx,xty,lam2 = lam2_seq[j])  
        lamsearch_list[[j]] <- res$relamb#用于网格搜索的lamb序列
        lamk_list[j] <- length(res$relamb)
        reb_list[[j]] <- res$reb
        mse_sum[[j]] <- rep(0,lamk_list[j])
      }
      for(i in 1:k){#k为交叉验证折数
        #i = 1
        ni <- sum(index==i)
        nii <- n-ni 
        xi <- x[index==i,,drop=FALSE]
        xmi <- apply( xi , 2, mean)
        yi <- y[index==i]#第i组数据后的均值
        ymi <- mean(yi) 
        xmi_ <- (n*xm-ni*xmi)/nii
        ymi_ <- (n*ym-ni*ymi)/nii 
        xixi <- sweep(xi,2,xmi,"-")
        yiyi <- yi - ymi
        xitxi <- xtx - t(xixi)%*%xixi - ni*n/nii*outer(xmi-xm,xmi-xm)#利用离差阵分解可推导此公式
        xityi <- xty - drop(t(xixi)%*%yiyi) - (ni*n/nii)*(ymi-ym)*(xmi-xm)
        for(j in 1:nlam){
          #j = 5
          lam_b <- elasticnetcpp(xitxi,xityi,lam2 = lam2_seq[j])#将去掉第i组数据的xtx,xty传入lars得到迭代轨迹数据list(relamb,reb) 
          ##lam_search <- lam_b$relamb#用于网格搜索的lambda
          ##r_b <- lam_b$reb 
          reb <- lam_b$reb #基准的参数路径
          relamb <- lam_b$relamb  #基准lamb序列
          reb0 <- ymi_ - drop(reb%*%xmi_)
          ###################
          lam_search <- lamsearch_list[[j]]#用于网格搜索的lamb序列
          lamk <- lamk_list[j]
          ###################
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
      #giveup <- length(lam_search)-length(sum_mse)
      mmse=lapply(mse_sum,'/',n)
      #mmse
      
      mmsea <- unlist(lapply(mmse,which.min))
      mmsev <- unlist(lapply(mmse,min))
      lambda1 <- numeric(nlam)
      bb <- matrix(nrow=nlam, ncol=p)
      
      for (i in 1:nlam){
        lambda1[i] <- lamsearch_list[[i]][mmsea[i]]
        bb[i,] <- reb_list[[i]][mmsea[i],]
      }
      #lambdaa=as.matrix(lambdaa,ncol=1)
      #msee=as.matrix(msee,ncol=1)
      #lam2_seq=as.matrix(lam2_seq,ncol=1)
      bb <- cbind(ym-drop(bb%*%xm),bb)
      #bb=cbind(bb1,bb)
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
      #return(list(comb=comb, lambdav=lambdav,b=b,mse=CV_mse))
      
    },
    output = function(){####输出
      cat('mse:',CV_mse,"\n")
      cat('b:',b,"\n")
      cat("lambda1",lambdav[1],"\n")
      cat("lambda2",lambdav[2],"\n")
      #com
      #if(plot_cv)plot(lam,CV_mse)
      #(cbind(lam, CV_mse))   
    }
  )
  #lassolarsCV的方法（可从外部调用）
}
