library(foreach)
library(doParallel)
library(glmnet)
library(splines)
library(tidyverse)
library(ggplot2)
# library(WGCNA)

truef = function(x, b0, ...){
  if(ncol(x)!= length(b0)) return("invalid dim")
  x%*%b0
}

if(FALSE){
  
  x = runif(n,-1,1)
  truey = .5 + x - x^2 + 1.5*x^3 + x^4
  
  x = runif(n,-1,1)
  truey = .5 + x - 2*x^2 + 1.5*x^3 - 2*x^4
  
}
genDatx = function(b0, n, p, 
                  sigma = 1/3, int = 1){
  
  x = runif(n,-1,1)
  x1 = matrix(runif(n*p,-1,1), nrow = n)
  xmat = cbind(x, x^2, x^3, x^4, x1)
  
  truey = 0.5 + xmat%*%b0
  y = truey + rnorm(n)* sigma
  
  dat = list(truey = truey,
             y = y,
             x = xmat)
}

genDatx = function(b0, n, p, 
                  sigma = 1/3, int = 1){
  
  x = runif(n,-1,1)
  x1 = matrix(runif(n*p,0,1), nrow = n)
  xmat = cbind(x, x^2, x^3, cos(x*10), x*x1[,1], x1)
  xog = cbind(x,x^2, x^3, x^4, x1, x*x1)
  
  truey = 0.5 + xmat%*%b0
  y = truey + rnorm(n)* sigma
  
  dat = list(truey = truey,
             y = y,
             x = xog,
             xmat = xmat)
}

genDat = function(b0, n, p, xmin=0,xmax=1,
                  sigma = 1/3, int = 0){
  
  # x = runif(n,-1,1)
  # x1 = matrix(runif(n*p,-1,1), nrow = n)
  # xmat = cbind(x, x^2, x^3, x^4, x1)
  xmat = matrix(runif(n*p, xmin,xmax), nrow = n)
  
  truey = int + xmat%*%b0
  y = truey + rnorm(n)* sigma
  
  dat = list(truey = truey,
             y = y,
             x = xmat)
}

boot_conform = function(n,n1, dattrain, dattest, 
                 replace = FALSE, nzero = NA,
                 oracle = FALSE, refit = FALSE, 
                 Int = TRUE){
  samp1 = sample(n, n1, replace = replace)
  
  tab1 = table(samp1)
  idx = rep(0,n)
  idx[as.numeric(names(tab1))] = tab1
  
  xmat = as.matrix(dattrain$x)
  
  if(oracle){
    m1 = glm(dattrain$y[samp1] ~ xmat[samp1,s1])
    coef1 = coef(m1)# [c(1,1+s1)]
    ypred = coef1[1] + dattest$x[,s1]%*% coef1[-1]
    
  }else {
    g1 = cv.glmnet(xmat[samp1,],dattrain$y[samp1])
    lam1 = ifelse(is.na(nzero), 
                  "lambda.1se",
                  g1$lambda[which(g1$nzero>=nzero)[1]])
    coef1 = coef(g1, s= lam1 )# [c(1,1+s1)]
    shat = which(coef1[-1]!=0)
    if(refit){
      samp2 = setdiff(1:n, samp1)
      g2 = lm(dattrain$y[samp2] ~ xmat[samp2,shat])
      coef2 = coef(g2)
      ypred = dattest$x[,shat]%*%coef2[-1] 
      if(Int) ypred = ypred + coef2[1]
      tab1 = table(samp2)
      idx = rep(0,n)
      idx[as.numeric(names(tab1))] = tab1
      
    }else{
      yval = predict(g1, xmat[-samp1,], s = lam1)
      if(!Int) yval = yval - coef(g1)[1]
      rval = abs(dattrain$y[-samp1] - yval)
      
      
      ypred = predict(g1, as.matrix(dattest$x), s= lam1)
      if(!Int) ypred = ypred - coef(g1)[1]
    }
    
  }
  
  return(list(J = idx,
              ypred = ypred,
              coeff =  as.numeric(coef1),
              rval = rval) 
         )
}

boot1 = function(n,n1, dattrain, dattest, 
                 replace = TRUE, nzero = NA,
                 oracle = FALSE, refit = FALSE, Int = TRUE){
  if(n1 == n-1){
    samp1 = setdiff(1:n,i)
  }else{
    samp1 = sample(n, n1, replace = replace)
    
  }
  tab1 = table(samp1)
  idx = rep(0,n)
  idx[as.numeric(names(tab1))] = tab1
  
  if(oracle){
    m1 = glm(dattrain$y[samp1] ~ dattrain$x[samp1,s1])
    coef1 = coef(m1)# [c(1,1+s1)]
    ypred = coef1[1] + dattest$x[,s1]%*% coef1[-1]
    
  }else {
    g1 = cv.glmnet(dattrain$x[samp1,],dattrain$y[samp1])
    lam1 = ifelse(is.na(nzero), 
                  "lambda.1se",
                  g1$lambda[which(g1$nzero>=nzero)[1]])
    coef1 = coef(g1, s= lam1 )# [c(1,1+s1)]
    shat = which(coef1[-1]!=0)
    if(refit){
      samp2 = setdiff(1:n, samp1)
      g2 = lm(dattrain$y[samp2] ~ dattrain$x[samp2,shat])
      coef2 = coef(g2)
      ypred = dattest$x[,shat]%*%coef2[-1] 
      if(Int) ypred = ypred + coef2[1]
      tab1 = table(samp2)
      idx = rep(0,n)
      idx[as.numeric(names(tab1))] = tab1
      
    }else{
      ypred = predict(g1, dattest$x, s= lam1)
      if(!Int) ypred = ypred - coef(g1)[1]
    }
    
  }
  
  return(list(J = idx,
              ypred = ypred,
              coeff =  as.numeric(coef1)) )
}

boot_scaled = function(n,n1, dattrain, dattest, 
                 replace = TRUE, nzero = NA,
                 oracle = FALSE,refit=FALSE, 
                 Int = FALSE, scaled = TRUE){
  if(n1 == n-1){
    samp1 = setdiff(1:n,i)
  }else{
    samp1 = sample(n, n1, replace = replace)
    
  }
  tab1 = table(samp1)
  idx = rep(0,n)
  idx[as.numeric(names(tab1))] = tab1
  
  if(oracle){
    m1 = glm(dattrain$y[samp1] ~ dattrain$x[samp1,s1])
    coef1 = coef(m1)# [c(1,1+s1)]
    ypred = coef1[1] + dattest$x[,s1]%*% coef1[-1]
    
  }else {
    g1 = cv.glmnet(dattrain$x[samp1,],dattrain$y[samp1])
    lam1 = ifelse(is.na(nzero), 
                  "lambda.1se",
                  g1$lambda[which(g1$nzero>=nzero)[1]])
    coef1 = coef(g1, s= lam1 )# [c(1,1+s1)]
    shat = which(coef1[-1]!=0)
    if(refit){
      samp2 = setdiff(1:n, samp1)
      g2 = lm(dattrain$y[samp2] ~ dattrain$x[samp2,shat])
      coef2 = coef(g2)
      ypred = dattest$x[,shat]%*%coef2[-1] 
      if(Int) ypred = ypred + coef2[1]
      tab1 = table(samp2)
      idx = rep(0,n)
      idx[as.numeric(names(tab1))] = tab1
      
    }else{
      ypred = predict(g1, dattest$x, s= lam1)
      if(!Int) ypred = ypred - coef(g1)[1]
    }
    
    ### error variance
    id_val = which(idx==0)
    yval = predict(g1, dattrain$x[id_val,], s= lam1)
    sige2val = mean((dattrain$y[id_val] - yval)^2)
    
    yfit = predict(g1, dattrain$x[-id_val,], s= lam1)
    sige2train = mean((dattrain$y[-id_val] - yfit)^2)
    
  }
  sige2scale = NA
  sige2scale2 = NA
  
  if(scaled){
    if(Int) ytrain = dattrain$y[samp1] else
      ytrain = dattrain$y[samp1] - int
    
    nzero = 50
    lam1 = ifelse(is.na(nzero), 
                  "lambda.1se",
                  g1$lambda[which(g1$nzero>=nzero)[1]])
    
    sfit = scalreg(dattrain$x[samp1,], ytrain,
                   lam0 = lam1, 
                   LSE=TRUE)
    
    sige2scale = sfit$hsigma
    sige2scale2 = sfit$lse$hsigma
    
  }
  
  return(list(J = idx,
              ypred = ypred,
              coeff =  as.numeric(coef1),
              sige2val = sige2val,
              sige2train = sige2train,
              sige2scale = sige2scale,
              sige2scale2 = sige2scale2
              ) )
}

regboot = function(n,n1, dattrain, dattest, 
                 replace = TRUE, nzero = NA,
                 oracle = FALSE){
  samp1 = sample(n, n1, replace = replace)
  tab1 = table(samp1)
  idx = rep(0,n)
  idx[as.numeric(names(tab1))] = tab1
  
  if(oracle){
    m1 = glm(dattrain$y[samp1] ~ dattrain$x[samp1,s1])
    coef1 = coef(m1)# [c(1,1+s1)]
    ypred = coef1[1] + dattest$x[,s1]%*% coef1[-1]
    
  }else {
    g1 = cv.glmnet(dattrain$x[samp1,],dattrain$y[samp1])
    lam1 = ifelse(is.na(nzero), 
                  "lambda.1se",
                  g1$lambda[which(g1$nzero>=nzero)[1]])
    ypred = predict(g1, dattest$x, s= lam1)
    coef1 = coef(g1, s= lam1 )# [c(1,1+s1)]
    
  }
  
  return(list(J = idx,
              ypred = ypred,
              coeff =  as.numeric(coef1)) )
}

extract<-function(fit1=list(),k,pos){
  temp1<-fit1[[k]][[pos]]
}

boottest = function(x, y, samp1){
  xs = bs(x[samp1])
  l1 = lm(y[samp1]~xs)
  coef1 = coef(l1)
  ypred = coef1[1] + xstest%*%coef1[-1]
  
  return(list("coef"=coef1,
              "ypred" = ypred))
}

varjack = function(predage,Ycount){
  n = ncol(Ycount)
  sum((apply((predage*(Ycount != 0)),2,
             function(a){mean(a[a!=0])}) - mean(predage))^2)*(n-1)/n
  
}
  
varfn<-function(predage,Ycount,n=NA,n1=NA,B=NA, 
                replace = TRUE, unbiased = FALSE,
                fac = 0.9){
  if(sum(is.na(c(n,n1,B))) > 0){
    n = ncol(Ycount)
    n1 = sum(Ycount[1,])
    B = nrow(Ycount)
  }
  if(replace){
    v1 = sum(cov(predage,Ycount)^2)
    ifelse(unbiased, v1 - n1/B*var(predage), v1)
    
  }else{
    v1 = (sum(cov(predage,Ycount)^2))*(n-1)*n/(n-n1)^2
    ifelse(unbiased, v1 - fac*n*n1/B/(n-n1)*var(predage), v1)
    
  }
  
}
