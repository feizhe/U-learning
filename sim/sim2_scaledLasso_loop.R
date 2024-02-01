# TrainData = list()

simname = paste0("sim_scaled_n=",n,"_n1=",n1,"_p=",p,
                 "_p0=",p0,"_replace=",replace,
                 "_scaled=",scaled)
print(simname)


nk = 1
predmat = coefmat = sampmat = sige2mat = NULL
nsim = 101
t1 = proc.time()
while (nk<= nsim) {
  if(nsc != 1) dattrain = TrainData[[nk]] else {
    dattrain = genDat(b0, n=n, p=p, sigma = sigm1,
                      xmin = -1, xmax = 1, int = int)
    TrainData[[nk]] = dattrain
  }
    
  
  fit1 <- foreach(i = 1:nboot, .packages=c("MASS","glmnet","scalreg")) %dopar% 
    boot_scaled(n,n1=n1, dattrain, dattest, 
                replace = replace, # nzero = 10,
                oracle = oracle, refit = refit, Int = Int, scaled = scaled)
  
  J <- sapply(1:nboot,extract,fit1=fit1,pos=1)
  ypred <- sapply(1:nboot,extract,fit1=fit1,pos=2)
  coeff <- sapply(1:nboot,extract,fit1=fit1,pos=3)
  
  sige2val = sapply(1:nboot,extract,fit1=fit1,pos=4)
  sige2train = sapply(1:nboot,extract,fit1=fit1,pos=5)
  sige2scale = sapply(1:nboot,extract,fit1=fit1,pos=6)
  sige2scale2 = sapply(1:nboot,extract,fit1=fit1,pos=7)
  
  coefmat = abind(coefmat, coeff, along = 3)
  predmat = abind(predmat, ypred, along = 3)
  sampmat = abind(sampmat, J, along = 3)
  
  sige2mat = abind(sige2mat, 
                   cbind(sige2val, sige2train, sige2scale, sige2scale2), 
                   along = 3)
  
  
  
  
  if(nk%%10==0){
    print(nk)
    print(dim(predmat))
    # print(dim(sampmat))
    print(proc.time() - t1)
    t1 = proc.time()
  }
  nk=nk+1
}
nk

coefmean = apply(coefmat,1,mean)
plot(b0, coefmean[-1])
abline(0,1)

if(nsc == 1) saveRDS(TrainData, paste0("TrainData_",simname, ".rds"))


# save.image(paste0(simname, ".RData"))


if(FALSE){
  #### smoothed variance
  m = nboot
  nsplit = dim(coefmat)[3]
  
  preds = SEs = SEs2 = SEs_jack =  
    bootSE = bootPred = COVER =  SIG_E2 = NULL
  for (k in 1:nsplit) {
    # idx = (1:m) + (k-1)*m
    
    sige2 = colMeans(sige2mat[,,k])
    
    pred1 = rowMeans(predmat[,,k])
    bootvar = apply(predmat[,,k], 1, var)
    bootpred = predmat[,,k][,sample(m,1)]
    
    vars = apply(predmat[,,k], 1, varfn, Ycount = t(sampmat[,,k]),
                 replace = replace, unbiased = TRUE)
    vars2 = apply(predmat[,,k], 1, varfn, Ycount = t(sampmat[,,k]),
                  replace = replace, unbiased = FALSE)
    
    varsjack = apply(predmat[,,k], 1, varjack, Ycount = t(sampmat[,,k]))
    
    
    preds = cbind(preds, pred1)
    SEs = cbind(SEs, sqrt(vars))
    SEs2 = cbind(SEs2, sqrt(vars2))
    SEs_jack =  cbind(SEs_jack, sqrt(varsjack))
    
    bootSE = cbind(bootSE, sqrt(bootvar))
    bootPred = cbind(bootPred, bootpred)
    SIG_E2 = rbind(SIG_E2, sige2)
    
    covs = NULL
    for (ci in 1:2) {
      for (cj in 1:3) {
        
        cov1 = as.numeric(
          (pred1 + 1.96*sqrt(switch(ci,vars, vars2) + sige2[cj]) - dattest$y )*
            (pred1 - 1.96*sqrt(switch(ci,vars, vars2) + sige2[cj]) - dattest$y ) < 0
        )
        
        covs = cbind(covs, cov1)
        
      }
    }
    # colMeans(covs)
    
    COVER = abind(COVER, covs, along = 3)
  }
  
  EmpSD = apply(preds, 1, sd,na.rm = TRUE)
  
  covprob = apply(COVER, c(1,2), mean)
  
  res2 = data.frame(
    ytrue = dattest$truey,
    ytest = dattest$y,
    pred = rowMeans(preds,na.rm = TRUE),
    bootpred = rowMeans(bootPred,na.rm = TRUE),
    EmpSD,
    bootEmpSD = apply(bootPred, 1, sd,na.rm = TRUE),
    SE = rowMeans(SEs,na.rm = TRUE), ## unbiased = TRUE
    SE2 = rowMeans(SEs2,na.rm = TRUE),  ## unbiased = FALSE
    SE_jack = rowMeans(SEs_jack,na.rm = TRUE),
    bootSE = rowMeans(bootSE,na.rm = TRUE),
    sig2val = colMeans(SIG_E2)[1],
    sig2train = colMeans(SIG_E2)[2],
    sig2scale = colMeans(SIG_E2)[3],
    cov1scaled = covprob[,3],
    cov2scaled = covprob[,6]
  )
  
  write.csv(res2, paste0(simname,".csv"))
  
  
  plist=list()
  for(i in 1:6){
    
    xr = range(res2$EmpSD,res2[,i+4]) + c(-0.01, 0.01)
    cor1 = round(cor(res2$EmpSD,res2[,i+4]),2)
    mae = round(mean(abs(res2$EmpSD - res2[,i+4])), 2)
    plist[[i]] = ggplot(res2, aes_string("EmpSD", colnames(res2)[i+4])) + 
      ggtitle(paste0("nboot=",nboot,", Cor = ",cor1,", MAE = ", mae)) + 
      xlab("Empirical SD") + ylab(colnames(res2)[i+4] ) + xlim(xr) + ylim(xr) +
      geom_point() + geom_abline(slope = 1, intercept = 0) +theme_bw()
    
  }
  # res2$pred = res2$pred -1
  xr = range(res2$ytrue,res2$pred) + c(-0.01, 0.01)
  cor1 = round(cor(res2$ytrue,res2$pred),2)
  mae = round(mean(abs(res2$ytrue - res2$pred)), 2)
  plist[[i+1]] = ggplot(res2, aes(ytrue, pred) ) + 
    ggtitle(paste0("nboot=",nboot,", Cor = ",cor1,", MAE = ", mae)) + 
    xlab("Truth") + ylab("Prediction" ) + xlim(xr) + ylim(xr) +
    geom_point() + geom_abline(slope = 1, intercept = 0) +theme_bw()
  
  
  pdf(paste0(simname, "_Var.pdf"))
  for(i in c(7,1:5)) print(plist[[i]])
  dev.off()
  
}

