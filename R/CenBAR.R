#' @export CenBAR
#' @import MASS
#' @import mvtnorm
#' @import glmnet
#' @import splines
#' @import survival
#' @import cvTools
#' @import foreach
#' @import parallel


library (MASS)
library (mvtnorm)
library (glmnet)
library (splines)
library (survival)
require(cvTools)
require(foreach)
require(parallel)

ridge_regression<- function() {
  # ToDo
}
ridge_opt_beta <- function(X, Y, Yt, lambda.path, lambda.path0) {
  beta0s <- L0_ridge_regression(X, Y,lambda.path)
  n = ncol(as.matrix(beta0s))
  ret = list()
  minMse = 0
  beta_idx_minMse = 0
  betas = c()
  if (n > 1) {
    for (i in 1:n) {
      fit = cv5.L0(X, Y, Yt, lambda.path0, beta0=beta0s[,i])
      if (i == 1) {
        minMse = fit$mse
        beta_idx_minMse = i
        betas = fit$beta
      } else {
        if (fit$mse < minMse) {
          minMse = fit$mse
          beta_idx_minMse = i
        }
        betas = rbind(betas, fit$beta)
      }
    }
    ret$beta0 = beta0s[,beta_idx_minMse]
    ret$r_opt_lam = lambda.path[beta_idx_minMse]
    ret$betas = betas
  } else {
    ret$beta0 = beta0s
    ret$r_opt_lam = lambda.path
    fit = cv5.L0(X, Y, Yt, lambda.path0, beta0s)
    ret$betas = fit$beta
  }
  ret$beta0s=beta0s
  return(ret)
}

cv5.L0 <- function(X, Y, Yt, lambda.path, beta0=NULL, F=NULL) {
  lam_len = length(lambda.path)
  lam_opt = 2
  opt_idx = 1
  if ( lam_len > 1) {
    Ys = list()
    len = length(Y)
    folds <- cvFolds(len, K = 5, type = "random")
    res = list()
    amse = c()
    aterr = c()
    for (i in 1:5) {
      trainIdx = folds$subsets[folds$which != i]
      Xtrain = X[trainIdx, ]
      Xtest = X[-trainIdx, ]
      Ytrain = as.matrix(Y[trainIdx,])
      Ytest = as.matrix(Y[-trainIdx,])
      if ( !is.null(F)) {
        beta <- L0_ridge_regression(Xtrain, Ytrain,lambda.path)
      } else {
        beta = L0regression(Xtrain , Ytrain, lambda.path, beta0)$res
      }
      tryCatch({
        Yeluvate = Xtest %*% beta
      }, error = function(ex) {
        warning("An error was detected.\n")
      })
      mse = apply(Yeluvate, 2, function(y){ return(mean((y - Ytest)^2))})
      amse = rbind(amse, mse)
    }

    Eall = apply(amse, 2, mean)
    mMse = min(Eall)[1]
    opt_idx = which(Eall == mMse)
    lams_opt = lambda.path[opt_idx]
    lam_opt = lams_opt[1]

  } else {
    lam_opt = lambda.path
    lam1_opt = lambda.path
    opt_idx = 1
  }

  if (!is.null(F)) {
    lbeta <- L0_ridge_regression(X, Y, lam_opt)
    beta <- as.matrix(lbeta)
    times=0
  } else {
    if (ncol(as.matrix(beta0)) == 1 ) {
      ret = L0regression(X, Y, lam_opt, beta0)
      beta = ret$res
      times = ret$times
    } else {
      ret = L0regression(X, Y, lam_opt, beta0[,opt_idx])
      beta = ret$res
      times = ret$times
    }
  }

  if (lam_len == 1) {
    Yeluvate = X %*% beta
    mMse = mean((Yeluvate - Y)^2)
  }
  res = list()
  res$beta=beta
  res$opt_lam=lam_opt
  res$mse = mMse
  res
}

L0_crossvalidate <- function(X, Y, Yt, lambda.path, lambda.path1, beta0.method = c("ridge_cv", "ridge_opt")) {
  betas=c()
  beta0s=c()
  if (missing(beta0.method) || beta0.method == "ridge_cv") {
    res <- cv5.L0(X, Y, Yt, lambda.path, F=ridge_regression)
    beta0 = res$beta
    r_opt_lam = 0
  } else {
    myfit = ridge_opt_beta(X, Y, Yt, lambda.path=lambda.path1, lambda.path0=lambda.path)
    beta0 = myfit$beta0
    r_opt_lam=myfit$r_opt_lam
    betas = myfit$betas
    beta0s = myfit$beta0s
  }
  myfit1 <- cv5.L0(X, Y, Yt, lambda.path, beta0 = beta0,F=NULL)
  ret = list()
  ret$ret = myfit1
  ret$r_opt_lam = r_opt_lam
  ret$betas=betas
  ret$beta0s=beta0s
  return(ret)
}

L0_EMAlgorithm_nxn <- function(X, Y, beta0, lam, epn, eta, p=0, itn=500) {
  n = nrow(X)
  p = ncol(X)
  J = 1
  cont = 1
  md_beta = 0
  beta = beta0
  fn_beta = paste("L0_beta_dz_lam_", lam, ".csv", sep="")
  cnames = c(1:p)
  colnames(X)= cnames
  Xnz=X
  Ynz=Y
  while (cont) {
    zero_features = which(beta==0)
    zf_len=length(zero_features)
    if (nrow(beta) == zf_len) {
      break
    } else if (zf_len ==0 ) {
      u = beta
      Xnz=Xnz
    } else {
      u = beta[-zero_features,]
      Xnz=Xnz[,-zero_features]
    }
    Xnz=as.matrix(Xnz)
    if(length(u) > 1) {
      up = diag(as.vector(u^2))
    } else {
      up = u^2
    }

    tryCatch({
      if (nrow(Xnz) > 1) {
        D<- diag(nrow(Xnz))
      } else {
        D=1
      }
      Ku = Xnz %*% (up %*% t(Xnz)) + lam * D
    }, error = function(ex) {
      warning("An error was detected in norm\n")
    })
    beta = up%*%t(Xnz)%*%solve(Ku) %*% Ynz
    beta[abs(beta) < eta] = 0

    J = J + 1
    if ( J > itn) {
      warning("warning: too much tries")
      break
    }

    md_beta= norm(beta-u,"F")
    tryCatch({
      if (md_beta < epn) {
        cont = 0
      }
    }, error = function(ex) {
      warning("An error was detected in norm\n")
    })
  }
  beta[abs(beta) < eta] = 0
  if (length(beta) != 1 && ncol(Xnz) != length(beta)) {
    warning("warning: error Xnz & beta number!")
  }
  fNZ=as.numeric(colnames(Xnz))
  out_beta=rep(0,p)
  out_beta[fNZ]=beta
  as.matrix(out_beta)
  res = c()
  res$beta=out_beta
  res$times=J
  return(res)
}

L0_EMAlgorithm_d <- function(X, Y, beta0, lam, epn, eta, p=0, itn=500) {
  n = nrow(X)
  p = ncol(X)
  J = 1
  cont = 1
  md_beta = 0
  beta = beta0
  fn_beta = paste("L0_beta_dz_lam_", lam, ".csv", sep="")
  cnames = c(1:p)
  colnames(X)= cnames
  Xnz=X
  Ynz=Y
  while (cont) {
    zero_features = which(beta==0)
    zf_len=length(zero_features)
    if (nrow(beta) == zf_len) {
      break
    } else if (zf_len ==0 ) {
      u = beta
      Xnz=Xnz
    } else {
      u = beta[-zero_features,]
      Xnz=Xnz[,-zero_features]
    }
    up = 1/(u^2)
      if (length(up) > 1 ) {
        D = diag(as.vector(up), ncol(Xnz))
      } else {
        D = up
      }
      Ku = t(Xnz) %*% Xnz + lam * D
      beta = solve(Ku) %*% t(Xnz) %*% Ynz
    beta[abs(beta) < eta] = 0
    J = J + 1
    if ( J > itn) {
      warning("warning: too much tries")
      break
    }

    md_beta= norm(beta-u,"F")
    tryCatch({
      if (md_beta < epn) {
        cont = 0
      }
    }, error = function(ex) {
      warning("An error was detected in norm\n")
      message(md_beta)
      message(epn)
      message(beta)
      message(u)
    })
  }
  beta[abs(beta) < eta] = 0
  if (length(beta) != 1 && ncol(Xnz) != length(beta)) {
    warning("warning: error Xnz & beta number!")
  }
  fNZ=as.numeric(colnames(Xnz))
  out_beta=rep(0,p)
  out_beta[fNZ]=beta
  as.matrix(out_beta)
  res = c()
  res$beta=out_beta
  res$times=J
  return(res)
}

L0_evaluation <- function(X, Y, lambda.path, p=0, itn=500, epn, eta, theta0) {
  n = nrow(X)
  m = ncol(X)
  theta0 = as.matrix(theta0)
  nlams = length(lambda.path)
  res = array(0, dim=c(nrow(theta0), nlams))
  times = array(0, dim=c(nrow(theta0), nlams))
  if (nlams > 1) {
    wrapper <- function(i) {
      if (ncol(theta0) != 1 ) {
        stheta0 = theta0[,i]
      } else {
        stheta0 = theta0
      }
      lam = lambda.path[i]
      if (nrow(X) > ncol(X)) {
        theta <- L0_EMAlgorithm_d(X, Y, stheta0, lam, epn, eta)
      } else {
        theta <- L0_EMAlgorithm_nxn(X, Y, stheta0, lam, epn, eta)
      }
      return(theta)
    }
    mc <- getOption("mc.cores", 16)
    reh <- mclapply(1:nlams, wrapper, mc.cores = 1)
    for (i in 1:nlams) {
      res[,i] = reh[[i]]$beta
      times[,i] = reh[[i]]$times
    }
  } else {
    if (nrow(X) > ncol(X)) {
      reh = L0_EMAlgorithm_d(X, Y, theta0, lambda.path, epn, eta)
    } else {
      reh = L0_EMAlgorithm_nxn(X, Y, theta0, lambda.path, epn, eta)
    }
    res = reh$beta
    times = reh$times
    res[abs(res) < eta] = 0
  }
  ret = c()
  ret$res=res
  ret$times=times
  return(ret)
}

L0regression <- function(X, Y, lambda.path, stheta) {
  L0time = system.time(theta <- L0_evaluation(X, Y, lambda.path, p=0, eta=1e-4, epn=1e-5, theta0 = stheta))
  theta
}

BJASS <- function(x,y,delta,beta_ini,d,S,beta_S){
  iter = 1
  maxiter = 500
  tol = 1;
  u = max(eigen(t(x)%*%x)$values)
  beta_old = beta_ini
  a = dim(x)
  n = a[1]
  p = a[2]
  gamma = matrix(0,p,1)
  sc= rep(0,n)
  while (iter<= maxiter && tol>1E-3){
    e = y-x[,S] %*% beta_S
    z = sort(e)
    index = order(e)
    deltaz = delta[index]
    sur = survfit(Surv(e,delta)~1)
    km = sur$surv
    km[n] = km[n-1]
    sc[index] = km
    epp = ep(e,1-delta)
    B = 1*(matrix(rep(e,n),nrow=n,byrow=TRUE)>=matrix(rep(e,n),nrow=n))
    KMI = rowSums(matrix(rep(epp,n),nrow=n,byrow=TRUE)*B*matrix(rep(e,n),nrow=n,byrow=TRUE))
    ty = delta*y + (1-delta)*(x[,S] %*% beta_S+KMI/sc)
    s = ty - x %*% beta_old
    lold = -t(s) %*% s
    gamma = beta_old + 1/u * t(x) %*% s
    tmp = sort(abs(gamma),decreasing = TRUE)
    beta_new = gamma * (abs(gamma)>= tmp[d])
    lnew = -t(ty-x %*% beta_new) %*% (ty-x %*% beta_new)
    if (lnew< lold){
      u = 2*u
    }
    tol = sqrt(sum((beta_new-beta_old)^2))
    beta_old = beta_new
    iter = iter+1
  }

  S = which(beta_old != 0)
  beta_S = beta_old[S]
  return(list(S=S,beta_S=beta_S))
}

KSV <- function(x,y,delta,beta_ini,d){
  iter = 1
  maxiter = 500
  tol = 1;
  u = max(eigen(t(x)%*%x)$values)
  beta_old = beta_ini
  a = dim(x)
  n = a[1]
  p = a[2]
  gamma = matrix(0,p,1)
  sc= rep(0,n)
  while (iter<= maxiter && tol>1E-3){
    z = sort(y)
    index = order(y)
    deltaz = delta[index]
    sur = survfit(Surv(y,1-delta)~1)
    km = sur$surv
    km[n] = km[n-1]
    sc[index] = km
    ty = delta/sc
    s = ty*y - x %*% beta_old
    lold = -t(s) %*% s
    gamma = beta_old + 1/u * t(x) %*% s
    tmp = sort(abs(gamma),decreasing = TRUE)
    beta_new = gamma * (abs(gamma)>= tmp[d])
    lnew = -t(ty*y-x %*% beta_new) %*% (ty*y-x %*% beta_new)
    if (lnew< lold){
      u = 2*u
    }
    tol = sqrt(sum((beta_new-beta_old)^2))
    beta_old = beta_new
    iter = iter+1
  }
  S = which(beta_old != 0)
  beta_S = beta_old[S]
  return(list(S=S,beta_S=beta_S))
}

ep <- function(y,delta)
{
  n <- length(y);
  z = sort(y)
  index = order(y)
  deltaz = delta[index]
  delta_cen <- 1-deltaz;
  tmp1 =  delta_cen/rowSums(matrix(rep(z,n),nrow=n,byrow=TRUE)>=matrix(rep(z,n),nrow=n))
  tmp2 = matrix(1,n,1)-tmp1
  tmp3 = c(1,tmp2[1:n-1])
  ep = tmp1 * cumprod(tmp3)
  ep[index] = ep
  return(ep)
}

lambdaMax <- function(X, Y) {
  lam = apply(X, 2, function(x) ((t(x)%*%Y)^2 / (4*t(x) %*% x)))
  max(lam, na.rm=TRUE)
}

L0_ridge_regression <- function(x, y, lambda.path) {
  colnames(y) <- c("y")
  data <- cbind(y, x)
  fit <- lm.ridge(y~. , as.data.frame(data), lambda=lambda.path)
  theta0 = fit$coef
  theta0
}

inputData <- function(X_raw,Y_raw,delta_raw,enableScreening) {
  X0=X_raw
  Y0=Y_raw
  delta = delta_raw
  samples_cnt=nrow(X0)
  survival<- data.frame(X0, Y0, delta)
  survival.sort0<- survival[order(survival$Y0),]
  my.fit0 <- survfit(Surv(survival.sort0$Y0,1-survival.sort0$delta) ~ 1, data = survival.sort0)
  f.new0 <- my.fit0$surv
  B<- diag(length(f.new0))
  B[lower.tri(B)]=1
  f00<- as.matrix(f.new0)
  f0<- c(1,f00[1:length(f00)-1])
  f1<-((f0)^-1)
  f1[!is.finite(f1)] = 0
  T0<- as.matrix(survival.sort0$Y0)
  TT<- diff(c(0,T0))
  Ynew<- B%*%(TT*f1)
  dn <-  ncol(survival.sort0)
  Xnew <- as.matrix(survival.sort0[, -c(dn, dn-1)])
  delta = survival.sort0$delta

  X0=Xnew
  Y0=log(Ynew)
  S=c(1:ncol(X0))
  if (enableScreening) {
    a = dim(X0)
    n = a[1]
    p = a[2]
    d = floor(2*(n/log(n)))
    mr = 1 - sum(delta)/n
    beta_ini = matrix(0, p, 1)
    result = KSV(X0,Y0,delta,beta_ini,d)
    S = result$S
    beta_S = result$beta_S
    for (j in 1:5) {
      result1 = BJASS(X0,Y0,delta,beta_ini,d,S,beta_S)
      S= result1$S
      beta_S = result1$beta_S
    }
    X0<- X0[,S]
  }

  Xnew=X0
  Ynew=Y0
  Ynew=as.matrix(Ynew)
  delta=as.matrix(delta)
  Xnew <- scale(Xnew)
  Ynew <- scale(Ynew,scale = F)
  res = list()
  res$Xt <- Xnew
  res$Yt <- Ynew
  res$delta = delta
  res
}

CBAR_AFT <- function(X,Y,delta,lambda.path=NULL, enableScreening=FALSE) {
  rawdata <- inputData(X,Y,delta,enableScreening)
  X = rawdata$Xt
  Y = rawdata$Yt
  dataname <-colnames(X)
  delta=rawdata$delta
  residual <- rep(0,nrow(X))
  T = Y
  c = T
  Y = Y
  survival.data<- data.frame(X, Y, T, c, residual, delta)
  survival.sort0<- survival.data[order(survival.data$residual),]
  dn <-  ncol(survival.sort0)

  if(is.null(lambda.path)) {
    lmin=1e-4
    lamdMax <- lambdaMax(X, Y)
    lam = apply(X, 2, function(x) ((t(x)%*%Y)^2 / (4*t(x) %*% x)))
    lamdMax = max(lam, na.rm=TRUE)
    lambda.path <- exp(seq(log(lmin), log(lamdMax),l=50))
  } else {
    lamdMax = max(lambda.path)
  }

  lambda.path1 <- seq(1, 20, l=8)
  res = list()
  Xt = X
  Yt = Y
  set.seed(1000)
  myfit_L0 <- L0_crossvalidate(Xt, Yt, Yt, lambda.path, lambda.path1, beta0.method = "ridge_opt")
  res$CBAR = myfit_L0
  res$lamMax = lamdMax
  res$dataname=dataname
  res
}

CenBAR <- function(X,Y,delta,lambda.path=NULL, enableScreening=FALSE) {
  res <- CBAR_AFT(X,Y,delta,lambda.path, enableScreening)
  return(res)
}


