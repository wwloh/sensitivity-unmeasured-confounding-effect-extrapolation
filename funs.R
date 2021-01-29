# function to calculate DR estimator of ATE
OneDRstd_Est <- function(Ls,mydata,Alogit=TRUE,varest=FALSE) {
  # PS model
  ps.Ls <- as.formula(paste0("treat~",paste(Ls,collapse="+")))
  alink <- ifelse(Alogit==TRUE,yes="logit",no="probit")
  ps.glm <- glm(ps.Ls,family=binomial(link = alink),data=mydata)  
  ps.glm.hat <- predict.glm(ps.glm,type="response")
  ipw.ps.glm <- mydata$treat/ps.glm.hat + (1-mydata$treat)/(1-ps.glm.hat)
  rm(ps.Ls,ps.glm,ps.glm.hat)
  
  # trim observations with extreme weights from sample
  ipw.ps.glm[ipw.ps.glm > 1/(.Machine$double.eps*1e1)] <- 0
  
  # Outcome model fitted via weighted regression
  out.Ls <- as.formula(paste0("Y~",paste(c("treat",Ls),collapse="+")))
  if (all(0<=mydata$Y & mydata$Y<=1)) {
    # logistic regression model
    out.lm <- glm(out.Ls,family=quasibinomial(link = "logit"),data=mydata,
                  weights=ipw.ps.glm)
  } else {
    # linear regression model
    out.lm <- lm(out.Ls,data=mydata,weights=ipw.ps.glm)
  }
  # for predicting counterfactuals
  mydata.A1 <- mydata.A0 <- mydata
  mydata.A1[,"treat"] <- 1
  mydata.A0[,"treat"] <- 0
  infn.est <- (2*mydata$treat-1)*ipw.ps.glm*
    (mydata$Y - predict(out.lm,type="response")) + # residuals
    predict(out.lm,newdata=mydata.A1,type="response") -
    predict(out.lm,newdata=mydata.A0,type="response")
  # one-step plug-in estimator
  ate.hat <- mean(infn.est)
  # individual influence functions
  infn.est <- infn.est - ate.hat
  
  if (varest==TRUE) {
    return( c("ate"=ate.hat,"se"=sqrt(var(infn.est)/length(infn.est))) )
  } else {
    return( list("ate"=ate.hat,"infn"=infn.est) )  
  }
}

# function to calculate DR estimator of ATE given treatment model coefficients
OneDRstd_Est_PSCoefs <- function(ps_coef,Ls,mydata,Alogit=TRUE,varest=FALSE) {
  # PS model
  a_x <- as.matrix(data.frame("(Intercept)"=1L,mydata[,Ls]))
  ps.glm.hat <- (a_x %*% ps_coef)[,1]
  if (Alogit==TRUE) {
    ps.glm.hat <- exp(ps.glm.hat)/(1+exp(ps.glm.hat))
  } else {
    ps.glm.hat <- pnorm(ps.glm.hat)
  }
  ipw.ps.glm <- mydata$treat/ps.glm.hat + (1-mydata$treat)/(1-ps.glm.hat)
  rm(a_x,ps.glm.hat)
  
  # trim observations with extreme weights from sample
  ipw.ps.glm[ipw.ps.glm > 1/(.Machine$double.eps*1e1)] <- 0
  
  # Outcome model fitted via weighted regression
  out.Ls <- as.formula(paste0("Y~",paste(c("treat",Ls),collapse="+")))
  if (all(0<=mydata$Y & mydata$Y<=1)) {
    # logistic regression model
    out.lm <- glm(out.Ls,family=quasibinomial(link = "logit"),data=mydata,
                  weights=ipw.ps.glm)
  } else {
    # linear regression model
    out.lm <- lm(out.Ls,data=mydata,weights=ipw.ps.glm)
  }
  # for predicting counterfactuals
  mydata.A1 <- mydata.A0 <- model.matrix(out.lm)
  mydata.A1[,"treat"] <- 1
  mydata.A0[,"treat"] <- 0
  # perturb outcome model coefficient values
  out_coef.mean <- out.lm$coefficients
  out_coef.sigma <- vcov(out.lm)
  out_betas <- rmvnorm(n=1,mean=out_coef.mean,sigma=out_coef.sigma)[1,]
  yhat <- (model.matrix(out.lm) %*% out_betas)[,1]
  y1hat <- (mydata.A1 %*% out_betas)[,1]
  y0hat <- (mydata.A0 %*% out_betas)[,1]
  if (all(0<=mydata$Y & mydata$Y<=1)) {
    yhat <- exp(yhat)/(1+exp(yhat))
    y1hat <- exp(y1hat)/(1+exp(y1hat))
    y0hat <- exp(y0hat)/(1+exp(y0hat))
  }
  infn.est <- (2*mydata$treat-1)*ipw.ps.glm*
    (mydata$Y - yhat) + # residuals
    y1hat - y0hat
  # one-step plug-in estimator
  ate.hat <- mean(infn.est)
  # individual influence functions
  infn.est <- infn.est - ate.hat
  
  if (varest==TRUE) {
    return( c("ate"=ate.hat,"se"=sqrt(var(infn.est)/length(infn.est))) )
  } else {
    return( list("ate"=ate.hat,"infn"=infn.est) )  
  }
}

BackwardSelect <- function(
  mydata, # observed data
  order_rule="min", # criterion for selecting confounder to be eliminated
  Alogit=TRUE, # logit or probit link for exposure model
  perturb.PS=FALSE, # whether to perturb exposure model or not
  return_est="standardized" # estimator to be returned
  ) {
  X.names <- grep("L",names(mydata),value=TRUE)
  p <- length(X.names)
  n <- nrow(mydata)
  X.curr <- X.names
  X.select <- NULL
  if (!grepl(pattern="perturb",x=order_rule)) {
    ate.curr <- OneDRstd_Est(Ls=X.curr,mydata,Alogit=Alogit,varest=FALSE)
    ate.none <- OneDRstd_Est(Ls="1",mydata,Alogit=Alogit,varest=FALSE)
  } else {
    # PS model
    alink <- ifelse(Alogit==TRUE,yes="logit",no="probit")
    ps.curr <- glm(as.formula(paste0("treat~",paste(X.curr,collapse="+"))),
                   family=binomial(link = alink),data=mydata)
    ps.none <- glm(treat~1,family=binomial(link = alink),data=mydata)
    if (perturb.PS==TRUE) {
      ## perturb PS model coef. values
      ps_betas <- rmvnorm(n=1,
                          mean=ps.curr$coefficients,
                          sigma=vcov(ps.curr))[1,]
      ps_nones <- rmvnorm(n=1,
                          mean=ps.none$coefficients,
                          sigma=vcov(ps.none))[1,]
    } else {
      ps_betas <- ps.curr$coefficients
      ps_nones <- ps.none$coefficients
    }
    # perturb ATE estimator
    ate.curr <- OneDRstd_Est_PSCoefs(
      ps_coef=ps_betas,Ls=X.curr,mydata,Alogit=Alogit,varest=FALSE)
    rm(ps_betas)
    ate.none <- OneDRstd_Est_PSCoefs(
      ps_coef=ps_nones,Ls=NULL,mydata,Alogit=Alogit,varest=FALSE)
  }
  ate.ordered <- list(ate.curr)
  while(length(X.curr) > 1) {
    # find covariate among those currently in the model
    ate.cands <- sapply(X.curr, function(X.cand) {
      X.nocand <- X.curr[!(X.curr %in% X.cand)]
      if (!grepl(pattern="perturb",x=order_rule)) {
        ate.cand <- OneDRstd_Est(Ls=X.nocand,mydata,Alogit=Alogit,varest=FALSE)
      } else {
        # PS model
        ps.xa <- glm(
          as.formula(paste0("treat~",paste(X.nocand,collapse="+"))),
          family=binomial(link = alink),data=mydata)
        if (perturb.PS==TRUE) {
          ## perturb PS model coef. values
          ps_betas <- rmvnorm(n=1,
                              mean=ps.xa$coefficients,
                              sigma=vcov(ps.xa))[1,]
        } else {
          ps_betas <- ps.xa$coefficients
        }
        # perturb ATE estimators
        ate.cand <- OneDRstd_Est_PSCoefs(
          ps_coef=ps_betas,Ls=X.nocand,mydata,Alogit=Alogit,varest=FALSE)
        rm(ps_betas)
        return( ate.cand )
      }
    })
    Xj.crit <- apply(ate.cands, 2, function(ate.cand) {
      # debiased effect estimator
      (ate.cand$ate-ate.curr$ate)^2-var(ate.cand$infn-ate.curr$infn)/n
    })
    # set NAs to Inf, and negative values to zero
    Xj.crit[is.na(Xj.crit)] <- Inf
    Xj.crit[Xj.crit<0] <- 0
    # randomly select one among all candidates that meet the same criterion
    if (grepl(pattern="max",x=order_rule)) {
      X.del <- sample(names(Xj.crit[Xj.crit==max(Xj.crit)]),1)
    } else if (grepl(pattern="min",x=order_rule)) {
      X.del <- sample(names(Xj.crit[Xj.crit==min(Xj.crit)]),1)
    }
    X.select <- c(X.del,X.select)
    X.curr <- X.curr[!(X.curr %in% X.del)]
    ate.curr <- ate.cands[,X.del]
    ate.ordered <- c(list(ate.curr),ate.ordered)
  }
  # append remaining covariate
  X.select <- c(X.curr,X.select)
  ate.ordered <- c(list(ate.none),ate.ordered)
  names(ate.ordered) <- c("none",X.select)
  ate.ordered <- do.call(rbind,lapply(ate.ordered, function(ate) {
    ate.se <- sqrt(var(ate$infn)/n)
    ate.std <- ate$ate/ate.se
    if (return_est=="standardized") {
      return( ate.std )
    } else if (return_est=="CI") {
      return( c(ate.std, ate$ate+c(0,-1,1)*qnorm(0.975)*ate.se) )
    }
  }))
  return(ate.ordered)
}

# Extrapolate standardized effect estimator for a given covariate sequence
Extrapolate_StdEst_Lsequence <- function(ates, ns.df, q_extrap=1) {
  p_ <- length(ates)-1
  extrap <- data.frame("m"=0:p_,"ate"=ates)
  stm <- tryCatch(system.time(
    extrap.lm <- lm(ate ~ ns(m, df=ns.df), data=extrap)
  )[3], error=function(cond) return(NA))
  if (!is.na(stm)) {
    pred.lm <- predict(extrap.lm, newdata=data.frame("m"=0:(p_+q_extrap)))  
  } else {
    pred.lm <- NA
  }
  return( as.numeric(pred.lm) )
}
