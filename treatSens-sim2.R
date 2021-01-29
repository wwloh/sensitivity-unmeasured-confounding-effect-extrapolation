rm(list=ls())
libraries_check <- c("data.table","xtable","splines","mvtnorm","treatSens",
                     "sensemakr")
for (libs in libraries_check) {
  # if(!libs %in% rownames(installed.packages())) {
  #   install.packages(libs,repos="http://lib.ugent.be/CRAN/")
  # }
  library(libs, character.only=TRUE)
}
rm(libraries_check, libs)

# library("remotes")
# remotes::install_github("vdorie/treatSens")
# ?treatSens::treatSens

# simulation settings
simsets <- expand.grid(n=2000,p=c(16),q=c(0,4,8),d=c(0,2),ldist=c(0,1))

# initialize for parallel MC jobs
args <- nrow(simsets)
if (!grepl("apple",sessionInfo()[[1]]$platform)) {
  args <- commandArgs(trailingOnly=TRUE) # for CMD BATCH '--args 1'  
  simsets <- simsets[rep(1:nrow(simsets),each=200),]
  nrow(simsets)
}
(seed <- as.integer(args[1]))

# load helper functions
source("funs.R")

N <- 50000 # population
(n <- simsets[seed,"n"]) # sample size
(p <- simsets[seed,"p"]) # number of observed confounders
(q <- simsets[seed,"q"]) # number of unobserved confounders
(d <- simsets[seed,"d"]) # ATE
(ldist <- simsets[seed,"ldist"]) # logit (0) or probit (1) for exposure

One_PopData <- function() {
  l <- matrix(rnorm(n=N*(p+q)),nrow=N,ncol=(p+q))
  # conditional probability of treatment: true propensity score
  g_raw <- runif(n=ncol(l),min=-0.25,max=0.25)
  linear.ps <- (l %*% g_raw)[,1]
  if (ldist==0) {
    pA1 <- exp(linear.ps)/(1+exp(linear.ps))
  } else {
    pA1 <- pnorm(linear.ps) # probit regression for PS model
  }
  a <- rbinom(N,1,pA1)
  rm(linear.ps,pA1)
  # outcomes
  b_raw <- runif(n=ncol(l),min=-4,max=4)
  linear.out <- d*a + (l %*% b_raw)[,1]
  y <- linear.out + rnorm(N,sd=sqrt(ncol(l)))
  
  truedata <- data.frame("i"=1:N,"L"=l,"treat"=a,"Y"=y)
  if (q>0) {
    # select q confounders to be unmeasured
    u.idx <- as.integer(unlist(lapply(strsplit(tail(x=row.names(
      BackwardSelect(mydata=truedata,order_rule="min",return_est="CI",
                     Alogit=(ldist==0))
      ),n=q),split="L."),"[",2)))
    u <- l[,u.idx,drop=FALSE]
    l <- l[,-u.idx]
    mydata <- data.frame("i"=1:N,"L"=l,"U"=u,"treat"=a,"Y"=y)
    if ("U" %in% colnames(mydata)) {
      colnames(mydata)[colnames(mydata)=="U"] <- "U.1"
    }
    zetaz <- g_raw[u.idx]
    zetay <- b_raw[u.idx]
  } else {
    mydata <- data.frame("i"=1:N,"L"=l,"treat"=a,"Y"=y)
    zetaz <- zetay <- NA
  }
  return( list("zetaz"=zetaz,"zetay"=zetay,"data"=mydata) )
}

One_treatSens <- function(onedata,nsim_) {
  tSdata <- onedata$data
  X <- as.matrix(tSdata[,paste0("L.",1:p)])
  colnames(X) <- NULL
  Z <- tSdata[,"treat"]
  Y <- tSdata[,"Y"]
  # calibrated using observed covariates
  (y_fit.X <- coef(lm(Y~Z+X))[-(1:2)])
  (z_fit.X <- coef(glm(Z~X,binomial(link="probit")))[-1])
  if (q==0) {
    (zetay.cali <- y_fit.X[which.min(abs(y_fit.X))[1]])
    (zetaz.cali <- z_fit.X[which.min(abs(z_fit.X))[1]])
  } else {
    (zetay.cali <- y_fit.X[which.max(abs(y_fit.X))[1]])
    (zetaz.cali <- z_fit.X[which.max(abs(z_fit.X))[1]])
  }
  rm(y_fit.X,z_fit.X)
  
  # sensitivity analysis
  treatSens_OK <- tryCatch(
    out.bin <- treatSens(Y~Z+X, trt.family = binomial(link="probit"), 
                         nsim = nsim_, 
                         spy.range = c(zetay.cali,zetay.cali), 
                         spz.range = c(zetay.cali,zetaz.cali),
                         grid.dim = c(1,1)),
    error=function(cond) return(NA))
  if(all(is.na(treatSens_OK))) {
    tauhat_est <- rep(NA,4)
    names(tauhat_est) <- c("pt","ci_l","ci_u","std")
    return( tauhat_est )
  } else {
    tauhat <- matrix(NA,nrow=dim(out.bin$tau)[1],ncol=dim(out.bin$tau)[2])
    row.names(tauhat) <- row.names(out.bin$tau)
    colnames(tauhat) <- colnames(out.bin$tau)
    tauhat_est <- list()
    tauhat_est[["std"]] <- tauhat_est[["ci_u"]] <- tauhat_est[["ci_l"]] <- 
      tauhat_est[["pt"]] <- tauhat
    for (i in 1:nrow(tauhat)) {
      for (j in 1:ncol(tauhat)) {
        # average over all draws 
        tauhat.k <- out.bin$tau[i,j,]
        sigmahat.k <- out.bin$se.tau[i,j,]
        # SE estimate using Rubin's rules
        tauhat.se <- sqrt( mean(sigmahat.k^2) + 
                             (1+1/length(tauhat.k))*var(tauhat.k) )
        tauhat_est[["pt"]][i,j] <- mean(tauhat.k)
        tauhat_est[["ci_l"]][i,j] <- tauhat_est$pt[i,j] - qnorm(.975)*tauhat.se
        tauhat_est[["ci_u"]][i,j] <- tauhat_est$pt[i,j] + qnorm(.975)*tauhat.se
        tauhat_est[["std"]][i,j] <- tauhat_est$pt[i,j]/tauhat.se
      }
    }
  }
  return( unlist(tauhat_est) )
}

One_sensemakr <- function(mydata) {
  data.sensemakr <<- mydata[,c("Y","treat",paste0("L.",1:p))]
  # runs regression model
  obsL.model <- lm(Y~.,data = data.sensemakr)
  # benchmark using observed covariates
  (y_fit.X <- coef(obsL.model)[paste0("L.",1:p)])
  if (q==0) {
    (benchmark.L <- names(y_fit.X[which.min(abs(y_fit.X))[1]]))
  } else {
    (benchmark.L <- names(y_fit.X[which.max(abs(y_fit.X))[1]]))
  }
  reduceH0 <- (abs(d)<.Machine$double.eps)

  # runs sensemakr for sensitivity analysis
  sensemakr.fit <- sensemakr(model = obsL.model,
                             treatment = "treat",
                             benchmark_covariates = benchmark.L,
                             kd = 1, ky = 1, reduce=reduceH0)
  sensemakr.res <- sensemakr.fit$bounds[
    c("adjusted_estimate","adjusted_lower_CI", "adjusted_upper_CI",
      "adjusted_t")]
  names(sensemakr.res) <- c("pt","ci_l","ci_u","std")
  return( unlist(sensemakr.res) )
}

# generate data for the superpopulation
## confounder that is hidden remains the same for all observed samples
set.seed(9000)
popdata <- One_PopData()
# true standardized treatment effect
if (q>0) {
  trueLs <- c(paste0("L.",1:p),paste0("U.",1:q))  
} else {
  trueLs <- paste0("L.",1:p)
}
popeffs <- d

One_sim <- function(mc) {
  # random samples from the superpopulation
  mydata <- popdata$data[sort(sample(N,n,replace=FALSE)),]
  mydata[, "i"] <- 1:n
  
  onedata <- popdata
  onedata$data <- mydata
  res.treatSens <- One_treatSens(onedata=onedata,nsim_=mc)
  rm(onedata)
  
  res.sensemakr <- One_sensemakr(mydata)
  
  # given confounders
  Lothers <- list(
    "none" = "1",
    "obs" = paste0("L.",1:p),
    "all.LU" = trueLs
  )
  # effect estimators
  res <- lapply(Lothers, function(l) {
    x <- OneDRstd_Est(Ls=l, mydata=mydata, Alogit=FALSE, varest=TRUE)
    ## standardized
    c("est"=x["ate"],
      "ci_l"=x["ate"] - qnorm(.975)*x["se"],
      "ci_u"=x["ate"] + qnorm(.975)*x["se"],
      "est_std"=as.numeric(x["ate"]/x["se"]))
  })
  res[["popeff"]] <- popeffs
  
  # different orderings of covariates
  Lordered <- BackwardSelect(mydata,order_rule="min",return_est="CI",
                             Alogit=FALSE)
  colnames(Lordered) <- c("ate_std","ate_hat","ci_l","ci_u")
  row.names(Lordered) <- NULL
  
  # perturbed
  Lordered.perturb.min <- replicate(mc,BackwardSelect(
    mydata,order_rule="perturb.min",return_est="CI",Alogit=FALSE))
  row.names(Lordered.perturb.min) <- NULL
  Lordered_list <- lapply(1:ncol(Lordered), function(ts) {
    Lordered.perturb.min.ts <- Lordered.perturb.min[,ts,]
    colnames(Lordered.perturb.min.ts) <- paste(colnames(Lordered)[ts],1:mc,sep=".")
    return( cbind("obs"=Lordered[,ts],Lordered.perturb.min.ts) )
  })
  names(Lordered_list) <- colnames(Lordered)
  rm(Lordered)
  
  # extrapolated effects
  q_ex <- max(q,2)
  res.Lordered <- lapply(Lordered_list, function(Lordered) {
    apply(Lordered, 2, Extrapolate_StdEst_Lsequence, 
          ns.df=nrow(Lordered)-1,q_extrap=q_ex)
  })
  res[["extrap"]] <- lapply(res.Lordered, function(x) x[(p+1)+(q_ex),"obs"])
  
  # summaries of extrapolations based on perturbed sequences
  res.extrap.perturbed <- lapply(res.Lordered, function(ts) {
    range(ts[p+q_ex,-1]) # prediction intervals
  })
  ## how often absolute value of standardized estimate exceeds 2
  res.extrap.perturbed[["ate_std.absabove2"]] <- 
    mean(abs(res.Lordered$ate_std[p+q_ex,-1])>qnorm(.975))
  ## how often lower (upper) endpoint of CI is above (below) zero
  res.extrap.perturbed[["ci_l.above0"]] <- 
    mean(res.Lordered$ci_l[p+q_ex,-1]>0)
  res.extrap.perturbed[["ci_u.below0"]] <- 
    mean(res.Lordered$ci_u[p+q_ex,-1]<0)
  
  res[["perturbed"]] <- res.extrap.perturbed
  
  return( c(unlist(simsets[seed,]), unlist(res), 
            "treatSens"=res.treatSens, "sensemakr"=res.sensemakr) )
  
}

n_sims <- 5
ptm=proc.time()[3]
sim_res <- NULL
# do not use replicate function!! treatSens internally resets to fixed seed
for (ss in 1:n_sims) {
  rm(.Random.seed)
  set.seed(NULL)
  sim_res[[ss]] <- One_sim(mc=500)
}
sim_res <- do.call(rbind,sim_res)
proc.time()[3]-ptm
# 15 mins per sim (with 500 posterior draws and 500 perturbations)

subfolder <- "treatSens-sim2/"
myfile <- "treatSens-sim-"
save(sim_res,file=paste0(subfolder,myfile,seed,".Rdata"))
cat(seed,"\n")
q()

# results =====================================================================
rm(list=ls())
library("data.table")
library("xtable")
subfolder <- "treatSens-sim2/"
myfiles <- list.files(subfolder)
myfiles <- myfiles[grep(pattern=".Rdata",myfiles)]
simres <- list()
for (ll in myfiles) {
  load(paste0(subfolder,ll))
  simres <- c(simres, list(sim_res))
  rm(sim_res)
  cat(ll,"\n")
}

sim_res <- data.table(do.call(rbind,simres))
setkey(sim_res)
# unique sims
sim_res[,length(unique(all.LU.est.ate)),by=list(n,p,q,d,ldist)]
# unique ATEs
setnames(sim_res,old="popeff",new="popeff.ate")
sim_res[,unique(popeff.ate),by=list(n,p,q,d,ldist)]
# unsuccessful calibrated estimates
sim_res[, "NAs.treatSens.pt" := is.na(treatSens.pt)]
sim_res[, lapply(.SD, mean), by=list(n,p,q,d,ldist), 
        .SDcols=c("NAs.treatSens.pt")]

sim_res[, ldist := ifelse(test=ldist==0,yes="logit",no="probit")]

# CIs and point estimates for effects =========================================
sim_res[, c("MSE.all.LU","MSE.obs","MSE.extrap",
            "MSE.treatSens","MSE.sensemakr",
            "CIcover.all.LU","CIcover.obs","CIcover.extrap",
            "CIcover.treatSens","CIcover.sensemakr",
            "PIcover","PIwidth") := list(
              sqrt((mean(all.LU.est.ate)-popeff.ate)^2 + var(all.LU.est.ate)),
              sqrt((mean(obs.est.ate)-popeff.ate)^2 + var(obs.est.ate)),
              sqrt((mean(extrap.ate_hat.obs)-popeff.ate)^2 + var(extrap.ate_hat.obs)),
              sqrt((mean(treatSens.pt,na.rm=TRUE)-popeff.ate)^2 + 
                     var(treatSens.pt,na.rm=TRUE)),
              sqrt((mean(sensemakr.pt)-popeff.ate)^2 + var(sensemakr.pt)),
              (all.LU.ci_l.ate <= popeff.ate) & (popeff.ate <= all.LU.ci_u.ate),
              (obs.ci_l.ate <= popeff.ate) & (popeff.ate <= obs.ci_u.ate),
              (extrap.ci_l.obs <= popeff.ate) & (popeff.ate <= extrap.ci_u.obs),
              ifelse(test=is.na(treatSens.ci_l) || is.na(treatSens.ci_u),
                     yes=FALSE, no=(treatSens.ci_l <= popeff.ate) & 
                       (popeff.ate <= treatSens.ci_u)),
              (sensemakr.ci_l <= popeff.ate) & (popeff.ate <= sensemakr.ci_u),
              (perturbed.ci_l1 <= popeff.ate) & (popeff.ate <= perturbed.ci_u2),
              perturbed.ci_u2-perturbed.ci_l1
            ), by=list(n,p,q,d,ldist)]
setkey(sim_res)

setkey(sim_res,q,ldist,d,p,n)
print(xtable(sim_res[, lapply(.SD, function(x) {
  if (length(unique(x))<=2) {
    format(round(mean(x,na.rm=TRUE), 2), nsmall = 2)
  } else {
    paste0(format(round(mean(x,na.rm=TRUE), 2), nsmall = 2)," (",
           format(round(sd(x,na.rm=TRUE), 2), nsmall = 2),")")  
  }
}), by=list(n,p,q,d,ldist)][
  ,c("q","ldist","popeff.ate",
     "all.LU.est.ate","obs.est.ate",
     "NAs.treatSens.pt","treatSens.pt","sensemakr.pt",
     "extrap.ate_hat.obs"
  ),
  with=FALSE],digits=0),
include.rownames=FALSE)
print(xtable(sim_res[, lapply(.SD, function(x) {
  if (length(unique(x))<=2) {
    format(round(mean(x,na.rm=TRUE), 2), nsmall = 2)
  } else {
    paste0(format(round(mean(x,na.rm=TRUE), 2), nsmall = 2)," (",
           format(round(sd(x,na.rm=TRUE), 2), nsmall = 2),")")  
  }
}), by=list(n,p,q,d,ldist)][
  ,c("q","ldist","popeff.ate",
     "MSE.all.LU","MSE.obs",
     "MSE.treatSens","MSE.sensemakr",
     "MSE.extrap"
  ),
  with=FALSE],digits=0),
include.rownames=FALSE)
print(xtable(sim_res[, lapply(.SD, function(x) {
  if (length(unique(x))<=2) {
    format(round(mean(x,na.rm=TRUE), 2), nsmall = 2)
  } else {
    paste0(format(round(mean(x,na.rm=TRUE), 2), nsmall = 2)," (",
           format(round(sd(x,na.rm=TRUE), 2), nsmall = 2),")")  
  }
}), by=list(n,p,q,d,ldist)][
  ,c("q","ldist","popeff.ate",
     "CIcover.all.LU","CIcover.obs",
     "CIcover.treatSens","CIcover.sensemakr",
     "PIcover"# predicted inference
  ),
  with=FALSE],digits=0),
include.rownames=FALSE)