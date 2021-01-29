rm(list=ls())
libraries_check <- c("data.table","xtable","splines","mvtnorm")
for (libs in libraries_check) {
  # if(!libs %in% rownames(installed.packages())) {
  #   install.packages(libs,repos="http://lib.ugent.be/CRAN/")
  # }
  library(libs, character.only=TRUE)
}
rm(libraries_check, libs)

# simulation settings
simsets <- expand.grid(n=1000,p=c(12,16),q=c(0,4,8),d=c(0,1))

# initialize for parallel MC jobs
args <- nrow(simsets)
if (!grepl("apple",sessionInfo()[[1]]$platform)) {
  args <- commandArgs(trailingOnly=TRUE) # for CMD BATCH '--args 1'  
  simsets <- simsets[rep(1:nrow(simsets),each=100),]
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

One_PopData <- function() {
  l <- matrix(rnorm(n=N*(p+q)),nrow=N,ncol=(p+q))
  # conditional probability of treatment: true propensity score
  g_raw <- runif(ncol(l),min=-1,max=1)
  linear.ps <- (l %*% g_raw)[,1]
  pA1 <- exp(linear.ps)/(1+exp(linear.ps))
  a <- rbinom(N,1,pA1)
  rm(linear.ps,pA1)
  # outcomes
  b_raw <- g_raw
  linear.out <- d*a + (l %*% b_raw)[,1]
  pY1 <- exp(linear.out)/(1+exp(linear.out))
  y <- rbinom(N,1,pY1)
  
  truedata <- data.frame("i"=1:N,"L"=l,"treat"=a,"Y"=y)
  if (q>0) {
    # select q confounders to be unmeasured
    u.idx <- as.integer(unlist(lapply(strsplit(tail(x=row.names(
      BackwardSelect(mydata=truedata,order_rule="min",return_est="CI")
      ),n=q),split="L."),"[",2)))
    u <- l[,u.idx,drop=FALSE]
    l <- l[,-u.idx]
    mydata <- data.frame("i"=1:N,"L"=l,"U"=u,"treat"=a,"Y"=y)
    if ("U" %in% colnames(mydata)) {
      colnames(mydata)[colnames(mydata)=="U"] <- "U.1"
    }
  } else {
    mydata <- data.frame("i"=1:N,"L"=l,"treat"=a,"Y"=y)
  }
  return(mydata)
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
popeffs <- OneDRstd_Est(Ls=trueLs, mydata=popdata, varest=TRUE)
popeffs
popeffs <- popeffs["ate"]

set.seed(NULL)
One_sim <- function(mc) {
  # random samples from the superpopulation
  mydata <- popdata[sort(sample(N,n,replace=FALSE)),]
  mydata[, "i"] <- 1:n
  
  # given confounders
  Lothers <- list(
    "none" = "1",
    "obs" = paste0("L.",1:p),
    "all.LU" = trueLs
  )
  # effect estimators
  res <- lapply(Lothers, function(l) {
    x <- OneDRstd_Est(Ls=l, mydata=mydata, varest=TRUE)
    ## standardized
    c("est"=x["ate"],
      "ci_l"=x["ate"] - qnorm(.975)*x["se"],
      "ci_u"=x["ate"] + qnorm(.975)*x["se"],
      "est_std"=as.numeric(x["ate"]/x["se"]))
  })
  res[["popeff"]] <- popeffs
  
  # different orderings of covariates
  Lordered <- BackwardSelect(mydata,order_rule="min",return_est="CI")
  colnames(Lordered) <- c("ate_std","ate_hat","ci_l","ci_u")
  row.names(Lordered) <- NULL
  
  # perturbed
  Lordered.perturb.min <- replicate(mc,BackwardSelect(
    mydata,order_rule="perturb.min",return_est="CI"))
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
  
  if (FALSE) {
    # 95% CI for effect
    pdf(paste0("sim1-example-plots-.pdf"),width=6,height=4)
    plot(res.Lordered$ate_hat[,"obs"],
         xlim=c(0,p+q_ex),ylim=range(res.Lordered[-1]),
         xlab="Number of confounders adjusted for",ylab="Effect estimate",
         main="Simulated data",type="n")
    abline(h=0,lty=3)
    for (ts in c("ci_l","ci_u")) {
      for (j in 2:(mc+1)) {
        lines(x=0:(p+q_ex),y=res.Lordered[[ts]][1:(p+1+q_ex),j],
              lwd=.01,col='grey')
      }
    }
    for (ts in c("ate_hat","ci_l","ci_u")) {
      points(x=0:p,y=res.Lordered[[ts]][1:(p+1),"obs"],col=1,
             pch=ifelse(ts=="ate_hat",1,ifelse(ts=="ci_u",2,6)))
      if (ts=="ate_hat") {
        lines(x=0:p,y=res.Lordered[[ts]][1:(p+1),"obs"])
        points(x=p+(q_ex),y=res$extrap[[ts]],pch=19)
      } else {
        for (i in q_ex) {
          points(x=p+i,
                 y=ifelse(ts=="ci_u",
                          yes=max(res.Lordered[[ts]][(p+1)+i,-1]),
                          no=min(res.Lordered[[ts]][(p+1)+i,-1])),
                 pch=ifelse(ts=="ci_u",24,25),bg=1)
        }
      }
    }
    points(x=p+q+0.5,y=res$all.LU["est.ate"],col=1,pch=21,bg="grey")
    points(x=p+q+0.5,y=res$all.LU["ci_u.ate"],col=1,pch=24,bg="grey")
    points(x=p+q+0.5,y=res$all.LU["ci_l.ate"],col=1,pch=25,bg="grey")
    dev.off()
  }
  
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
  
  return( c(unlist(simsets[seed,]), unlist(res)) )
  
}

n_sims <- 10
ptm=proc.time()[3]
sim_res <- replicate(n=n_sims,expr=One_sim(mc=500),simplify=FALSE)
proc.time()[3]-ptm
# 13 mins per sim (with 500 perturbations)

subfolder <- "sens-extrap-sim2/"
myfile <- "sens-extrap-sim-"
save(sim_res,file=paste0(subfolder,myfile,seed,".Rdata"))
cat(seed,"\n")
q()

# results =====================================================================
rm(list=ls())
library("data.table")
library("xtable")
subfolder <- "sens-extrap-sim2/"
myfiles <- list.files(subfolder)
myfiles <- myfiles[grep(pattern=".Rdata",myfiles)]
simres <- list()
for (ll in myfiles) {
  load(paste0(subfolder,ll))
  simres <- c(simres, sim_res)
  rm(sim_res)
  cat(ll,"\n")
}

sim_res <- data.table(do.call(rbind,simres))
setkey(sim_res)
# unique sims
sim_res[,length(unique(all.LU.est.ate)),by=list(n,p,q,d)]
# unique ATEs
sim_res[,unique(popeff.ate),by=list(n,p,q,d)]

# CIs and point estimates for effects =========================================
sim_res[, c("MSE.all.LU","MSE.obs","MSE.extrap",
            "CIcover.all.LU","CIcover.obs","CIcover.extrap",
            "PIcover","PIwidth") := list(
              sqrt((mean(all.LU.est.ate)-popeff.ate)^2 + var(all.LU.est.ate)),
              sqrt((mean(obs.est.ate)-popeff.ate)^2 + var(obs.est.ate)),
              sqrt((mean(extrap.ate_hat.obs)-popeff.ate)^2 + var(extrap.ate_hat.obs)),
              (all.LU.ci_l.ate <= popeff.ate) & (popeff.ate <= all.LU.ci_u.ate),
              (obs.ci_l.ate <= popeff.ate) & (popeff.ate <= obs.ci_u.ate),
              (extrap.ci_l.obs <= popeff.ate) & (popeff.ate <= extrap.ci_u.obs),
              (perturbed.ci_l1 <= popeff.ate) & (popeff.ate <= perturbed.ci_u2),
              perturbed.ci_u2-perturbed.ci_l1
  ), by=list(n,p,q,d)]
setkey(sim_res)

setkey(sim_res,q,d,p,n)
print(xtable(sim_res[, lapply(.SD, function(x) {
  if (length(unique(x))<=2) {
    format(round(mean(x), 2), nsmall = 2)
  } else {
    paste0(format(round(mean(x), 2), nsmall = 2)," (",
           format(round(sd(x), 2), nsmall = 2),")")  
  }
}), by=list(n,p,q,d)][
  ,c("q","popeff.ate","p",
     "all.LU.est.ate","obs.est.ate",
     "extrap.ate_hat.obs", # predicted effects
     "MSE.all.LU","MSE.obs","MSE.extrap",
     "CIcover.all.LU","CIcover.obs",
     "PIcover"# predicted inference
     ),
  with=FALSE],digits=0),
include.rownames=FALSE)
