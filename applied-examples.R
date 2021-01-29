rm(list=ls())
libraries_check <- c("data.table","speff2trial","CovSel","locfit","xtable",
                     "splines","BMA","bacr","mvtnorm")
for (libs in libraries_check) {
  if(!libs %in% rownames(installed.packages())) {
    install.packages(libs,repos="http://lib.ugent.be/CRAN/")
  }
  library(libs, character.only=TRUE)
}
rm(libraries_check, libs)

# load helper functions
source("funs.R")

# load processed data
args <- 1
args <- commandArgs(trailingOnly=TRUE) # for CMD BATCH '--args 1'
(seed <- as.integer(args[1]))
rm(args)
set.seed(9000*seed)

if (seed==1) {
  data.name="rhc"
  plot_main="RHC"
} else if (seed==2) {
  data.name="speff2"
  plot_main="ACTG175"
}
data.name

source(paste0("https://raw.githubusercontent.com/wwloh/stability-confounder-select/master/data-prep-",
              data.name,".R"))

# mydata: data.frame with 
## "i"=observation number
## "L"=covariates labelled L.1, ... L.p
## "treat"=observed treatment
## "Y"=observed outcome

(n <- nrow(mydata)) # sample size
(p <- sum(grepl("L",colnames(mydata)))) # number of observed confounders
(q <- ceiling(p/2)) # number of unobserved confounders
mc <- 500

# given confounders
Lothers <- list(
  "none" = "1",
  "obs" = paste0("L.",1:p)
)
# effect estimators
res <- lapply(Lothers, function(l) {
  x <- OneDRstd_Est(Ls=l, mydata=mydata, varest=TRUE)
  ## standardized
  as.numeric(x["ate"]/x["se"])
})

image_file <- paste0("applied-examples-res-",data.name,".Rdata")
if (file.exists(image_file)) {
  load(file=image_file)
} else {
  # different orderings of covariates
  Lordered <- BackwardSelect(mydata,order_rule="min",return_est="CI")
  colnames(Lordered) <- c("ate_std","ate_hat","ci_l","ci_u")
  ## save order in which covariates were eliminated
  eliminated_L <- list("obs"=row.names(Lordered))
  row.names(Lordered) <- NULL
  
  # perturbed
  Lordered.perturb.min <- lapply(1:mc, BackwardSelect, 
                                 mydata=mydata,
                                 order_rule="perturb.min",
                                 return_est="CI")
  eliminated_L <- c(eliminated_L,lapply(Lordered.perturb.min, row.names))
  Lordered_list <- lapply(1:ncol(Lordered), function(ts) {
    Lordered.perturb.min.ts <- do.call(cbind,lapply(Lordered.perturb.min,
                                                    function(x) x[,ts]))
    row.names(Lordered.perturb.min.ts) <- NULL
    colnames(Lordered.perturb.min.ts) <- paste(colnames(Lordered)[ts],1:mc,sep=".")
    return( cbind("obs"=Lordered[,ts],Lordered.perturb.min.ts) )
  })
  names(Lordered_list) <- colnames(Lordered)
  rm(Lordered)
  
  # extrapolated effects
  if (data.name=="rhc") {
    # use cross-validation to select number of knots to reduce overfitting ====
    ns_df.cv <- lapply(Lordered_list, function(Lordered) {
      nknots_criteria <- sapply(1:(nrow(Lordered)-1), function(nk) {
        extrap.pred <- apply(Lordered[,-1], 2, function(ates) {
          mean((predict(lm(ate~ns(m, df=nk), data=data.frame("m"=0:p,"ate"=ates)))-
                  Lordered[,1])^2)
        })
        return(mean(extrap.pred))
      })
      return(which.min(nknots_criteria))
    })
    print(ns_df.cv)
    ns_df <- ns_df.cv$ate_hat
  } else {
    ns_df <- nrow(Lordered)-1
  }
  res.Lordered <- lapply(Lordered_list, function(Lordered) {
    apply(Lordered, 2, Extrapolate_StdEst_Lsequence, 
          ns.df=ns_df,q_extrap=q)
  })
  res[["extrap"]] <- lapply(res.Lordered, function(x) x[(p+1)+(1:q),"obs"])
}

save.image(file=image_file)

# number of unmeasured confounders to plot  
q_ex <- ifelse(data.name=="rhc",yes=round(q/4),no=q)

pdf(paste0("applied-examples-plots-",data.name,".pdf"),
    width=ifelse(data.name=="rhc",9,6),height=ifelse(data.name=="rhc",6,4))
# 95% CI for effect
plot(res.Lordered$ate_hat[0:(p+q_ex),"obs"],
     xlim=c(0,p+q_ex),
     ylim=range(lapply(res.Lordered[-1], function(x) range(x[1:(p+1+q_ex),-1]))),
     xlab="Number of confounders adjusted for",ylab="Effect estimate",
     main=paste("Uncertainty intervals for",plot_main,"data"),
     type="n")
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
    points(x=p+(1:q_ex),y=res$extrap[[ts]][1:q_ex],pch=19)
    for (i in 1:q_ex) {
      lines(x=rep(p+i,2),y=quantile(res.Lordered$ate_hat[(p+1)+i,-1],
                                    probs=c(.025,.975),na.rm=TRUE))
    }
  } else {
    for (i in 1:q_ex) {
      points(x=p+i,
             y=ifelse(ts=="ci_u",
                      yes=quantile(res.Lordered[[ts]][(p+1)+i,-1],
                                   probs=.975,na.rm=TRUE),
                      no=quantile(res.Lordered[[ts]][(p+1)+i,-1],
                                  probs=.025,na.rm=TRUE)),
             pch=ifelse(ts=="ci_u",24,25),bg=1)
    }
  }
}
dev.off()
q()

# estimates adjusting for measured covariates
round(unlist(lapply(Lordered_list, function(x) x[p+1,"obs"])),digits=2)
rowMeans(sapply(1:(mc+1), function(j) {
  # CI excluding zero
  pmax((Lordered_list$ci_u[1:(p+1),j]<0),
       (Lordered_list$ci_l[1:(p+1),j]>0),na.rm=TRUE)
}))

# closer inspection of predicted effects and CIs based on perturbed sequences
## predicted effects that are negative
rowMeans(sapply(2:(mc+1), function(j) {
  res.Lordered$ate_hat[(p+1)+(1:q_ex),j]<0
}))
## lower endpoint of CI negative
rowMeans(sapply(2:(mc+1), function(j) {
  res.Lordered$ci_l[(p+1)+(1:q_ex),j]<0
}))
## CI includes zero
rowMeans(sapply(2:(mc+1), function(j) {
  pmax((res.Lordered$ci_u[(p+1)+(1:q_ex),j]<0),
       (res.Lordered$ci_l[(p+1)+(1:q_ex),j]>0),na.rm=TRUE)
}))

# check for systematic patterns between predicted effects and orbits ==========
## negative predicted effect
neg_pred <- (res.Lordered$ci_l[p+1+1,-1]<0)*1L
mean(neg_pred)
# orbit that each covariate was deleted in
orbit_eliminated <- do.call(cbind,lapply(paste0("L.",1:p), function (lj) {
  del_in_orbit_ <- as.factor(unlist(lapply(eliminated_L[-1], function (eL) {
    which(lj == eL)-1
  })))
}))
colnames(orbit_eliminated) <- paste0("L.",1:p)
orbit_eliminated <- data.frame(cbind(orbit_eliminated,neg_pred))
neg_pred_orbit_eliminated <- glm(neg_pred~., family=poisson("log"), 
                                 data=orbit_eliminated)
summary(neg_pred_orbit_eliminated)$coef
