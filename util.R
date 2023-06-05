library(huge)
library(pcalg)

#######################################
## GNARI with intercept
#######################################

simulate.GNARfit0 <- function(object, nsim=object$frbic$time.in, seed=NULL,
                              future=TRUE, set.noise=NULL, allcoefs=FALSE, ...){
  stopifnot(is.GNARfit(object))
  stopifnot(floor(nsim)==nsim)
  stopifnot(nsim>0)
  if(!is.null(seed)){
    stopifnot(floor(seed)==seed)
    set.seed(seed)
  }
  if(!is.null(object$frbic$fact.var)){
    stop("fact.var not currently supported with simulate")
  }
  dotarg <- list(...)
  if(length(dotarg)!=0){
    if(!is.null(names(dotarg))){
      warning("... not used here, input(s) ", paste(names(dotarg), collapse=", "), " ignored")
    }else{
      warning("... not used here, input(s) ", paste(dotarg, collapse=", "), " ignored")
    }
  }
  if(!is.null(set.noise)){
    sig <- set.noise
  }else{
    sig <- sigma(object$mod)
  }
  if(!allcoefs){
    nas <- is.na(object$mod$coefficients)
    pvs <- summary(object$mod)$coefficients[,4] < 0.05
    vals <- rep(0, length(pvs))
    vals[pvs] <- summary(object$mod)$coefficients[pvs,1]
    coefvec <- rep(0, length(nas))
    coefvec[(!nas)] <- vals
    
  }else{
    coefvec <- object$mod$coefficients
    coefvec[is.na(coefvec)] <- 0
  }
  
  #use GNARsim
  if(object$frbic$globalalpha){
    #global alpha has one alpha per time lag
    alphaout <-  vector(mode="list", length=object$frbic$alphas.in)
    betaout <- as.list(rep(0,length=object$frbic$alphas.in))
    count <- 1
    for(jj in 1:object$frbic$alphas.in){
      alphaout[[jj]] <- rep(coefvec[count], object$frbic$nnodes)
      if(object$frbic$betas.in[jj]>0){
        betaout[[jj]] <- coefvec[(count+1):(count+object$frbic$betas.in[jj])]
      }
      count <- count + object$frbic$betas.in[jj] + 1
    }
    
  }else{
    #multiple alphas per time lag
    alphaout <-  vector(mode="list", length=object$frbic$alphas.in)
    betaout <- as.list(rep(0,length=object$frbic$alphas.in))
    count <- 1
    for(jj in 1:object$frbic$alphas.in){
      alphaout[[jj]] <- coefvec[count:(count+object$frbic$nnodes-1)]
      if(object$frbic$betas.in[jj]>0){
        betaout[[jj]] <- coefvec[(count+object$frbic$nnodes):(count+
                                                                object$frbic$nnodes+object$frbic$betas.in[jj]-1)]
      }
      count <- count + object$frbic$nnodes + object$frbic$betas.in[jj]
    }
  }
  if(!future){
    newser <- GNARsim(n=nsim, net = object$frbic$net.in,
                      alphaParams = alphaout, betaParams = betaout,
                      sigma=sig)
  }else{
    nnodes <- object$frbic$nnodes
    max.nei <- max(unlist(lapply(betaout, length)))
    nei.mats <- vector(mode="list", length=max.nei)
    #create weight matrices for neighbours
    #flip network so that NofNeighbours gives into node information
    netmat <- as.matrix(object$frbic$net.in, normalise=FALSE)
    if(!isSymmetric(netmat)){
      net <- as.GNARnet(t(netmat))
    }else{
      net <- object$frbic$net.in
    }
    for(ii in 1:max.nei){
      nei.mats[[ii]] <- as.matrix(x=net, stage=ii, normalise=TRUE)
      if(sum(nei.mats[[ii]])==0){
        warning("beta order too large for network, neighbour set ",ii," is empty")
      }
    }
    
    xx.init <- object$frbic$final.in
    ntimedep <- object$frbic$alphas.in
    stopifnot(nrow(xx.init)==ntimedep)
    
    xx.gen <- matrix(NA, nrow=nsim+ntimedep, ncol=nnodes)
    xx.gen[1:ntimedep,] <- xx.init
    
    for(tt in (ntimedep+1):(nsim+ntimedep)){
      
      for(ii in 1:ntimedep){
        if(ii==1){
          time.forecast <- alphaout[[ii]]*xx.gen[(tt-ii),]
        }else{
          tmpi <- alphaout[[ii]]*xx.gen[(tt-ii),]
          time.forecast <- time.forecast + tmpi
        }
        
        
      }
      
      nei.forecast <- 0
      beta.pos <- NULL
      for(aa in 1:ntimedep){
        bb <- length(betaout[[aa]])
        if(bb>0){
          for(dd in 1:bb){
            nei.forecast <- nei.forecast + betaout[[aa]][dd]*xx.gen[tt-aa,]%*%t(nei.mats[[dd]])
          }
        }
      }
      xx.gen[tt,] <- time.forecast+nei.forecast+rnorm(n=object$frbic$nnodes, mean=0, sd=sig) + 
        object$mod$coefficients[(length(object$mod$coefficients)-object$frbic$nnodes+1):length(object$mod$coefficients)]
    }
    if(nsim==1){
      newser <- as.ts(t(xx.gen[(ntimedep+1):(nsim+ntimedep),]), start=1, end=nsim)
    }else{
      newser <- as.ts(xx.gen[(ntimedep+1):(nsim+ntimedep),], start=1, end=nsim)
    }
  }
  
  return(newser)
}

predict.GNARfit0 <- function(object, n.ahead=1, ...){
  stopifnot(is.GNARfit(object))
  if(!is.null(object$frbic$fact.var)){
    stop("fact.var not currently supported with predict")
  }
  return(simulate.GNARfit0(object, nsim=n.ahead, future=TRUE, set.noise=0, ...))
}



GNARfit0 <- function (vts = GNAR::fiveVTS, net = GNAR::fiveNet, alphaOrder = 2, 
                      betaOrder = c(1, 1), fact.var = NULL, globalalpha = TRUE, 
                      tvnets = NULL, netsstart = NULL, ErrorIfNoNei = TRUE) 
{
  stopifnot(is.GNARnet(net))
  stopifnot(ncol(vts) == length(net$edges))
  stopifnot(alphaOrder > 0)
  stopifnot(floor(alphaOrder) == alphaOrder)
  stopifnot(length(betaOrder) == alphaOrder)
  stopifnot(floor(betaOrder) == betaOrder)
  if (!is.null(fact.var)) {
    stopifnot(length(fact.var) == length(net$edges))
  }
  stopifnot(is.matrix(vts))
  stopifnot(is.logical(globalalpha))
  if (!is.null(tvnets)) {
    cat("Time-varying networks not yet supported")
  }
  stopifnot(is.null(tvnets))
  useNofNei <- 1
  frbic <- list(nnodes = length(net$edges), alphas.in = alphaOrder, 
                betas.in = betaOrder, fact.var = fact.var, globalalpha = globalalpha, 
                xtsp = tsp(vts), time.in = nrow(vts), net.in = net, 
                final.in = vts[(nrow(vts) - alphaOrder + 1):nrow(vts), 
                ])
  dmat <- GNARdesign(vts = vts, net = net, alphaOrder = alphaOrder, 
                     betaOrder = betaOrder, fact.var = fact.var, globalalpha = globalalpha, 
                     tvnets = tvnets, netsstart = netsstart)
  if (ErrorIfNoNei) {
    if (any(apply(dmat == 0, 2, all))) {
      parname <- strsplit(names(which(apply(dmat == 0, 
                                            2, all)))[1], split = NULL)[[1]]
      betastage <- parname[(which(parname == ".") + 1):(length(parname))]
      stop("beta order too large for network, use max neighbour set smaller than ", 
           betastage)
    }
  }
  predt <- nrow(vts) - alphaOrder
  yvec <- NULL
  for (ii in 1:length(net$edges)) {
    yvec <- c(yvec, vts[((alphaOrder + 1):(predt + alphaOrder)), 
                        ii])
  }
  tmp <- NULL
  tmp2 <- nrow(dmat) / length(net$edges)
  for(l in 1:length(net$edges)){
    tmp3 <- rep(0, nrow(dmat))
    tmp3[(tmp2*(l-1)+1) : (tmp2*l)] <- 1
    tmp <- cbind(tmp, tmp3)
  }
  colnames(tmp) <- paste("intercept", 1:length(net$edges), sep="")
  
  if (sum(is.na(yvec)) > 0) {
    yvec2 <- yvec[!is.na(yvec)]
    dmat2 <- dmat[!is.na(yvec), ]
    modNoIntercept <- lm(yvec2 ~ dmat2 + tmp + 0)
  }
  else {
    modNoIntercept <- lm(yvec ~ dmat + tmp + 0)
  }
  out <- list(mod = modNoIntercept, y = yvec, dd = dmat, frbic = frbic)
  class(out) <- "GNARfit"
  return(out)
}

NARIMAselect0 <- function(vts, net, max.alpha, max.beta, globalalpha = TRUE){
  minbic <- 10^5
  minalpha <- max.alpha
  for(i in 1:max.alpha){
    tryCatch({
      skip_to_next <<- FALSE
      modelbic <- BIC(GNARfit0(vts = vts, net = net, alphaOrder = i,
                               betaOrder = rep(0,i), globalalpha = globalalpha))
    }, error=function(e){skip_to_next <<- TRUE})
    if(skip_to_next) { next }
    
    if(!(is.nan(modelbic))){
      if(is.nan(minbic)){
        minbic <- modelbic
        minalpha <- i
      } else{
        if(minbic > modelbic){
          minbic <- modelbic
          minalpha <- i
        }
      }
    }
  }
  
  blist <- expand.grid(rep(list(0:max.beta), minalpha))
  minbic <- 10^5
  minbeta <- as.numeric(blist[nrow(blist),])
  
  for(i in 1:nrow(blist)){
    tryCatch({
      skip_to_next <<- FALSE
      modelbic <- BIC(GNARfit0(vts = vts,
                               net = net, alphaOrder = minalpha, betaOrder = as.numeric(blist[i,]), globalalpha = globalalpha))
    }, error=function(e){skip_to_next <<- TRUE})
    if(skip_to_next) { next }
    
    if(!(is.nan(modelbic))){
      if(is.nan(minbic)){
        minbic <- modelbic
        minbeta <- as.numeric(blist[i,])
      } else{
        if(minbic > modelbic){
          minbic <- modelbic
          minbeta <- as.numeric(blist[i,])
        }
      }
    }
  }
  
  res <- list()
  res$alphaOrder <- minalpha
  res$betaOrder <- minbeta
  
  return(res)
}



forecast_narima0 <- function(vts, h, N, net, max.alpha, max.beta, globalalpha = TRUE, centering=TRUE, strd=TRUE){
  if(centering){
    vmean <- rowMeans(t(vts)) 
    vts <- t(t(vts)-vmean)
  } 
  if(strd){
    vstd <- apply(vts,2,max)
    vts <- t(t(vts) / vstd)
  }
  narimaorder.gnar <- NARIMAselect0(vts = vts, net = net, max.alpha = max.alpha, max.beta = max.beta, globalalpha = globalalpha)
  print(narimaorder.gnar)
  model.fit <- GNARfit0(vts = vts, net = net, alphaOrder = narimaorder.gnar$alphaOrder, 
                        betaOrder = narimaorder.gnar$betaOrder, globalalpha = globalalpha)
  print(model.fit)
  pred.gnar <- predict.GNARfit0(model.fit, n.ahead=h)
  res <- t(pred.gnar)
  if(strd){
    res <- res * vstd
  }
  if(centering){
    res <- res + vmean
  } 
  return(res)
}

rmse <- function(x,y){
  return(sqrt(mean((x-y)^2)))
}


# Define CPDAG function ##
construct_CPDAG <- function(table, title, lambda=0.3, isnpn=TRUE,test_alpha=0.1,verbose=FALSE){
  if (isnpn){
    # Nonparanormal transform.
    table <- huge.npn(table, npn.func = "truncation", verbose = verbose)
  }
  district_graph <- huge(as.matrix(table), method = "glasso", lambda = lambda,verbose=verbose)
  
  district_graph_select <- huge.select(district_graph,verbose=verbose)
  
  
  g_suff = list(C=district_graph_select$opt.icov, n=nrow(table))
  pc_g.fit <- pc(suffStat = g_suff,
                 indepTest = gaussCItest, ## indep.test: partial correlations
                 alpha=test_alpha, labels = colnames(table), verbose = verbose,
                 solve.confl = TRUE,
                 fixedEdges = district_graph_select$refit
  )
  pc_g.fit
  par(family='AppleGothic',bg='white')
  plot(pc_g.fit, main=title,cex=1.5)
  
  return(pc_g.fit)
}
