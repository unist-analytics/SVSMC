soft_thresholding <- function(par,listdata){
  #G <- c()
  #for(i in 1:5){
  
  betatrain<-listdata$betatrain
  betavalidation<-listdata$betavalidation
  
  tuningsvd <- function(betatrain,BETA) {
    trainset.nl <- betatrain[,-classattribute]
    train.svd <- svd(trainset.nl)
    D <- diag(train.svd$d)
    
    diag(D) <- diag(D) - BETA
    D[which(D<0)] <- 0
    
    
    SVD_train <- train.svd$u %*% D %*% (t(train.svd$v))
    SVD_train <- as.data.frame(SVD_train)
    
    
    return(SVD_train)
  }
  
  tuningsvd_app <- function(par){
    svd.lev <- c()
    for (s in 1:classnumber){
      svd.level <- tuningsvd(subset(betatrain,betatrain[,classattribute]==s),par[s])
      svd.lev <- rbind(svd.lev,svd.level)
    }
    return(svd.lev)
  }
  
  syntrain.kknn <- c()
  syntrain.kknn <- betatrain
  syntrain.kknn[,-classattribute] <- tuningsvd_app(par)
  syn.kknn <- kknn.ordinal(targetCandidate~., syntrain.kknn,betavalidation[,-classattribute],distance =2,k=7,kernel = "triangular",param=0.5)
  fit.syn <- fitted(syn.kknn)
  conf.syn <-table(Actual=betavalidation[,classattribute],"Predictions"= fit.syn)
  G_syn <- Gmetric(conf.syn)
  
  #G <- c(G,G_syn)
  #}
  return(G_syn)
}

dataapproximation <- function(dataframe,boundary) {
  trainset.nl <- dataframe[,!colnames(dataframe)%in%"targetCandidate"]
  train.svd <- svd(trainset.nl)
  D <- diag(train.svd$d)
  
  diag(D) <- diag(D) - boundary
  D[which(D<0)] <- 0
  
  
  SVD_train <- train.svd$u %*% D %*% (t(train.svd$v))
  SVD_train <- as.data.frame(SVD_train)
  
  return(SVD_train)
}

SVDtruncation <- function(bound){
  
  svd.lev <- c()
  for (i in 1:5){
    svd.level <- dataapproximation(train.kknn[train.kknn$targetCandidate==i,],bound[i])
    svd.lev <- rbind(svd.lev,svd.level)
  }
  
  return(svd.lev)
}

Gmetric<-function(x){
  d<-nrow(x)
  (prod(diag(x)/rowSums(x)))^(1/d)
}


tuningbeta_whole <- function(train.kknn,classattribute){
  beta <- c()
  for(i in 1:5){
    resample <- createDataPartition(train.kknn[,classattribute], p = .1)$Resample1
    betatrain <- train.kknn[-resample,]
    betavalidation <- train.kknn[resample,]
    
    #SVD
    tuningsvd <- function(betatrain,BETA) {
      trainset.nl <- betatrain[,-classattribute]
      train.svd <- svd(trainset.nl)
      D <- diag(train.svd$d)
      BETA <- round(BETA)
      
      if(BETA < dim(D)[1]) diag(D)[BETA:dim(D)[1]] <- diag(D)[BETA:dim(D)[1]]*0.0
      
      SVD_train <- train.svd$u %*% D %*% (t(train.svd$v))
      SVD_train <- as.data.frame(SVD_train)
      
      
      return(SVD_train)
    }
    Beta <- dim(train.kknn)[2]
    beta0 <- 0

    while (sum(abs(Beta-beta0))>1){
     
        beta0 <- Beta
        
        G <- c()
        for (vv in 1:(dim(train.kknn)[2]-2)){
          Beta <-  vv + 1 
          syntrain.kknn <- c()
          syntrain.kknn <- betatrain
          syntrain.kknn[,-classattribute] <-tuningsvd(syntrain.kknn,Beta)
          syn.kknn <- kknn.ordinal(targetCandidate~., syntrain.kknn,betavalidation[,-classattribute],distance =2,k=7,kernel = "triangular",param=0.5)
          fit.syn <- fitted(syn.kknn)
          conf.syn <-table(Actual=betavalidation[,classattribute],"Predictions"= fit.syn)
          G_syn <- Gmetric(conf.syn)*1 + sum(diag(conf.syn))/nrow(betavalidation)*0.0
          
          G <- c(G,G_syn)
        }
        ratio <-seq(dim(train.kknn)[2]-2)+1
        plot(ratio,G,type='l')
        Beta <- ratio[which.max(G)] 
      }  
    
    
    beta[i] <- list(Beta)
    
    
  }
  
  sumbeta <- 0
  for (k in 1:5){
    sumbeta <- sumbeta + unlist(beta[k])
  }
  
  return(round(sumbeta/5))
}




kknn.ordinal<-function (formula = formula(train), train, test, na.action = na.omit(), param=0.5,
          k = 7, distance = 2, kernel = "optimal", ykernel = NULL, 
          scale = TRUE, contrasts = c(unordered = "contr.dummy", ordered = "contr.ordinal")) {
  if (is.null(ykernel)) 
    ykernel = 0
  
  weight.y = function(l = 1, diff = 0) {
    k = diff + 1 # The k whoses value is one is the diffrent parameter from the given parameter k=7
    result = matrix(0, l, l)
    diag(result) = k
    for (i in 1:(k - 1)) {
      for (j in 1:(l - i)) {
        result[j, j + i] = k - i
        result[j + i, j] = k - i
      }
    }
    result
  }
  kernel <- match.arg(kernel, c("rectangular", "triangular", 
                                "epanechnikov", "biweight", "triweight", "cos", "inv", 
                                "gaussian", "rank", "optimal"), FALSE)
  ca <- match.call()
  response = NULL
  old.contrasts <- getOption("contrasts")
  options(contrasts = contrasts)
  formula = as.formula(formula)
  mf <- model.frame(formula, data = train)
  mt <- attr(mf, "terms")
  mt2 <- delete.response(mt)
  cl <- model.response(mf)
  d <- sum(attr(mt, "order"))
  if (is.ordered(cl)) {
    response <- "ordinal"
    lev <- levels(cl)
  }
  if (is.numeric(cl)) 
    response <- "continuous"
  if (is.factor(cl) & !is.ordered(cl)) {
    response <- "nominal"
    lev <- levels(cl)
  }
  if (distance <= 0) 
    stop("distance must >0")
  if (k <= 0) 
    stop("k must >0")
  learn <- model.matrix(mt, mf)
  valid <- model.matrix(mt2, test)
  m <- dim(learn)[1]
  p <- dim(valid)[1]
  q <- dim(learn)[2]
  ind <- attributes(learn)$assign
  d.sd <- numeric(length(ind)) + 1
  we <- numeric(length(ind)) + 1
  d.sd = apply(learn, 2, stats::var)
  for (i in unique(ind)) {
    d.sd[ind == i] = sqrt(mean(d.sd[ind == i]))
    we[ind == i] = 1/sum(ind == i)
  }
  we[d.sd == 0] = 0
  d.sd[d.sd == 0] = 1
  if (scale) {
    learn <- sweep(learn, 2L, d.sd, "/", check.margin = FALSE)
    valid <- sweep(valid, 2L, d.sd, "/", check.margin = FALSE)
  }
  ord = order(we * apply(learn, 2, sd), decreasing = TRUE)
  we = we[ord]
  learn = learn[, ord, drop = FALSE]
  valid = valid[, ord, drop = FALSE]
  Euclid <- FALSE
  if (distance == 2) 
    Euclid <- TRUE
  if (Euclid) 
    dmtmp <- .C("dmEuclid", as.double(learn), as.double(valid), 
                as.integer(m), as.integer(p), as.integer(q), dm = double((k + 
                                                                            1L) * p), cl = integer((k + 1L) * p), k = as.integer(k + 
                                                                                                                                   1), as.double(distance), as.double(we), dup = FALSE, 
                PACKAGE = "kknn")
  else dmtmp <- .C("dm", as.double(learn), as.double(valid), 
                   as.integer(m), as.integer(p), as.integer(q), dm = double((k + 
                                                                               1L) * p), cl = integer((k + 1L) * p), k = as.integer(k + 
                                                                                                                                      1), as.double(distance), as.double(we), dup = FALSE, 
                   PACKAGE = "kknn")
  D <- matrix(dmtmp$dm, nrow = p, ncol = k + 1)
  C <- matrix(dmtmp$cl, nrow = p, ncol = k + 1)
  maxdist <- D[, k + 1]
  maxdist[maxdist < 1e-06] <- 1e-06
  D <- D[, 1:k]
  C <- C[, 1:k] + 1
  #print(cl[C])         #####HERE
  CL <- matrix(cl[C], nrow = p, ncol = k)
  if (response != "continuous") {
    l <- length(lev) ## l is the number of levels of data in total
    weightClass <- matrix(0, p, l)
  }
  if (response == "continuous") {
    weightClass <- NULL
  }
  W <- D/maxdist
  W <- pmin(W, 1 - (1e-06))
  W <- pmax(W, 1e-06)
  if (kernel == "rank") 
    W <- (k + 1) - t(apply(as.matrix(D), 1, rank))
  if (kernel == "inv") 
    W <- 1/W
  if (kernel == "rectangular") 
    W <- matrix(1, nrow = p, ncol = k)
  if (kernel == "triangular") 
    W <- 1 - W
  if (kernel == "epanechnikov") 
    W <- 0.75 * (1 - W^2)
  if (kernel == "biweight") 
    W <- dbeta((W + 1)/2, 3, 3)
  if (kernel == "triweight") 
    W <- dbeta((W + 1)/2, 4, 4)
  if (kernel == "cos") 
    W <- cos(W * pi/2)
  if (kernel == "triweights") 
    W <- 1
  if (kernel == "gaussian") {
    alpha = 1/(2 * (k + 1))
    qua = abs(qnorm(alpha))
    W = W * qua
    W = dnorm(W, sd = 1)
  }
  if (kernel == "optimal") {
    W = rep(optKernel(k, d = d), each = p)
  }
  W <- matrix(W, p, k)
  if (response != "continuous") {
    for (i in 1:l) {
      weightClass[, i] <- rowSums(W * (CL == lev[i]))
    }
    weightClass <- weightClass/rowSums(weightClass)
    colnames(weightClass) <- lev

  }
  if (response == "ordinal") {
    blub = length(lev)
    weightClass = weightClass %*% weight.y(blub, ykernel)
    weightClass <- weightClass/rowSums(weightClass)
    weightClass <- t(apply(weightClass, 1, cumsum))
    colnames(weightClass) <- lev
    fit <- numeric(p)
    for (i in 1:p) fit[i] <- min((1:l)[weightClass[i, ] >= param])
    fit <- ordered(fit, levels = 1:l, labels = lev)
  }

  if (response == "continuous") 
    fit <- rowSums(W * CL)/pmax(rowSums(W), 1e-06)
  
  options(contrasts = old.contrasts)
  result <- list(fitted.values = fit, CL = CL, W = W, D = D, 
                 C = C, prob = weightClass, response = response, distance = distance, 
                 call = ca, terms = mt)
  class(result) = "kknn"
  result
}

kknn.cl<-function (formula = formula(train), train, test, na.action = na.omit(), param=0.5,
                        k = 7, distance = 2, kernel = "optimal", ykernel = NULL, 
                        scale = TRUE, contrasts = c(unordered = "contr.dummy", ordered = "contr.ordinal")) {
  if (is.null(ykernel)) 
    ykernel = 0
  
  weight.y = function(l = 1, diff = 0) {
    k = diff + 1 # The k whoses value is one is the diffrent parameter from the given parameter k=7
    result = matrix(0, l, l)
    diag(result) = k
    for (i in 1:(k - 1)) {
      for (j in 1:(l - i)) {
        result[j, j + i] = k - i
        result[j + i, j] = k - i
      }
    }
    result
  }
  kernel <- match.arg(kernel, c("rectangular", "triangular", 
                                "epanechnikov", "biweight", "triweight", "cos", "inv", 
                                "gaussian", "rank", "optimal"), FALSE)
  ca <- match.call()
  response = NULL
  old.contrasts <- getOption("contrasts")
  options(contrasts = contrasts)
  formula = as.formula(formula)
  mf <- model.frame(formula, data = train)
  mt <- attr(mf, "terms")
  mt2 <- delete.response(mt)
  cl <- model.response(mf)
  d <- sum(attr(mt, "order"))
  if (is.ordered(cl)) {
    response <- "ordinal"
    lev <- levels(cl)
  }
  if (is.numeric(cl)) 
    response <- "continuous"
  if (is.factor(cl) & !is.ordered(cl)) {
    response <- "nominal"
    lev <- levels(cl)
  }
  if (distance <= 0) 
    stop("distance must >0")
  if (k <= 0) 
    stop("k must >0")
  learn <- model.matrix(mt, mf)
  valid <- model.matrix(mt2, test)
  m <- dim(learn)[1]
  p <- dim(valid)[1]
  q <- dim(learn)[2]
  ind <- attributes(learn)$assign
  d.sd <- numeric(length(ind)) + 1
  we <- numeric(length(ind)) + 1
  d.sd = apply(learn, 2, stats::var)
  for (i in unique(ind)) {
    d.sd[ind == i] = sqrt(mean(d.sd[ind == i]))
    we[ind == i] = 1/sum(ind == i)
  }
  we[d.sd == 0] = 0
  d.sd[d.sd == 0] = 1
  if (scale) {
    learn <- sweep(learn, 2L, d.sd, "/", check.margin = FALSE)
    valid <- sweep(valid, 2L, d.sd, "/", check.margin = FALSE)
  }
  ord = order(we * apply(learn, 2, sd), decreasing = TRUE)
  we = we[ord]
  learn = learn[, ord, drop = FALSE]
  valid = valid[, ord, drop = FALSE]
  Euclid <- FALSE
  if (distance == 2) 
    Euclid <- TRUE
  if (Euclid) 
    dmtmp <- .C("dmEuclid", as.double(learn), as.double(valid), 
                as.integer(m), as.integer(p), as.integer(q), dm = double((k + 
                                                                            1L) * p), cl = integer((k + 1L) * p), k = as.integer(k + 
                                                                                                                                   1), as.double(distance), as.double(we), dup = FALSE, 
                PACKAGE = "kknn")
  else dmtmp <- .C("dm", as.double(learn), as.double(valid), 
                   as.integer(m), as.integer(p), as.integer(q), dm = double((k + 
                                                                               1L) * p), cl = integer((k + 1L) * p), k = as.integer(k + 
                                                                                                                                      1), as.double(distance), as.double(we), dup = FALSE, 
                   PACKAGE = "kknn")
  D <- matrix(dmtmp$dm, nrow = p, ncol = k + 1)
  C <- matrix(dmtmp$cl, nrow = p, ncol = k + 1)
  maxdist <- D[, k + 1]
  maxdist[maxdist < 1e-06] <- 1e-06
  D <- D[, 1:k]
  C <- C[, 1:k] + 1
  #print(cl[C])         #####HERE
  CL <- matrix(cl[C], nrow = p, ncol = k)
  if (response != "continuous") {
    l <- length(lev) ## l is the number of levels of data in total
    weightClass <- matrix(0, p, l)
  }
  if (response == "continuous") {
    weightClass <- NULL
  }
  W <- D/maxdist
  W <- pmin(W, 1 - (1e-06))
  W <- pmax(W, 1e-06)
  if (kernel == "rank") 
    W <- (k + 1) - t(apply(as.matrix(D), 1, rank))
  if (kernel == "inv") 
    W <- 1/W
  if (kernel == "rectangular") 
    W <- matrix(1, nrow = p, ncol = k)
  if (kernel == "triangular") 
    W <- 1 - W
  if (kernel == "epanechnikov") 
    W <- 0.75 * (1 - W^2)
  if (kernel == "biweight") 
    W <- dbeta((W + 1)/2, 3, 3)
  if (kernel == "triweight") 
    W <- dbeta((W + 1)/2, 4, 4)
  if (kernel == "cos") 
    W <- cos(W * pi/2)
  if (kernel == "triweights") 
    W <- 1
  if (kernel == "gaussian") {
    alpha = 1/(2 * (k + 1))
    qua = abs(qnorm(alpha))
    W = W * qua
    W = dnorm(W, sd = 1)
  }
  if (kernel == "optimal") {
    W = rep(optKernel(k, d = d), each = p)
  }
  W <- matrix(W, p, k)
  if (response != "continuous") {
    for (i in 1:l) {
      weightClass[, i] <- rowSums(W * (CL == lev[i]))
    }
    weightClass <- weightClass/rowSums(weightClass)
    colnames(weightClass) <- lev
    
  }
  if (response == "ordinal") {
    blub = length(lev)
    weightClass = weightClass %*% weight.y(blub, ykernel)
    weightClass <- weightClass/rowSums(weightClass)
    weightClass <- t(apply(weightClass, 1, cumsum))
    colnames(weightClass) <- lev
    fit <- numeric(p)
    for (i in 1:p) fit[i] <- min((1:l)[weightClass[i, ] >= param])
    fit <- ordered(fit, levels = 1:l, labels = lev)
  }
  
  if (response == "continuous") 
    fit <- rowSums(W * CL)/pmax(rowSums(W), 1e-06)
  
  options(contrasts = old.contrasts)
  result <- list(fitted.values = fit, CL = CL, W = W, D = D, 
                 C = C, prob = weightClass, response = response, distance = distance, 
                 call = ca, terms = mt)
  class(result) = "kknn"
  #result
  return(cl[C])
}


bounds_rm <- function(nrow,ncol,alpha,variance){
  #alpha from 0.01 to 0.05
  #variance 0.1 0.01 0.001
  sigma <- (variance)^.5
  #first bounds
  u = seq(0,3.09,by=0.01)
  p =pnorm(u)
  k <- sigma*min(which(p>(1-(alpha)/2)))*0.01
  first <- c(k,k*(nrow*ncol)^0.5)
  
  #second bounds
  second <- c(sigma*(qchisq(1-alpha,df=nrow))^.5,sigma*(ncol*qchisq(1-alpha,df=nrow))^.5)
  
  #Third bounds
  third <- c(sigma*(qchisq(1-alpha,df=nrow))^.5,sigma*(qchisq(1-alpha,df=nrow*ncol))^.5)
 

  
  return(c(k*(nrow*ncol)^0.5,sigma*(ncol*qchisq(1-alpha,df=nrow))^.5,sigma*(qchisq(1-alpha,df=nrow*ncol))^.5))
}


