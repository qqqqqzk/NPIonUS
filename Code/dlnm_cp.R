mkXpred <-
  function(type, basis, at, predvar, predlag, cen) {
    #
    ################################################################################
    # CREATE THE MATRIX OF TRANSFORMED CENTRED VARIABLES (DEPENDENT ON TYPE)
    #
    # CREATE VECTORIZED LAGGED VALUES
    varvec <- if(is.matrix(at)) as.numeric(at) else rep(at,length(predlag))
    lagvec <- rep(predlag,each=length(predvar))
    #
    if(type=="cb") {
      # IF STANDARD CROSS-BASIS, CREATE MARGINAL BASIS AND CALL TENSOR
      # NB: ORDER OF BASIS MATRICES IN TENSOR CHANGED SINCE VERSION 2.2.4
      # CENTERING APPLIED ONLY MARGINALLY TO VAR DIMENSION
      basisvar <- do.call("onebasis",c(list(x=varvec),attr(basis,"argvar")))
      basislag <- do.call("onebasis",c(list(x=lagvec),attr(basis,"arglag")))
      if(!is.null(cen)) {
        basiscen <- do.call("onebasis",c(list(x=cen),attr(basis,"argvar")))
        basisvar <- scale(basisvar,center=basiscen,scale=FALSE)
      }
      Xpred <- tensor.prod.model.matrix(list(basisvar,basislag))
    } else if(type=="one") {
      # IF ONEBASIS, SIMPLY CALL THE FUNCTION WITH PROPER ARGUMENTS
      ind <- match(c("fun",names(formals(attr(basis,"fun")))),
                   names(attributes(basis)),nomatch=0)
      basisvar <- do.call("onebasis",c(list(x=varvec),attributes(basis)[ind]))
      if(!is.null(cen)) {
        basiscen <- do.call("onebasis",c(list(x=cen),attributes(basis)[ind]))
        basisvar <- scale(basisvar,center=basiscen,scale=FALSE)
      }
      Xpred <- basisvar
    } else {
      # FINALLY, IF GAM, CALL PredictMat WITH PROPER DATA
      # CENTERING APPLIED TO THE TENSOR PRODUCT (NOT EFFICIENT BUT EASIER)
      data <- list(varvec,lagvec)
      names(data) <- basis$term
      Xpred <- PredictMat(basis,data,n=length(varvec))
      if(!is.null(cen)) {
        data[[1]] <- rep(cen,length(varvec))
        cbcen <- PredictMat(basis,data,n=length(varvec))
        Xpred <- Xpred-cbcen
      }
    }
    #
    return(Xpred)
  }

mkcen <- function(cen, type, basis, range) {
  #
  ################################################################################
  #
  # IF NULL, TRY TO EXTRACT IT FROM BASIS
  if (nocen <- is.null(cen))
    cen <- switch(
      type,
      cb = attributes(basis)$argvar$cen,
      one = attributes(basis)$cen,
      gam = NULL
    )
  #
  # DEPENDING ON FUNCTION
  fun <- switch(
    type,
    cb = attributes(basis)$argvar$fun,
    one = attributes(basis)$fun,
    gam = basis$margin[[1]]$fun
  )
  # SET CENTERING DEPENDING ON FUNCTION AND cen TYPE
  if (!is.null(fun) && fun %in% c("thr", "strata", "integer", "lin")) {
    if (is.logical(cen)) cen <- NULL
  } else {
    # IF NULL OR TRUE, SET TO (APPROXIMATELY) MID-RANGE AND MESSAGE
    if (is.null(cen) || (is.logical(cen) && cen)) cen <- median(pretty(range))
    # IF FALSE, NULL
    if (is.logical(cen) && !cen) cen <- NULL
  }
  #
  # HOWEVER, IF INTERCEPT IS PRESENT, SET TO NULL
  int <- switch(
    type,
    cb = attributes(basis)$argvar$intercept,
    one = attributes(basis)$intercept,
    gam = basis$margin[[1]]$intercept
  )
  if (is.logical(int) && int) cen <- NULL
  #
  # MESSAGE
  if(nocen && !is.null(cen))
    message("centering value unspecified. Automatically set to ",cen)
  #
  return(cen)
}

seqlag <-
  function(lag,by=1) seq(from=lag[1],to=lag[2],by=by)


mkat <-
  function(at, from, to, by, range, lag, bylag) {
    #
    ################################################################################
    #
    # IF at IS NULL
    if(is.null(at)) {
      if(is.null(from)) from <- range[1]
      if(is.null(to)) to <- range[2]
      nobs <- ifelse(is.null(by),50,max(1,diff(range)/by))
      pretty <- pretty(c(from,to),n=nobs)
      pretty <- pretty[pretty>=from&pretty<=to]	
      at <- if(is.null(by)) pretty else seq(from=min(pretty),
                                            to=to,by=by)
      # IF at IS A MATRIX, CHECK AND NAME ROWS
    } else if(is.matrix(at)) {
      if(dim(at)[2]!=diff(lag)+1L)
        stop("matrix in 'at' must have ncol=diff(lag)+1")
      if(bylag!=1) stop("'bylag!=1 not allowed with 'at' in matrix form")
      if(is.null(rownames(at))) rownames(at) <- seq(nrow(at))
      # IF at IS A VECTOR
    } else at <- sort(unique(at))
    #
    return(at)
  }



getlink <-
  function(model, class, model.link=NULL) {
    #
    ################################################################################
    #
    # IF PROVIDED, JUST RETURN
    if(!is.null(model.link)) return(model.link)
    #
    # OTHERWISE, EXTRACT FROM MODEL (IF AVAILABLE)
    link <- if(all(class%in%c("lm")) || all(class%in%c("lme")) ||
               any(class%in%"nlme") || any(class%in%"lmerMod")) "identity" else 
                 if(any(class %in% c("clogit"))) "logit" else
                   if(any(class %in% c("coxph"))) "log" else
                     if(any(class %in% c("glm")) || any(class %in% c("glmmPQL")))
                       model$family$link else if(any(class %in% c("glmerMod")))
                         model@resp$family$link else NA
    #
    return(link)
  }

getvcov <-
  function(model, class) {
    #
    ################################################################################
    #
    # EXTRACT VCOV
    # NB: gam, gee AND geeglm HAVE CLASS glm AS WELL
    vcov <- if(any(class%in%c("lm","glm","lme","coxph")) &&
               !any(class%in%c("gee"))) vcov(model) else if(identical(class,c("gee","glm")))
                 model$robust.variance else if(any(class%in%c("geeglm")))
                   summary(model)$cov.scaled else if(any(class%in%c("lmerMod","glmerMod","lmerModLmerTest")))
                     as.matrix(vcov(model)) else tryCatch(vcov(model),error=function(w) "error")
    if(identical(vcov,"error")) stop("methods for coef() and vcov() must ",
                                     "exist for the class of object 'model'. If not, extract them manually and ",
                                     "use the arguments 'coef' and 'vcov'")
    #
    return(vcov)
  }

getcoef <-
  function(model, class) {
    #
    ################################################################################
    #
    # EXTRACT COEF
    # NB: gam, gee AND geeglm HAVE CLASS glm AS WELL
    coef <- if(any(class%in%c("glm","gam","coxph"))) coef(model) else
      if(any(class%in%c("lme","lmerMod","glmerMod","lmerModLmerTest"))) fixef(model) else
        tryCatch(coef(model),error=function(w) "error")
    if(identical(coef,"error")) stop("methods for coef() and vcov() must ",
                                     "exist for the class of object 'model'. If not, extract them manually and ",
                                     "use the arguments 'coef' and 'vcov'")
    #
    return(coef)
  }


crosspred <-
  function(basis, model=NULL, coef=NULL, vcov=NULL, model.link=NULL, at=NULL,
           from=NULL, to=NULL, by=NULL, lag, bylag=1, cen=NULL, ci.level=0.95,
           cumul=FALSE) {
    #
    ################################################################################
    # DETERMINE THE TYPE OF MODEL AND CHECKS
    #
    # TYPE OF PREDICTION: CROSSBASIS, ONEBASIS, OR PENALIZED GAM
    type <- if(any(class(basis)%in%"crossbasis")) "cb" else 
      if(any(class(basis)%in%"onebasis")) "one" else "gam"
    #
    # CHECKS ON TYPE, AND SET name, basis AND RESET type
    errormes <- "arguments 'basis' and 'model' not consistent. See help(crosspred)"
    if(type=="gam") {
      if(!is.character(basis) || length(basis)>1L) stop(errormes)
      if(is.null(model) || !any(class(model)%in%"gam")) stop(errormes)
      name <- basis
      sterms <- sapply(model$smooth, function(x) x$term[1])
      if(name%in%sterms) basis <- model$smooth[[which(sterms==name)[1]]] else 
        stop(errormes)
      if(length(which(sterms==name)) > 1)
        warning(paste(name,"included in multiple smoothers, only the first one taken"))
      if(!"cb.smooth"%in%class(basis) && basis$dim>1L)
        stop("predictions not provided for multi-dimensional smoothers other than 'cb'")
    } else name <- deparse(substitute(basis))
    #
    #  EXTRACT ORIGINAL lag (DEPENDENT ON TYPE)
    origlag <- switch(type,
                      cb = attr(basis,"lag"),
                      one = c(0,0),
                      gam = if(is.null(basis$lag)) c(0,0) else basis$lag
    )
    lag <- if(missing(lag)) origlag else mklag(lag)
    #
    # CHECKS ON lag AND bylag
    if(!all.equal(lag,origlag) && cumul)
      stop("cumulative prediction not allowed for lag sub-period")
    lagfun <- switch(type, cb=attr(basis,"arglag")$fun, one=NULL,
                     gam=if(basis$dim==1L) NULL else basis$margin[[2]]$fun)
    if(bylag!=1L && !is.null(lagfun) && lagfun=="integer")
      stop("prediction for non-integer lags not allowed for type 'integer'")
    #
    # OTHER COHERENCE CHECKS
    if(is.null(model) && (is.null(coef) || is.null(vcov)))
      stop("At least 'model' or 'coef'-'vcov' must be provided")
    if(!is.numeric(ci.level) || ci.level>=1 || ci.level<=0)
      stop("'ci.level' must be numeric and between 0 and 1")
    #
    ################################################################################
    # SET COEF, VCOV CLASS AND LINK FOR EVERY TYPE OF MODELS
    #
    # WRITE CONDITIONS (DEPENDENT ON TYPE AND IF MATRIX/VECTOR)
    cond <- if(type=="gam") with(basis,first.para:last.para) else
      if(ncol(basis)==1L) name else
        if(type=="one") paste(name,"b[0-9]{1,2}",sep="") else
          paste(name,"v[0-9]{1,2}\\.l[0-9]{1,2}",sep="")
    #
    # IF MODEL PROVIDED, EXTRACT FROM HERE, OTHERWISE DIRECTLY FROM COEF AND VCOV
    if(!is.null(model)) {
      model.class <- class(model)
      coef <- getcoef(model, model.class)
      vcov <- getvcov(model, model.class)
      indcoef <- if(type=="gam") cond else grep(cond,names(coef))
      indvcov <- if(type=="gam") cond else grep(cond,rownames(vcov))
      coef <- coef[indcoef]
      #print(indcoef)
      #print(cond)
      vcov <- vcov[indvcov,indvcov,drop=FALSE]
      model.link <- getlink(model, model.class, model.link)
    } else model.class <- NA
    #
    # CHECK
    npar <- if(type=="gam") length(indcoef) else ncol(basis)
    if(length(coef)!=npar || length(coef)!=dim(vcov)[1] || any(is.na(coef)) ||
       any(is.na(vcov)))
      stop("coef/vcov not consistent with basis matrix. See help(crosspred)")
    #
    ##########################################################################
    # AT, PREDVAR, PREDLAG AND CENTERING
    #
    # RANGE
    range <- if(type=="gam") range(model$model[[basis$term[1]]]) else
      attr(basis,"range")
    #
    # SET at, predvar AND predlag
    at <- mkat(at,from,to,by,range,lag,bylag)
    predvar <- if(is.matrix(at)) rownames(at) else at
    predlag <- seqlag(lag,bylag)
    #
    # DEFINE CENTERING VALUE (NULL IF UNCENTERED), AND REMOVE INFO FROM BASIS
    cen <- mkcen(cen, type, basis, range)
    if(type=="one") attributes(basis)$cen <- NULL
    if(type=="cb") attributes(basis)$argvar$cen <- NULL
    #
    ################################################################################
    # PREDICTION OF LAG-SPECIFIC EFFECTS
    #
    # CREATE THE MATRIX OF TRANSFORMED CENTRED VARIABLES (DEPENDENT ON TYPE)
    Xpred <- mkXpred(type,basis,at,predvar,predlag,cen)
    #
    # CREATE LAG-SPECIFIC EFFECTS AND SE
    matfit <- matrix(Xpred%*%coef, length(predvar), length(predlag)) 
    matse <- matrix(sqrt(pmax(0,rowSums((Xpred%*%vcov)*Xpred))), length(predvar),
                    length(predlag)) 
    #
    # NAMES
    rownames(matfit) <- rownames(matse) <- predvar
    colnames(matfit) <- colnames(matse) <- outer("lag",predlag,paste,sep="")
    #
    ################################################################################
    # PREDICTION OF OVERALL+CUMULATIVE EFFECTS
    #
    # RE-CREATE LAGGED VALUES (NB: ONLY LAG INTEGERS)
    predlag <- seqlag(lag)
    #
    # CREATE THE MATRIX OF TRANSFORMED VARIABLES (DEPENDENT ON TYPE)
    Xpred <- mkXpred(type,basis,at,predvar,predlag,cen)
    #
    # CREATE OVERALL AND (OPTIONAL) CUMULATIVE EFFECTS AND SE
    Xpredall <- 0
    if(cumul) {
      cumfit <- cumse <- matrix(0,length(predvar),length(predlag))
    }
    for (i in seq(length(predlag))) {
      ind <- seq(length(predvar)) + length(predvar)*(i-1)
      Xpredall <- Xpredall + Xpred[ind,,drop=FALSE]
      if(cumul) {
        cumfit[, i] <- Xpredall %*% coef
        cumse[, i] <- sqrt(pmax(0,rowSums((Xpredall%*%vcov)*Xpredall)))
      }
    }
    allfit <- as.vector(Xpredall %*% coef)
    allse <- sqrt(pmax(0,rowSums((Xpredall%*%vcov)*Xpredall)))
    #
    # NAMES
    names(allfit) <- names(allse) <- predvar
    if(cumul) {
      rownames(cumfit) <- rownames(cumse) <- predvar
      colnames(cumfit) <- colnames(cumse) <- outer("lag",seqlag(lag),paste,sep="")
    }
    #
    ################################################################################
    # CREATE THE OBJECT
    #
    # INITIAL LIST, THEN ADD COMPONENTS
    list <- list(predvar=predvar)
    if(!is.null(cen)) list$cen <- cen
    list <- c(list, list(lag=lag, bylag=bylag, coefficients=coef, vcov=vcov,
                         matfit=matfit, matse=matse, allfit=allfit, allse=allse))
    if(cumul) list <- c(list, list(cumfit=cumfit, cumse=cumse))
    #
    # MATRICES AND VECTORS WITH EXPONENTIATED EFFECTS AND CONFIDENCE INTERVALS
    z <- qnorm(1-(1-ci.level)/2)
    if(!is.null(model.link) && model.link %in% c("log","logit")) {
      list$matRRfit <- exp(matfit)
      list$matRRlow <- exp(matfit-z*matse)
      list$matRRhigh <- exp(matfit+z*matse)
      list$allRRfit <- exp(allfit)
      list$allRRlow <- exp(allfit-z*allse)
      names(list$allRRlow) <- names(allfit)
      list$allRRhigh <- exp(allfit+z*allse)
      names(list$allRRhigh) <- names(allfit)
      if(cumul) {
        list$cumRRfit <- exp(cumfit)
        list$cumRRlow <- exp(cumfit-z*cumse)
        list$cumRRhigh <- exp(cumfit+z*cumse)
      }
    } else {
      list$matlow <- matfit-z*matse
      list$mathigh <- matfit+z*matse
      list$alllow <- allfit-z*allse
      names(list$alllow) <- names(allfit)
      list$allhigh <- allfit+z*allse
      names(list$allhigh) <- names(allfit)
      if(cumul) {
        list$cumlow <- cumfit-z*cumse
        list$cumhigh <- cumfit+z*cumse
      }
    }
    #
    list$ci.level <- ci.level
    list$model.class <- model.class
    list$model.link <- model.link
    #
    class(list) <- "crosspred"
    #
    return(list)
  }

