Emt <- function(xyt,t.region,t.lambda,dt,kt="epanech",ht,correction="none",approach="simplified"){
  
  verifyclass(xyt,"stpp")
  
  correc <- c("none","isotropic","border","modified.border","translate","setcovf")
  id <- match(correction,correc,nomatch=NA)
  if (any(nbg <- is.na(id))){
    messnbg <- paste("unrecognised correction method:",paste(dQuote(correction[nbg]),collapse=","))
    stop(messnbg,call.=FALSE)
  }
  id <- unique(id)	
  correc2 <- rep(0,6)
  correc2[id] <- 1	
  
  appro <- c("simplified","standardised")
  im <- match(approach,appro,nomatch=NA)
  if (any(nbm <- is.na(im))){
    messnbm <- paste("unrecognised type of estimator:",paste(dQuote(approach[nbm]),collapse=","))
    stop(messnbm,call.=FALSE)
  }
  im <- unique(im)
  appro2 <- rep(0,2)
  appro2[im] <- 1
  
  ker <- c("box","epanech","biweight")
  ik <- match(kt,ker,nomatch=NA)
  if (any(nbk <- is.na(ik))){
    messnbk <- paste("unrecognised kernel function:",paste(dQuote(kt[nbk]),collapse=","))
    stop(messnbk,call.=FALSE)
  }
  ik <- unique(ik)
  ker2 <- rep(0,3)
  ker2[ik] <- 1
  
  dup <- duplicated(data.frame(xyt[,1],xyt[,2],xyt[,3]),fromLast = TRUE)[1]
  if (dup == TRUE){
    messnbd <- paste("spatio-temporal data contain duplicated points")
    warning(messnbd,call.=FALSE)
  }
  
  if (missing(ht)){
    d <- dist(xyt[,3])
    ht <- dpik(d,kernel=kt,range.x=c(min(d),max(d)))
  }
  
  if (missing(t.region)){
    tr <- range(xyt[,3],na.rm=TRUE)
    tw <- diff(tr)
    t.region <- c(tr[1]-0.01*tw,tr[2]+0.01*tw)
  }
  
  bsupt <- max(t.region)
  binft <- min(t.region)
  W <- sbox(xyt[,1:2], xfrac=0.01, yfrac=0.01)
  a <- diff(range(W[,1]))
  b <- diff(range(W[,2]))
  
  if (missing(dt)) {
    maxt <- (bsupt-binft)/4
    dt <- seq(ht,maxt,len=100)[-1]
  }
  if(dt[1]==0){
    dt <- dt[-1]
  }
  
  kernel <- c(kt=kt,ht=ht)
  Emttheo <- ((a^5+b^5-sqrt(a^2+b^2)*(a^4-3*(b*a)^2+b^4))/(15*(a*b)^2))+
    (((b^2)*log((a+sqrt(a^2+b^2))/b))/(6*a))+(((a^2)*log((b+sqrt(a^2+b^2))/a))/(6*b))
  
  pts <- xyt[,1:2]
  xytimes <- xyt[,3]
  ptsx <- pts[,1]
  ptsy <- pts[,2]
  ptst <- xytimes
  npt <- length(ptsx)
  ndt <- length(dt)
  snorm <- apply(pts,MARGIN=1,FUN=norm,type="2")
  eEmt <- rep(0,ndt)
  
  storage.mode(eEmt) <- "double"
  
  if (appro2[1]==1){
    Emtout <- .Fortran("emtcore",as.double(ptsx),as.double(ptsy),as.double(ptst),as.integer(npt),as.double(dt),
                       as.integer(ndt),as.integer(ker2),as.double(ht),(eEmt),PACKAGE="msfstpp")
    
    eEmt <- Emtout[[9]]
    Emt0 <- mean(snorm)

    invisible(return(list(eEmt=eEmt,Emt0=Emt0,dt=dt,kernel=kernel,Emttheo=Emttheo)))
  } else {
    
    if(missing(t.lambda)){
      misl <- 1
      t.lambda <- rep(npt/(bsupt-binft),npt)
    } else {
      misl <- 0
      if (length(t.lambda)==1){
        t.lambda <- rep(t.lambda,npt)
      }
    }
    
    wrt <- array(0,dim=c(npt,npt))
    wtt <- array(0,dim=c(npt,npt))
    wbit <- array(0,dim=c(npt,ndt))
    wbimodt <- array(0,dim=c(npt,ndt))
    wst <- rep(0,ndt)
    
    # correction="isotropic"
    if(correction=="isotropic"){
      wist <- tedgeRipley(ptst,binft,bsupt)
      wrt <- 1/wist
    }
    
    # correction="translate"
    if(correction=="translate"){
      wtrat <- tedgeTrans(ptst,t.region)
      wtt <- 1/wtrat
    }
    
    #  correction=="border" or "modified border"
    if(any(correction=="border")|any(correction=="modified.border")){
      bj = .bdist.times(xytimes, t.region)
      for(j in 1:ndt) { 
        wbit[,j] <- (bj>dt[j])/sum((bj>dt[j])/t.lambda)
        wbimodt[,j] <- (bj>dt[j])/.eroded.areat(t.region,dt[j])
      }
      wbit[is.na(wbit)] <- 0
      wbimodt[is.na(wbimodt)] <- 0
    }
    
    # correction="setcovf"
    if(correction=="setcovf"){
      wsett <- tsetcovf(dt,ndt,bsupt-binft)
      wst <- 1/wsett
    }
    
    Emtout <- .Fortran("emtcoreinh",as.double(ptsx),as.double(ptsy),as.double(ptst),as.integer(npt),
                       as.double(dt),as.integer(ndt),as.double(t.lambda),as.integer(ker2),
                       as.double(ht),as.double(wrt),as.double(wtt),as.double(wbit),
                       as.double(wbimodt),as.double(wst),as.integer(correc2),(eEmt),
                       PACKAGE="msfstpp")
    
    eEmt <- Emtout[[16]]
    Emt0 <- mean(snorm)
    
    invisible(return(list(eEmt=eEmt,Emt0=Emt0,dt=dt,kernel=kernel,Emttheo=Emttheo)))
  }
}