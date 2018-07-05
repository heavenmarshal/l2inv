exactinfo <- function(py,a,b,c,u1,u2)
{
    mu <- py$mean
    s2 <- py$s2
    sigma <- sqrt(s2)
    a1 <- a*s2
    b1 <- b*sigma-mu*sigma*a
    c1 <- c-a*mu^2+2*b*mu
    u1star <- (u1-mu)/sigma
    u2star <- (u2-mu)/sigma
    part1 <- (c1-a1)*(pnorm(u2star)-pnorm(u1star))
    part2 <- a1*(u2star*dnorm(u2star)-u1star*dnorm(u1star))
    part3 <- 2*b1*(dnorm(u2star)-dnorm(u1star))
    info <- part1+part2-part3
}
approxinfo <- function(py,barval,d2,cxi)
{
    n <- length(py$s2)
    p <- 1
    sig2 <- py$s2*d2
    mu2 <- d2*(py$mean-cxi)^2
    barval <- rep(barval,n)
    musum <- mu2+sig2
    mumk <- musum-barval
    upb <- 0.5/sig2
    out <- .C("saddleapprox_R",as.double(sig2),as.double(mu2),
              as.double(barval),as.double(mumk),as.double(upb),
              as.integer(n), as.integer(p), ans=double(n),
              PACKAGE="l2inv")
    info <- out$ans-mumk
}
mcmvinfo <- function(py,yobs,nmc,barval)
{
    basis <- py$basis
    tlen <- nrow(basis)
    numbas <- ncol(basis)
    coeffsd <- sqrt(py$coeffs2)
    sdhat <- sqrt(py$varres)
    nfea <- try(ncol(coeffsd))
    if(inherits(nfea,"try-error"))
        nfea <- length(coeffsd)
    out <- .C("illumcei_R",as.double(py$coeff),as.double(coeffsd),
              as.double(basis), as.double(yobs), as.double(barval),
              as.double(sdhat), as.integer(nmc), as.integer(nfea),
              as.integer(tlen), as.integer(numbas),
              ans = double(nfea), PACKAGE="l2inv")
    return(out$ans)
}
illuExactEI <- function(design,resp,yobs,feasible,
                        mtype=c("zmean","cmean","lmean"),
                        d=NULL,g=0.001)
{
    lbasis <- buildBasis(resp,numbas=1)
    a <- drop(crossprod(lbasis$basis))
    b <- drop(crossprod(lbasis$basis,yobs))
    kmin <- min(apply((resp-yobs)^2,2,sum))
    c <- kmin-drop(crossprod(yobs))
    delta <- b^2+a*c
    if(delta<0) stop("the determinant is smaller than 0!")
    delta <- sqrt(delta)
    u1 <- (b-delta)/a
    u2 <- (b+delta)/a
    coeff <- lbasis$coeff
    gpobj <- if(mtype=="zmean") gpsepms(coeff,design,d,g) else gpseplmms(coeff,design,mtype,d,g)
    pred <- predict(gpobj,feasible)
    delete(gpobj)
    info <- exactinfo(pred,a,b,c,u1,u2)
    return(info)
}

illuApproxEI <- function(design,resp,yobs,feasible,
                         mtype=c("zmean","cmean","lmean"),
                         d=NULL,g=0.001)
{
    lbasis <- buildBasis(resp,numbas=1)
    d2 <- lbasis$redd[1]^2
    basis <- lbasis$basis
    cxi <- drop(crossprod(yobs,basis))/d2
    offset <- drop(crossprod(yobs-basis*cxi))
    kmin <- min(apply((resp-yobs)^2,2,sum))
    barval <- kmin-offset
    coeff <- lbasis$coeff
    gpobj <- if(mtype=="zmean") gpsepms(coeff,design,d,g) else gpseplmms(coeff,design,mtype,d,g)
    pred <- predict(gpobj,feasible)
    delete(gpobj)
    info <- approxinfo(pred,barval,d2,cxi)
    return(info)
}

illuMCMVEI <- function(design,resp,yobs,feasible,nmc,frac=.95,
                       mtype=c("zmean","cmean","lmean"),
                       d=NULL,g=0.001)
{
    py <- svdgpsepms(feasible,design,resp,frac,mtype=mtype,d=d,g=g)
    barval <- min(apply((resp-yobs)^2,2,sum))
    tm <- system.time(info <- mcmvinfo(py,yobs,nmc,barval))
    ret <- list(info=info,time=tm[3])
    return(ret)
}
illuApproxMVEI <- function(design,resp,yobs,feasible,frac=.95,
                           mtype=c("zmean","cmean","lmean"),
                           d=NULL,g=0.001)
{
    py <- svdgpsepms(feasible,design,resp,frac,mtype=mtype,d=d,g=g)
    barval <- min(apply((resp-yobs)^2,2,sum))
    valist <- list(yobs=yobs,barval=barval)
    tm <- system.time(info <- oeiinfo(py,NULL,NULL,valist))
    ret <- list(info=info,time=tm[3])
    return(ret)
}
