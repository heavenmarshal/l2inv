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
mcmvinfo <- function(py,yobs,nmc,barval,ompthread)
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
              as.integer(tlen), as.integer(numbas), as.integer(ompthread),
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
                       d=NULL,g=0.001,ompthread=4)
{
    py <- svdgpsepms(feasible,design,resp,frac,mtype=mtype,d=d,g=g)
    barval <- min(apply((resp-yobs)^2,2,sum))
    info <- mcmvinfo(py,yobs,nmc,barval,ompthread)
    return(info)
}
illuApproxMVEI <- function(design,resp,yobs,feasible,frac=.95,
                           mtype=c("zmean","cmean","lmean"),
                           d=NULL,g=0.001)
{
    py <- svdgpsepms(feasible,design,resp,frac,mtype=mtype,d=d,g=g)
    barval <- min(apply((resp-yobs)^2,2,sum))
    valist <- list(yobs=yobs,barval=barval)
    info <- oeiinfo(py,NULL,NULL,valist)
    return(info)
}

illuoeiFixVarres <- function(coeff,coeffs2,varres,py,yobs,barval)
{
    d2 <- py$d2
    p <- length(d2)
    n <- ncol(coeff)
    if(nrow(coeff) != p || nrow(coeffs2) != p || ncol(coeffs2) !=n)
        stop("incorrect dimension of inputs!")
    L <- length(yobs)
    bz <- drop(t(py$basis)%*%yobs)
    z2 <- drop(crossprod(yobs))
    d2c <- coeff*d2
    scalec <- apply(coeff*d2c,2,sum)
    bzc <- apply(bz*coeff,2,sum)
    scales2 <- d2*coeffs2
    amat <- scales2+varres     #input
    sand <- coeffs2/amat
    sbz <- sand*bz
    sd2c <- sand*d2c
    sbz2 <- apply(sbz*bz,2,sum)
    sd2c2 <- apply(sd2c*d2c,2,sum)
    sbzd2c <- apply(sand*d2c*bz,2,sum)
    iomemu2 <- z2+scalec-2*bzc
    iomemu2 <- iomemu2-sbz2-sd2c2+2*sbzd2c
    iomemu2 <- iomemu2/varres    #input
    mubstar <- sbz-sd2c
    mub2star <- (bz-d2c)*mubstar        #input
    bound <- apply(amat,2,max)
    bound <- 0.5/bound
    mumk <- z2+scalec-2*bzc+apply(scales2,2,sum)+varres*L
    mumk <- mumk - barval
    ret <- .C("illuoeiFixVarres_R", as.integer(n), as.integer(p),
              as.integer(L), as.double(varres),
              as.double(barval),as.double(iomemu2),
              as.double(bound), as.double(amat),
              as.double(mub2star), as.double(mumk),
              info = double(n), spoints = double(n),
              lambda3 = double(n), wval = double(n),
              qval = double(n), dk2val = double(n),
              PACKAGE="l2inv")
    info <- ret$info-mumk
    out <- list(info=info,spoints=ret$spoints,lambda3=ret$lambda3,
                wval=ret$wval,qval=ret$qval,dk2val=ret$dk2val,
                mumk=mumk)
    return(out)
}

illuoeiFixCoef <- function(coeff,coeffs2,varres,py,yobs,barval)
{
    d2 <- py$d2
    p <- length(d2)
    if(length(coeff) != p || length(coeffs2) != p)
        stop("incorrect dimension of inputs!")
    L <- length(yobs)
    n <- length(varres)
    bz <- drop(t(py$basis)%*%yobs)
    z2 <- drop(crossprod(yobs))
    d2c <- coeff*d2
    scalec <- sum(coeff*d2c)## apply(coeff*d2c,2,sum)
    bzc <- sum(bz*coeff) ## apply(bz*coeff,2,sum)
    scales2 <- d2*coeffs2
    amat <- outer(scales2,varres,"+")     #input p x n matrix
    sand <- coeffs2/amat
    sbz <- sand*bz
    sd2c <- sand*d2c
    sbz2 <- apply(sbz*bz,2,sum)
    sd2c2 <- apply(sd2c*d2c,2,sum)
    sbzd2c <- apply(sand*d2c*bz,2,sum)
    iomemu2 <- z2+scalec-2*bzc
    iomemu2 <- iomemu2-sbz2-sd2c2+2*sbzd2c
    iomemu2 <- iomemu2/varres    #input length n vector
    mubstar <- sbz-sd2c
    mub2star <- (bz-d2c)*mubstar        #input p x n matrix
    bound <- apply(amat,2,max)
    bound <- 0.5/bound
    mumk <- z2+scalec-2*bzc+sum(scales2)+varres*L
    mumk <- mumk - barval
    ret <- .C("illuoeiMultVarres_R", as.integer(n), as.integer(p),
              as.integer(L), as.double(varres),
              as.double(barval),as.double(iomemu2),
              as.double(bound), as.double(amat),
              as.double(mub2star), as.double(mumk),
              info = double(n), spoints = double(n),
              lambda3 = double(n), wval = double(n),
              qval = double(n), dk2val = double(n),
              PACKAGE="l2inv")
    info <- ret$info-mumk
    out <- list(info=info,spoints=ret$spoints,lambda3=ret$lambda3,
                wval=ret$wval,qval=ret$qval,dk2val=ret$dk2val,
                mumk=mumk)

    return(out)
}

illuoei <- function(coeff,coeffs2,varres,py,yobs,barval)
{
    d2 <- py$d2
    p <- length(d2)
    n <- ncol(coeff)
    if(nrow(coeff) != p || nrow(coeffs2) != p || ncol(coeffs2) !=n || length(varres) != n)
        stop("incorrect dimension of inputs!")
    L <- length(yobs)
    bz <- drop(t(py$basis)%*%yobs)
    z2 <- drop(crossprod(yobs))
    d2c <- coeff*d2
    scalec <- apply(coeff*d2c,2,sum)
    bzc <- apply(bz*coeff,2,sum)
    scales2 <- d2*coeffs2
    amat <- t(t(scales2)+varres)     #input
    sand <- coeffs2/amat
    sbz <- sand*bz
    sd2c <- sand*d2c
    sbz2 <- apply(sbz*bz,2,sum)
    sd2c2 <- apply(sd2c*d2c,2,sum)
    sbzd2c <- apply(sand*d2c*bz,2,sum)
    iomemu2 <- z2+scalec-2*bzc
    iomemu2 <- iomemu2-sbz2-sd2c2+2*sbzd2c
    iomemu2 <- iomemu2/varres    #input
    mubstar <- sbz-sd2c
    mub2star <- (bz-d2c)*mubstar        #input
    bound <- apply(amat,2,max)
    bound <- 0.5/bound
    mumk <- z2+scalec-2*bzc+apply(scales2,2,sum)+varres*L
    mumk <- mumk - barval
    ret <- .C("illuoeiMultVarres_R", as.integer(n), as.integer(p),
              as.integer(L), as.double(varres),
              as.double(barval),as.double(iomemu2),
              as.double(bound), as.double(amat),
              as.double(mub2star), as.double(mumk),
              info = double(n), spoints = double(n),
              lambda3 = double(n), wval = double(n),
              qval = double(n), dk2val = double(n),
              PACKAGE="l2inv")
    info <- ret$info-mumk
    info <- ret$info-mumk
    out <- list(info=info,spoints=ret$spoints,lambda3=ret$lambda3,
                wval=ret$wval,qval=ret$qval,dk2val=ret$dk2val,
                mumk=mumk)
    return(out)
}
