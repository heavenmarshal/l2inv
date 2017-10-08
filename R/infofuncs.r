naiveinfo <- function(py,alpha,cht,valist)
{
    infomat <- alpha*py$coeffs2-(py$coeff-cht)^2
    info <- infomat*py$d2
    info <- apply(info,2,sum)
}
mvconinfo <- function(py,alpha,cht,valist)
{
    barrier <- alpha * apply(py$coeffs2*py$d2,2,sum)
    numbas <- nrow(py$coeff)
    nfea <- ncol(py$coeff)
    out <- .C("mvcon_R", as.double(py$coeff), as.double(sqrt(py$coeffs2)),
              as.double(barrier), as.double(cht), as.double(py$d2),
              as.integer(nfea), as.integer(valist$nmc), as.integer(numbas),
              ans=double(nfea), PACKAGE="l2inv")
    info <- out$ans
}
mvappinfo <- function(py,alpha,cht,valist)
{
    sig2 <- py$coeffs2*py$d2
    n <- ncol(sig2)
    p <- nrow(sig2)
    mu2 <- py$d2*(py$coeff-cht)^2
    barval <- alpha*apply(sig2,2,sum)
    musum <- apply(mu2+sig2,2,sum)
    mumk <- musum-barval
    upb <- apply(0.5/sig2,2,min)
    out <- .C("saddleapprox_R",as.double(sig2), as.double(mu2),
              as.double(barval), as.double(mumk), as.double(upb),
              as.integer(n), as.integer(p), ans=double(n))
    info <- out$ans-mumk
    return(info)
}
## alpha is not useful
mveiinfo <- function(py,alpha,cht,valist)
{
    sig2 <- py$coeffs2*py$d2
    n <- ncol(sig2)
    p <- nrow(sig2)
    mu2 <- py$d2*(py$coeff-cht)^2
    barval <- rep(valist$mindist,n)
    musum <- apply(mu2+sig2,2,sum)
    mumk <- musum-barval
    upb <- apply(0.5/sig2,2,min)
    out <- .C("saddleapprox_R",as.double(sig2), as.double(mu2),
              as.double(barval), as.double(mumk), as.double(upb),
              as.integer(n), as.integer(p), ans=double(n))
    info <- out$ans-mumk
    return(info)
}
oeiinfo <- function(py,alpha,cht,valist)
{
    yobs <- valist$yobs
    d2 <- py$d2
    L <- length(yobs)
    p <- length(d2)
    n <- ncol(py$coeff)
    bz <- drop(t(py$basis)%*%yobs)
    z2 <- drop(crossprod(yobs))
    d2c <- py$coeff*d2
    scalec <- apply(py$coeff*d2c,2,sum)
    bzc <- apply(bz*py$coeff,2,sum)
    scales2 <- d2*py$coeffs2
    amat <- scales2+py$varres     #input
    sand <- py$coeffs2/amat
    sbz <- sand*bz
    sd2c <- sand*d2c
    sbz2 <- apply(sbz*bz,2,sum)
    sd2c2 <- apply(sd2c*d2c,2,sum)
    sbzd2c <- apply(sand*d2c*bz,2,sum)
    iomemu2 <- z2+scalec-2*bzc
    iomemu2 <- iomemu2-sbz2-sd2c2+2*sbzd2c
    iomemu2 <- iomemu2/py$varres    #input
    mubstar <- sbz-sd2c
    mub2star <- (bz-d2c)*mubstar        #input
    bound <- apply(amat,2,max)
    bound <- 0.5/bound
    mumk <- z2+scalec-2*bzc+apply(scales2,2,sum)+py$varres*L
    mumk <- mumk - valist$barval
    ret <- .C("oeiinfo_R", as.integer(n), as.integer(p),
              as.integer(L), as.double(py$varres),
              as.double(valist$barval),as.double(iomemu2),
              as.double(bound), as.double(amat),
              as.double(mub2star), as.double(mumk),
              info = double(n), PACKAGE="l2inv")
    info <- ret$info-mumk
    return(info)
}
mininfo <- function(py,cmin)
{
    ave <- py$mean
    se <- sqrt(py$s2)
    norm <- (cmin-ave)/se
    info <- ifelse(is.infinite(norm),0,
                   se*(norm*pnorm(norm)+dnorm(norm)))
}
## coneinfo <- function(py,alpha,cht,valist)
## {
##     mse <- py$coeffs2
##     se <- sqrt(mse)
##     norm <- (cht-py$coeff)/se
##     u2 <- norm+alpha
##     u1 <- norm-alpha
##     du2 <- dnorm(u2)
##     du1 <- dnorm(u1)
##     part1 <- mse*(alpha*alpha-norm*norm-1)*(pnorm(u2)-pnorm(u1))
##     part2 <- mse*(u2*du2-u1*du1)
##     part3 <- -2*norm*mse*(du2-du1)
##     info <- part1+part2+part3
##     info <- ifelse(is.na(info) | is.infinite(info),0,info)
##     info <- apply(info*py$d2,2,sum)
## }
## mvextinfo <- function(py,alpha,cht,valist)
## {
##     chts2 <- valist$chts2
##     sig2 <- py$coeffs2*py$d2
##     n <- ncol(sig2)
##     p <- nrow(sig2)
##     mu2 <- py$d2*(py$coeff-cht)^2
##     barval <- alpha*apply(sig2,2,sum)
##     musum <- apply(mu2+sig2,2,sum)+p*chts2
##     mumk <- musum-barval
##     s2ext <- sig2+chts2
##     upb <- apply(0.5/s2ext,2,min)
##     out <- .C("saddleapprox_R",as.double(s2ext), as.double(mu2),
##               as.double(barval), as.double(mumk), as.double(upb),
##               as.integer(n), as.integer(p), ans=double(n))
##     info <- out$ans-mumk
##     return(info)
## }
