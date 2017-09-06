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
mininfo <- function(py,cmin)
{
    ave <- py$mean
    se <- sqrt(py$s2)
    norm <- (cmin-ave)/se
    info <- ifelse(is.infinite(norm),0,
                   se*(norm*pnorm(norm)+dnorm(norm)))
}
