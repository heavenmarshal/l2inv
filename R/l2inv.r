l2inv <- function(design,resp,yobs,feasible,frac=.95,
                  mtype=c("zmean","cmean","lmean"),
                  maxiter=50,tol=1e-3,d=NULL,g=0.001)
{
    mtype <- match.arg(mtype)
    design <- as.matrix(design)
    lbasis <- buildBasis(resp,frac)
    if(any(is.na(yobs)))
    {
        yimp <- imputels(yobs,lbasis,maxiter,tol)
        cht <- yimp$cht
    }
    else
        cht <- drop(t(lbasis$basis)%*%yobs/lbasis$redd^2)

    numbas <- lbasis$numbas
    coeff <- lbasis$coeff
    reddsq <- lbasis$redd^2
    nfea <- nrow(feasible)
    criter <- rep(0,nfea)
    for(i in 1:numbas)
    {
        ccoeff <- coeff[,i]
        gpobj <- if(mtype=="zmean") gpsepms(ccoeff,design,d,g) else gpseplmms(ccoeff,design,mtype,d,g)
        pred <- predict(gpobj,feasible)
        delete(gpobj)
        criter <- criter+reddsq[i]*(pred$s2+(pred$mean-cht[i])^2)
    }
    xopt <- feasible[which.min(criter),]
    return(xopt)
}
ojsinv <- function(design,resp,yobs,feasible,mtype=c("zmean","cmean","lmean"),
                   q=.5,d=NULL,g=0.001)
{
    mtype <- match.arg(mtype)
    design <- as.matrix(design)
    isnmis <- !is.na(yobs)
    yobsred <- yobs[isnmis]
    respred <- resp[isnmis,]
    dev <- sqrt(apply((respred-yobsred)^2,2,sum))
    gpobj <- if(mtype=="zmean") gpsepms(dev,design,d,g) else gpseplmms(dev,design,mtype,d,g)
    py <- predict(gpobj,feasible)
    criter <- py$mean + qnorm(q)*sqrt(py$s2)
    delete(gpobj)
    xopt <- feasible[which.min(criter),]
    return(xopt)
}
evallik_zm <- function(resp,design,mtype,d,g)
{
    gpobj <- gpsepms(resp,design,d,g)
    llik <- getlik(gpobj)
    delete(gpobj)
    return(llik)
}
evallik_lm <- function(resp,design,mtype,d,g)
{
    gpobj <- gpseplmms(resp,design,mtype,d,g)
    llik <- getlik(gpobj)
    delete(gpobj)
    return(llik)
}
lrinv <- function(design,resp,yobs,timepoints,feasible,
                  mtype=c("zmean","cmean","lmean"),
                  q=0.5,d=NULL,g=0.001,gl=0.1,nthread=4)
{
    mtype <- match.arg(mtype)
    tlen <- length(yobs)
    delta <- resp-yobs
    conlik <- -0.5*tlen*log(colMeans(delta^2))
    tinput <- as.matrix(timepoints)
    suffix <- if(mtype=="zmean")"zm" else "lm"
    evallik <- get(paste("evallik",suffix,sep="_"))

    cl <- parallel::makeCluster(nthread)
    uclik <- tryCatch(parallel::parApply(cl,delta,2,evallik,tinput,mtype,d,gl),
                      finally=parallel::stopCluster(cl))
    likratio <- -2*(conlik-uclik)
    gpobj <- if(mtype=="zmean") gpsepms(likratio,design,d,g) else gpseplmms(likratio,design,mtype,d,g)
    py <- predict(gpobj,feasible)
    delete(gpobj)
    criter <- py$mean + qnorm(q)*sqrt(py$s2)
    xopt <- feasible[which.min(criter),]
    return(xopt)
}
