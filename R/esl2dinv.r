esl2dinv <- function(design,resp,yobs,feasible,frac=.95,
                     mtype=c("zmean","cmean","lmean"),
                     d=NULL,g=0.001)
{
    mtype <- match.arg(mtype)
    design <- as.matrix(design)
    lbasis <- buildBasis(resp,frac)
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
naiveinv <- function(design,resp,yobs,feasible,frac=.95,
                     mtype=c("zmean","cmean","lmean"),
                     d=NULL,g=0.001)
{
    mtype <- match.arg(mtype)
    design <- as.matrix(design)
    lbasis <- buildBasis(resp,frac)
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
        criter <- criter+reddsq[i]*(pred$mean-cht[i])^2
    }
    xopt <- feasible[which.min(criter),]
    return(xopt)
}

sl2inv <- function(design,resp,yobs,feasible,mtype=c("zmean","cmean","lmean"),
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
