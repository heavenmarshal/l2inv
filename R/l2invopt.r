esl2d <- function(x,gpobjs,lbasis,cht)
{
    x <- matrix(x,nrow=1)
    nbas <- lbasis$numbas
    pmean <- ps2 <- vector("double",nbas)
    for(i in 1:nbas)
    {
        pred <- predict(gpobjs[[i]],x)
        pmean[i] <- pred$mean
        ps2[i] <- pred$s2
    }
    reddsq <- lbasis$redd^2
    crit <- sum(reddsq*(ps2+(pmean-cht)^2))
}
naive <- function(x,gpobjs,lbasis,cht)
{
    x <- matrix(x,nrow=1)
    nbas <- lbasis$numbas
    pmean <- ps2 <- vector("double",nbas)
    for(i in 1:nbas)
    {
        pred <- predict(gpobjs[[i]],x)
        pmean[i] <- pred$mean
        ps2[i] <- pred$s2
    }
    reddsq <- lbasis$redd^2
    crit <- sum(reddsq*(pmean-cht)^2)
}
l2invopt <- function(design,resp,yobs,frac=.95,
                     mtype=c("zmean","cmean","lmean"),
                     critype=c("esl2d","naive"),
                     d=NULL,g=0.001)
{
    ndim <- ncol(design)
    mtype <- match.arg(mtype)
    critype <- match.arg(critype)
    design <- as.matrix(design)
    lbasis <- buildBasis(resp,frac)
    cht <- drop(t(lbasis$basis)%*%yobs/lbasis$redd^2)
    numbas <- lbasis$numbas
    critfun <- get(critype)
    coeff <- lbasis$coeff
    gpobjs <- vector("list",numbas)
    for(i in 1:numbas)
    {
        ccoeff <- coeff[,i]
        gpobjs[[i]] <- if(mtype=="zmean") gpsepms(ccoeff,design,d,g)
                       else gpseplmms(ccoeff,design,mtype,d,g)
    }
    opt <- genoud(critfun,ndim,Domains=cbind(rep(0,ndim),rep(1,ndim)),
                  boundary.enforcement=2,print.level=0,
                  gpobjs=gpobjs, lbasis=lbasis, cht=cht)
    xopt <- opt$par
    for(i in 1:numbas) delete(gpobjs[[i]])
    return(xopt)
}
