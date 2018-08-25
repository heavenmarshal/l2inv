esl2dinvseq <- function(xi,yi,yobs,nadd,feasible,grid,func,...,
                        mtype=c("zmean","cmean","lmean"),
                        frac=.95,d=NULL,g=0.001,
                        nthread=4)
{
    xi <- as.matrix(xi)
    mtype <- match.arg(mtype)
    tlen <- length(yobs)
    maxinfo <- rep(0,nadd)
    for(i in 1:nadd)
    {
        lbasis <- buildBasis(yi,frac)
        barval <-  min(apply((yobs-yi)^2,2,sum))
        py <- svdgpsepms(feasible,xi,yi,frac,mtype=mtype,nthread=nthread)
        info <- oeiinfo(py,yobs,barval)
        mm <- max(info,na.rm=TRUE)
        maxinfo[i] <- mm
        newidx <- which.max(info)
        newx <- feasible[newidx,]
        newy <- func(newx,...)
        xi <- rbind(xi,newx)
        yi <- cbind(yi,newy)
        feasible <- feasible[-newidx,,drop=FALSE]
    }
    xopt <- esl2dinv(xi,yi,yobs,grid,frac,d=d,g=g)
    ret <- list(xx=xi,yy=yi,xopt=xopt,maxinfo=maxinfo)
    return(ret)
}
sl2invseq <- function(xi,yi,yobs,nadd,feasible,grid,mtype=c("zmean","cmean","lmean"),
                      func,...,d=NULL,g=0.001)
{
    xi <- as.matrix(xi)
    lw <- sqrt(apply((yi-yobs)^2,2,sum))
    wmin <- min(lw)
    maxinfo <- NULL
    for(i in 1:nadd)
    {
        gpobj <- if(mtype=="zmean") gpsepms(lw,xi,d,g) else gpseplmms(lw,xi,mtype,d,g)
        py <- predict(gpobj,feasible)
        delete(gpobj)
        info <- mininfo(py,wmin)
        mm <- max(info,na.rm=TRUE)
        maxinfo <- c(maxinfo,mm)
        newidx <- which.max(info)
        newx <- feasible[newidx,]
        newy <- func(newx,...)
        xi <- rbind(xi,newx)
        yi <- cbind(yi,newy)
        feasible <- feasible[-newidx,,drop=FALSE]
        newlw <- sqrt(sum((newy-yobs)^2))
        wmin <- min(wmin,newlw)
        lw <- c(lw,newlw)
    }
    xopt <- sl2inv(xi,yi,yobs,grid,mtype,d=d,g=g)
    ret <- list(xx=xi,yy=yi,xopt=xopt,maxinfo=maxinfo)
    return(ret)
}
