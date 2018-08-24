esl2dinvseqopt <- function(xi,yi,yobs,nadd,feasible,func,...,
                           mtype=c("zmean","cmean","lmean"),
                           frac=.95,d=NULL,g=0.001,gactl=list(),
                           nthread=4)
{
    xi <- as.matrix(xi)
    mtype <- match.arg(mtype)
    tlen <- length(yobs)
    maxinfo <- rep(0,nadd)
    for(i in 1:nadd)
    {
        lbasis <- buildBasis(yi,frac)
        cht <- drop(t(lbasis$basis)%*%yobs/lbasis$redd^2)
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
    xopt <- esl2dinvoptC(xi,yi,yobs,frac,d=d,g=g,gactl=gactl)
    ret <- list(xx=xi,yy=yi,xopt=xopt,maxinfo=maxinfo)
    return(ret)
}
sl2invseqopt <- function(xi,yi,yobs,nadd,feasible,mtype=c("zmean","cmean","lmean"),
                         func,...,relthres=0,d=NULL,g=0.001,gactl=list())
{
    xi <- as.matrix(xi)
    lw <- sqrt(apply((yi-yobs)^2,2,sum))
    wmin <- min(lw)
    maxinfo <- NULL
    thres <- 0
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
    xopt <- sl2invopt(xi,yi,yobs,mtype,d=d,g=g,gactl=gactl)
    ret <- list(xx=xi,yy=yi,xopt=xopt,maxinfo=maxinfo)
    return(ret)
}

