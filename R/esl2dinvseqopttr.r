esl2dinvseqopttr <- function(xi,yi,yobs,nadd,feasible, func,...,
                             mtype=c("zmean","cmean","lmean"),
                             frac=.95,d=NULL,g=0.001,
                             gactl=list(),nthread=4)
{
    xi <- as.matrix(xi)
    mtype <- match.arg(mtype)
    tlen <- length(yobs)
    xoptr <- NULL
    maxinfo <- rep(0,nadd)
    val <- var(yobs)*(tlen-1)
    for(i in 1:nadd)
    {
        cxopt <- esl2dinvopt(xi,yi,yobs,frac,mtype,d,g,gactl)
        xoptr <- rbind(xoptr,cxopt)
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
    xopt <- esl2dinvopt(xi,yi,yobs,frac,mtype,d,g,gactl)
    xoptr <- rbind(xoptr,xopt)
    ret <- list(xx=xi,yy=yi,xopt=xopt,xoptr=xoptr,maxinfo=maxinfo)
    return(ret)
}
sl2invseqopttr <- function(xi,yi,yobs,nadd,feasible,mtype=c("zmean","cmean","lmean"),
                           func,...,d=NULL,g=0.001,gactl=list())
{
    xi <- as.matrix(xi)
    lw <- sqrt(apply((yi-yobs)^2,2,sum))
    wmin <- min(lw)
    xoptr <- NULL
    maxinfo <- NULL
    for(i in 1:nadd)
    {
        cxopt <- sl2invopt(xi,yi,yobs,mtype,d=d,g=g,gactl=gactl)
        xoptr <- rbind(xoptr,cxopt)
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
    xoptr <- rbind(xoptr,xopt)
    ret <- list(xx=xi,yy=yi,xopt=xopt,xoptr=xoptr,maxinfo=maxinfo)
    return(ret)
}
