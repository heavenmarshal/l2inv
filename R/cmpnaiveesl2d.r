cmpnaiveesl2d <- function(xi,yi,yobs,nadd,feasible,grid,func,...,
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
    xoptesl2d <- esl2dinv(xi, yi, yobs, grid,frac, d=d, g=g)
    xoptnaive <- naiveinv(xi, yi, yobs, grid, frac, d=d, g=g)
    ret <- list(xx=xi,yy=yi,xoptesl2d=xoptesl2d,
                xoptnaive=xoptnaive,maxinfo=maxinfo)
    return(ret)
}
cmpnaiveesl2dopt <- function(xi,yi,yobs,nadd,feasible,func,...,
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
    xoptesl2d <- esl2dinvoptR(xi,yi,yobs,frac,critype="esl2d",d=d,g=g)
    xoptnaive <- esl2dinvoptR(xi,yi,yobs,frac,critype="naive",d=d,g=g)
    ret <- list(xx=xi,yy=yi,xoptesl2d=xoptesl2d,
                xoptnaive=xoptnaive,maxinfo=maxinfo)
    return(ret)
}
