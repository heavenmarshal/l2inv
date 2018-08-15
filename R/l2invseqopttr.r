l2invseqopttr <- function(xi,yi,yobs,nadd,feasible,alpha,
                          func,...,type=c("mvapp","mvei","projei","oei"),
                          mtype=c("zmean","cmean","lmean"),
                          frac=.95,d=NULL,g=0.001,
                          valist=list(),gactl=list(),nthread=4)
{
    xi <- as.matrix(xi)
    type <- match.arg(type)
    mtype <- match.arg(mtype)
    infoname <- paste(type,"info",sep="")
    tlen <- length(yobs)
    infofun <- get(infoname)
    xoptr <- NULL
    valist.default <- list(nmc=500)
    remnames <- setdiff(names(valist.default),names(valist))
    valist <- c(valist,valist.default[remnames])
    valist$yobs <- yobs
    maxinfo <- rep(0,nadd)
    val <- var(yobs)*(tlen-1)
    for(i in 1:nadd)
    {
        cxopt <- l2invoptC(xi,yi,yobs,frac,mtype,d,g,gactl)
        xoptr <- rbind(xoptr,cxopt)
        lbasis <- buildBasis(yi,frac)
        cht <- drop(t(lbasis$basis)%*%yobs/lbasis$redd^2)
        valist$chts2 <- drop(crossprod(yobs-lbasis$basis%*%cht))/tlen
        valist$mindist <- min(apply(lbasis$redd^2*(t(lbasis$coeff)-cht)^2,2,sum))
        valist$barval <-  min(apply((yobs-yi)^2,2,sum))
        py <- svdgpsepms(feasible,xi,yi,frac,mtype=mtype,nthread=nthread)
        info <- infofun(py,alpha,cht,valist)
        mm <- max(info,na.rm=TRUE)
        maxinfo[i] <- mm
        newidx <- which.max(info)
        newx <- feasible[newidx,]
        newy <- func(newx,...)
        xi <- rbind(xi,newx)
        yi <- cbind(yi,newy)
        feasible <- feasible[-newidx,,drop=FALSE]
    }
    xopt <- l2invoptC(xi,yi,yobs,frac,mtype,d,g,gactl)
    xoptr <- rbind(xoptr,xopt)
    ret <- list(xx=xi,yy=yi,xopt=xopt,xoptr=xoptr,maxinfo=maxinfo)
    return(ret)
}
ojsinvseqopttr <- function(xi,yi,yobs,nadd,feasible,mtype=c("zmean","cmean","lmean"),
                           func,...,d=NULL,g=0.001,gactl=list())
{
    xi <- as.matrix(xi)
    lw <- sqrt(apply((yi-yobs)^2,2,sum))
    wmin <- min(lw)
    xoptr <- NULL
    maxinfo <- NULL
    thres <- 0
    for(i in 1:nadd)
    {
        cxopt <- ojsinvopt(xi,yi,yobs,mtype,d=d,g=g,gactl=gactl)
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
    xopt <- ojsinvopt(xi,yi,yobs,mtype,d=d,g=g,gactl=gactl)
    xoptr <- rbind(xoptr,xopt)
    ret <- list(xx=xi,yy=yi,xopt=xopt,xoptr=xoptr,maxinfo=maxinfo)
    return(ret)
}
