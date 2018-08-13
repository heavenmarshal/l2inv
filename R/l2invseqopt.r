l2invseqopt <- function(xi,yi,yobs,nadd,feasible,alpha,func,...,
                        type=c("mvapp","mvei","projei","oei"),
                        mtype=c("zmean","cmean","lmean"),
                        relthres=0,frac=.95,d=NULL,g=0.001,
                        valist=list(),nthread=4)
{
    xi <- as.matrix(xi)
    type <- match.arg(type)
    mtype <- match.arg(mtype)
    infoname <- paste(type,"info",sep="")
    tlen <- length(yobs)
    infofun <- get(infoname)
    valist.default <- list(nmc=500)
    remnames <- setdiff(names(valist.default),names(valist))
    valist <- c(valist,valist.default[remnames])
    valist$yobs <- yobs
    maxinfo <- rep(0,nadd)
    thres <- 0
    for(i in 1:nadd)
    {
        lbasis <- buildBasis(yi,frac)
        cht <- drop(t(lbasis$basis)%*%yobs/lbasis$redd^2)
        valist$chts2 <- drop(crossprod(yobs-lbasis$basis%*%cht))/tlen
        valist$mindist <- min(apply(lbasis$redd^2*(t(lbasis$coeff)-cht)^2,2,sum))
        valist$barval <-  min(apply((yobs-yi)^2,2,sum))
        py <- svdgpsepms(feasible,xi,yi,frac,mtype=mtype,nthread=nthread)
        info <- infofun(py,alpha,cht,valist)
        mm <- max(info,na.rm=TRUE)
        maxinfo[i] <- mm
        if(i==1) thres <- mm*relthres
        if(mm<thres) break
        newidx <- which.max(info)
        newx <- feasible[newidx,]
        newy <- func(newx,...)
        xi <- rbind(xi,newx)
        yi <- cbind(yi,newy)
        feasible <- feasible[-newidx,,drop=FALSE]
    }
    xopt <- l2invopt(xi,yi,yobs,frac,d=d,g=g)
    ret <- list(xx=xi,yy=yi,xopt=xopt,maxinfo=maxinfo)
    return(ret)
}
ojsinvseqopt <- function(xi,yi,yobs,nadd,feasible,mtype=c("zmean","cmean","lmean"),
                         func,...,relthres=0,d=NULL,g=0.001)
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
        if(i == 1) thres <- mm*relthres
        if(mm<thres) break
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
    xopt <- ojsinvopt(xi,yi,yobs,mtype,d=d,g=g)
    ret <- list(xx=xi,yy=yi,xopt=xopt,maxinfo=maxinfo)
    return(ret)
}

