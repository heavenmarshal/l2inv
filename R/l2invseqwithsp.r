l2invseqwithsp <- function(xi,yi,yobs,nadd,feasible,grid,alpha,func,...,
                           mtype=c("zmean","cmean","lmean"),
                           relthres=0,frac=.95,d=NULL,g=0.001,
                           valist=list(),nthread=4)
{
    xi <- as.matrix(xi)
    mtype <- match.arg(mtype)
    tlen <- length(yobs)
    valist.default <- list(nmc=500)
    remnames <- setdiff(names(valist.default),names(valist))
    valist <- c(valist,valist.default[remnames])
    valist$yobs <- yobs
    maxinfo <- rep(0,nadd)
    thres <- 0
    saddlepoints <- list()
    for(i in 1:nadd)
    {
        lbasis <- buildBasis(yi,frac)
        cht <- drop(t(lbasis$basis)%*%yobs/lbasis$redd^2)
        valist$chts2 <- drop(crossprod(yobs-lbasis$basis%*%cht))/tlen
        valist$mindist <- min(apply(lbasis$redd^2*(t(lbasis$coeff)-cht)^2,2,sum))
        valist$barval <-  min(apply((yobs-yi)^2,2,sum))
        py <- svdgpsepms(feasible,xi,yi,frac,mtype=mtype,nthread=nthread)
        rlist <- oeiinfowithsp(py,alpha,cht,valist)
        info <- rlist$info
        saddlepoints[[i]] <- rlist$saddlepoints
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
    xopt <- l2inv(xi,yi,yobs,grid,frac,d=d,g=g)
    ret <- list(xx=xi,yy=yi,xopt=xopt,maxinfo=maxinfo,saddlepoints=saddlepoints)
    return(ret)
}
