## trace the update for the sequential design, for experiment not application.
esl2dinvseqtr <- function(xi,yi,yobs,nadd,feasible,grid, func,...,
                          mtype=c("zmean","cmean","lmean"),
                          frac=.95,d=NULL,g=0.001,nthread=4)
{
    xi <- as.matrix(xi)
    nd <- ncol(xi)
    mtype <- match.arg(mtype)
    tlen <- length(yobs)
    nfea <- nrow(feasible)
    xoptr <- NULL
    maxinfo <- NULL
    for(i in 1:nadd)
    {
        tfea <- list(feasible,grid)
        lbasis <- buildBasis(yi,frac)
        cht <- drop(t(lbasis$basis)%*%yobs/lbasis$redd^2)
        barval <-  min(apply((yobs-yi)^2,2,sum))
        py <- svdgpsepmslist(tfea,xi,yi,frac,mtype=mtype,nthread=nthread)
        pyfea <- py[[1]]
        pygrid <- py[[2]]
        numbas <- lbasis$numbas
        info <- oeiinfo(pyfea,yobs,barval)
        mm <- max(info,na.rm=TRUE)
        maxinfo <- c(maxinfo,mm)
        if(numbas==1)
            criter <- pygrid$d2*(pygrid$coeffs2+(pygrid$coeff-cht)^2)
        else
            criter <- apply(pygrid$d2*(pygrid$coeffs2+(pygrid$coeff-cht)^2),2,sum)
        cxopt <- grid[which.min(criter),]
        xoptr <- rbind(xoptr,cxopt)
        newidx <- which.max(info)
        newx <- feasible[newidx,]
        newy <- func(newx,...)
        xi <- rbind(xi,newx)
        yi <- cbind(yi,newy)
        feasible <- feasible[-newidx,,drop=FALSE]
        nfea <- nfea-1
    }
    xopt <- esl2dinv(xi,yi,yobs,grid,frac,d=d,g=g)
    xoptr <- rbind(xoptr,xopt)
    ret <- list(xx=xi,yy=yi,xopt=xopt,xoptr=xoptr,maxinfo=maxinfo)
    return(ret)
}
sl2invseqtr <- function(xi,yi,yobs,nadd,feasible,grid,
                        mtype=c("zmean","cmean","lmean"),
                        func,...,d=NULL,g=0.001)
{
    xi <- as.matrix(xi)
    nd <- ncol(xi)
    lw <- sqrt(apply((yi-yobs)^2,2,sum))
    tlen <- length(yobs)
    wmin <- min(lw)
    nfea <- nrow(feasible)
    xoptr <- NULL
    maxinfo <- NULL
    for(i in 1:nadd)
    {
        tfea <- rbind(feasible,grid)
        gpobj <- if(mtype=="zmean") gpsepms(lw,xi,d,g) else gpseplmms(lw,xi,mtype,d,g)
        py <- predict(gpobj,tfea)
        delete(gpobj)
        pyfea <- list(mean=py$mean[1:nfea], s2=py$s2[1:nfea])
        info <- mininfo(pyfea,wmin)
        mm <- max(info,na.rm=TRUE)
        maxinfo <- c(maxinfo,mm)
        lwhat <- py$mean[-(1:nfea)]
        cxopt <- grid[which.min(lwhat),]
        xoptr <- rbind(xoptr,cxopt)
        newidx <- which.max(info)
        newx <- feasible[newidx,]
        newy <- func(newx,...)
        xi <- rbind(xi,newx)
        yi <- cbind(yi,newy)

        feasible <- feasible[-newidx,,drop=FALSE]
        newlw <- sqrt(sum((newy-yobs)^2))
        wmin <- min(wmin,newlw)
        lw <- c(lw,newlw)
        nfea <- nfea-1
    }
    xopt <- sl2inv(xi,yi,yobs,grid,mtype,d=d,g=g)
    xoptr <- rbind(xoptr,xopt)
    ret <- list(xx=xi,yy=yi,xopt=xopt,xoptr=xoptr,maxinfo=maxinfo)
    return(ret)
}
