## trace the update for the sequential design, for experiment not application.
l2invseqtr <- function(xi,yi,yobs,nadd,feasible,grid,alpha,
                       func,...,type=c("mvapp","mvei","projei","oei"),
                       mtype=c("zmean","cmean","lmean"),relthres=0,
                       difthres=0,difstep=1,dxithres=0,dxistep=1,
                       frac=.95,d=NULL,g=0.001,
                       valist=list(nmc=500),nthread=4)
{
    xi <- as.matrix(xi)
    nd <- ncol(xi)
    type <- match.arg(type)
    mtype <- match.arg(mtype)
    infoname <- paste(type,"info",sep="")
    tlen <- length(yobs)
    infofun <- get(infoname)
    nfea <- nrow(feasible)
    xoptr <- NULL
    valist.default <- list(nmc=500)
    remnames <- setdiff(names(valist.default),names(valist))
    valist <- c(valist,valist.default[remnames])
    valist$yobs <- yobs
    maxinfo <- NULL
    thres <- 0
    difqueue <- initQueue(difstep)
    dxiqueue <- initQueue(dxistep)
    stopflag <- 0
    val <- var(yobs)*(tlen-1)
    for(i in 1:nadd)
    {
        tfea <- list(feasible,grid)
        lbasis <- buildBasis(yi,frac)
        cht <- drop(t(lbasis$basis)%*%yobs/lbasis$redd^2)
        valist$chts2 <- drop(crossprod(yobs-lbasis$basis%*%cht))/tlen
        valist$mindist <- min(apply(lbasis$redd^2*(t(lbasis$coeff)-cht)^2,2,sum))
        valist$barval <-  min(apply((yobs-yi)^2,2,sum))
        py <- svdgpsepmslist(tfea,xi,yi,frac,mtype=mtype,nthread=nthread)
        pyfea <- py[[1]]
        pygrid <- py[[2]]
        numbas <- lbasis$numbas
        info <- infofun(pyfea,alpha,cht,valist)
        mm <- max(info,na.rm=TRUE)
        maxinfo <- c(maxinfo,mm)
        if(i == 1) thres <- relthres*mm
        if(mm<thres)
        {
            stopflag <- 1
            break
        }
        difratio <- evalRatio(difqueue,mm)
        if(isFull(difqueue) && difratio < difthres)
        {
            stopflag <- 2
            break
        }
        difqueue <- enQueue(difqueue,mm)
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
        ## decide if early stop should be applied
        cdev <- sum((yobs-func(cxopt,...))^2)/val
        dxiratio <- evalRatio(dxiqueue,cdev)
        if(isFull(dxiqueue) && dxiratio < dxithres)
        {
            stopflag <- 3
            break
        }
        dxiqueue <- enQueue(dxiqueue,cdev)
        feasible <- feasible[-newidx,,drop=FALSE]
        nfea <- nfea-1
    }
    xopt <- l2inv(xi,yi,yobs,grid,frac,d=d,g=g)
    xoptr <- rbind(xoptr,xopt)
    ret <- list(xx=xi,yy=yi,xopt=xopt,xoptr=xoptr,maxinfo=maxinfo,stopflag=stopflag)
    return(ret)
}
ojsinvseqtr <- function(xi,yi,yobs,nadd,feasible,grid,mtype=c("zmean","cmean","lmean"),
                        func,...,relthres=0,difthres=0,difstep=1,dxithres=0,dxistep=1,
                        d=NULL,g=0.001)
{
    xi <- as.matrix(xi)
    nd <- ncol(xi)
    lw <- sqrt(apply((yi-yobs)^2,2,sum))
    tlen <- length(yobs)
    wmin <- min(lw)
    nfea <- nrow(feasible)
    xoptr <- NULL
    maxinfo <- NULL
    thres <- 0
    difqueue <- initQueue(difstep)
    dxiqueue <- initQueue(dxistep)
    val <- var(yobs)*(tlen-1)
    stopflag <- 0
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
        if(i == 1) thres <- relthres*mm
        if(mm<thres)
        {
            stopflag <- 1
            break
        }
        difratio <- evalRatio(difqueue,mm)
        if(isFull(difqueue) && difratio<difthres)
        {
            stopflag <- 2
            break
        }
        difqueue <- enQueue(difqueue,mm)

        lwhat <- py$mean[-(1:nfea)]
        cxopt <- grid[which.min(lwhat),]
        xoptr <- rbind(xoptr,cxopt)
        newidx <- which.max(info)
        newx <- feasible[newidx,]
        newy <- func(newx,...)
        xi <- rbind(xi,newx)
        yi <- cbind(yi,newy)
        cdev <- sum((yobs-func(cxopt,...))^2)/val
        dxiratio <- evalRatio(dxiqueue,cdev)
        if(isFull(dxiqueue) && dxiratio < dxithres)
        {
            stopflag <- 3
            break
        }
        dxiqueue <- enQueue(dxiqueue,cdev)

        feasible <- feasible[-newidx,,drop=FALSE]
        newlw <- sqrt(sum((newy-yobs)^2))
        wmin <- min(wmin,newlw)
        lw <- c(lw,newlw)
        nfea <- nfea-1
    }
    xopt <- ojsinv(xi,yi,yobs,grid,mtype,d=d,g=g)
    xoptr <- rbind(xoptr,xopt)
    ret <- list(xx=xi,yy=yi,xopt=xopt,xoptr=xoptr,maxinfo=maxinfo,stopflag=stopflag)
    return(ret)
}
lrinvseqtr <- function(xi,yi,yobs,timepoints,nadd,feasible,grid,
                       mtype=c("zmean","cmean","lmean"),
                       func,...,relthres=0,difthres=0, difstep=1,
                       dxithres=0, dxistep=1,d=NULL,g=0.001,gl=0.1,nthread=4)
{
    mtype <- match.arg(mtype)
    tlen <- length(yobs)
    xi <- as.matrix(xi)
    nd <- ncol(xi)
    delta <- yi-yobs
    conlik <- -0.5*tlen*log(colMeans(delta^2))
    tinput <- as.matrix(timepoints)
    fname <- paste("evallik",if(mtype=="zmean")"zm" else "lm",sep="_")
    evallik <- get(fname)
    cl <- parallel::makeCluster(nthread)
    uclik <- tryCatch(parallel::parApply(cl,delta,2,evallik,tinput,mtype,d,gl),
                      finally=parallel::stopCluster(cl))
    likratio <- -2*(conlik-uclik)
    wmin <- min(likratio)
    nfea <- nrow(feasible)
    xoptr <- NULL
    maxinfo <- NULL
    thres <- 0
    difqueue <- initQueue(difstep)
    dxiqueue <- initQueue(dxistep)
    val <- var(yobs)*(tlen-1)
    stopflag <- 0
    for(i in 1:nadd)
    {
        tfea <- rbind(feasible,grid)
        gpobj <- if(mtype=="zmean") gpsepms(likratio,xi,d,g) else gpseplmms(likratio,xi,mtype,d,g)
        py <- predict(gpobj,tfea)
        delete(gpobj)
        pyfea <- list(mean=py$mean[1:nfea],s2=py$s2[1:nfea])
        info <- mininfo(pyfea,wmin)
        mm <- max(info,na.rm=TRUE)
        maxinfo <- c(maxinfo,mm)
        if(i == 1) thres <- relthres*mm
        if(mm<thres)
        {
            stopflag <- 1
            break
        }
        difratio <- evalRatio(difqueue,mm)
        if(isFull(difqueue) && diffratio<difthres)
        {
            stopflag <- 2
            break
        }
        difqueue <- enQueue(difqueue,mm)

        lrhat <- py$mean[-(1:nfea)]
        cxopt <- grid[which.min(lrhat),]
        xoptr <- rbind(xoptr,cxopt)
        newidx <- which.max(info)
        newx <- feasible[newidx,]
        newy <- func(newx,...)
        xi <- rbind(xi,newx)
        yi <- cbind(yi,newy)
        cdev <- sum((yobs-func(cxopt,...))^2)/val
        dxiratio <- evalRatio(dxiqueue,cdev)
        if(isFull(dxiqueue) && dxiratio < dxithres)
        {
            stopflag <- 3
            break
        }
        dxiqueue <- enQueue(dxiqueue,cdev)

        feasible <- feasible[-newidx,,drop=FALSE]
        newdelta <- newy-yobs
        newuclik <- evallik(newdelta,tinput,mtype,d,gl)
        newconlik <- -0.5*tlen*log(mean(newdelta^2))
        newlikratio <- -2*(newconlik-newuclik)
        wmin <- min(wmin,newlikratio)
        likratio <- c(likratio,newlikratio)
        nfea <- nfea-1
    }
    xopt <- lrinv(xi,yi,yobs,timepoints,grid,mtype,d=d,g=g,gl=gl)
    xoptr <- rbind(xoptr,xopt)
    ret <- list(xx=xi,yy=yi,xopt=xopt,xoptr=xoptr,maxinfo=maxinfo,stopflag=stopflag)
    return(ret)
}
