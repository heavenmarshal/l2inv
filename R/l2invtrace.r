## trace the update for the sequential design, for experiment not application.
l2invseqtr <- function(xi,yi,yobs,nadd,feasible,grid,alpha,
                       func,...,type=c("mvapp","mvei","projei","oei"),
                       mtype=c("zmean","cmean","lmean"),relthres=0,
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
    xoptr <- matrix(nrow=nadd+1,ncol=nd)
    valist.default <- list(nmc=500)
    remnames <- setdiff(names(valist.default),names(valist))
    valist <- c(valist,valist.default[remnames])
    valist$yobs <- yobs
    maxinfo <- NULL
    thres <- 0
    for(i in 1:nadd)
    {
        tfea <- rbind(feasible,grid)
        lbasis <- buildBasis(yi,frac)
        cht <- drop(t(lbasis$basis)%*%yobs/lbasis$redd^2)
        valist$chts2 <- drop(crossprod(yobs-lbasis$basis%*%cht))/tlen
        valist$mindist <- min(apply(lbasis$redd^2*(t(lbasis$coeff)-cht)^2,2,sum))
        valist$barval <-  min(apply((yobs-yi)^2,2,sum))
        py <- svdgpsepms(tfea,xi,yi,frac,mtype=mtype,nthread=nthread)
        ## extract the part for feasible
        pyfea <- list(d2=py$d2, coeffs2=py$coeffs2[,1:nfea],coeff=py$coeff[,1:nfea],basis=py$basis,varres=py$varres)
        pygrid <- list(d2=py$d2, coeffs2=py$coeffs2[,-(1:nfea)],coeff=py$coeff[,-(1:nfea)])
        numbas <- lbasis$numbas
        info <- infofun(pyfea,alpha,cht,valist)
        mm <- max(info,na.rm=TRUE)
        maxinfo <- c(maxinfo,mm)
        if(i == 1) thres <- relthres*mm
        if(mm<thres) break
        if(numbas==1)
            criter <- pygrid$d2*(pygrid$coeffs2+(pygrid$coeff-cht)^2)
        else
            criter <- apply(pygrid$d2*(pygrid$coeffs2+(pygrid$coeff-cht)^2),2,sum)
        xoptr[i,] <- grid[which.min(criter),]
        newidx <- which.max(info)
        newx <- feasible[newidx,]
        newy <- func(newx,...)
        xi <- rbind(xi,newx)
        yi <- cbind(yi,newy)
        feasible <- feasible[-newidx,,drop=FALSE]
        nfea <- nfea-1
    }
    xopt <- l2inv(xi,yi,yobs,grid,frac,d=d,g=g)
    xoptr[nadd+1,] <- xopt
    ret <- list(xx=xi,yy=yi,xopt=xopt,xoptr=xoptr,maxinfo=maxinfo)
    return(ret)
}
ojsinvseqtr <- function(xi,yi,yobs,nadd,feasible,grid,mtype=c("zmean","cmean","lmean"),
                        func,...,relthres=0,d=NULL,g=0.001)
{
    xi <- as.matrix(xi)
    nd <- ncol(xi)
    lw <- sqrt(apply((yi-yobs)^2,2,sum))
    wmin <- min(lw)
    nfea <- nrow(feasible)
    xoptr <- matrix(nrow=nadd+1,ncol=nd)
    maxinfo <- NULL
    thres <- 0
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
        if(info<thres) break
        lwhat <- py$mean[-(1:nfea)]
        xoptr[i,] <- grid[which.min(lwhat),]
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
    xopt <- ojsinv(xi,yi,yobs,grid,mtype,d=d,g=g)
    xoptr[nadd+1,] <- xopt
    ret <- list(xx=xi,yy=yi,xopt=xopt,xoptr=xoptr,maxinfo=maxinfo)
    return(ret)
}
lrinvseqtr <- function(xi,yi,yobs,timepoints,nadd,feasible,grid,
                     mtype=c("zmean","cmean","lmean"),
                     func,...,relthres=0,d=NULL,g=0.001,gl=0.1,nthread=4)
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
    xoptr <- matrix(nrow=nadd+1,ncol=nd)
    maxinfo <- NULL
    thres <- 0
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
        if(mm<thres) break
        lrhat <- py$mean[-(1:nfea)]
        xoptr[i,] <- grid[which.min(lrhat),]
        newidx <- which.max(info)
        newx <- feasible[newidx,]
        newy <- func(newx,...)
        xi <- rbind(xi,newx)
        yi <- cbind(yi,newy)
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
    xoptr[nadd+1,] <- xopt
    ret <- list(xx=xi,yy=yi,xopt=xopt,xoptr=xoptr,maxinfo=maxinfo)
    return(ret)
}
