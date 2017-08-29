l2invseqset <- function(xi,yi,yobs,nadd,feasible,grid,fearesp,
                        alpha,type=c("naive","norm","mvcon","cone"),
                        frac=.95,d=NULL,g=0.001,nmc=500,nthread=4)
{
    xi <- as.matrix(xi)
    type <- match.arg(type)
    infoname <- paste(type,"info",sep="")
    infofun <- get(infoname)
    valist <- list(nmc=nmc)
    for(i in 1:nadd)
    {
        lbasis <- buildBasis(yi,frac)
        cht <- drop(t(lbasis$basis)%*%yobs/lbasis$redd^2)
        py <- svdgpsepms(feasible,xi,yi,frac,nthread=nthread)
        info <- infofun(py,alpha,cht,valist)
        newidx <- which.max(info)
        newx <- feasible[newidx,]
        newy <- fearesp[,newidx]
        xi <- rbind(xi,newx)
        yi <- cbind(yi,newy)
        feasible <- feasible[-newidx,]
        fearesp <- fearesp[,-newidx]
    }
    xopt <- l2inv(xi,yi,yobs,grid,frac,d=d,g=g)
    ret <- list(xx=xi,yy=yi,xopt=xopt)
    return(ret)
}
ojsinvseqset <- function(xi,yi,yobs,nadd,feasible,grid,fearesp,
                         mtype=c("zmean","cmean","lmean"),
                         d=NULL,g=0.001)
{
    xi <- as.matrix(xi)
    lw <- sqrt(apply((yi-yobs)^2,2,sum))
    wmin <- min(lw)
    for(i in 1:nadd)
    {
        gpobj <- if(mtype=="zmean") gpsepms(lw,xi,d,g) else gpseplmms(lw,xi,mtype,d,g)
        py <- predict(gpobj,feasible)
        delete(gpobj)
        info <- mininfo(py,wmin)
        newidx <- which.max(info)
        newx <- feasible[newidx,]
        newy <- fearesp[,newidx]
        xi <- rbind(xi,newx)
        yi <- cbind(yi,newy)
        feasible <- feasible[-newidx,]
        fearesp <- fearesp[,-newidx]
        newlw <- sqrt(sum((newy-yobs)^2))
        wmin <- min(wmin,newlw)
        lw <- c(lw,newlw)
    }
    xopt <- ojsinv(xi,yi,yobs,grid,mtype,d=d,g=g)
    ret <- list(xx=xi,yy=yi,xopt=xopt)
    return(ret)
}
lrinvseqset <- function(xi,yi,yobs,timepoints,nadd,feasible,grid,
                        fearesp,mtype=c("zmean","cmean","lmean"),
                        d=NULL,g=0.001,gl=0.1,nthread=4)
{
    mtype <- match.arg(mtype)
    tlen <- length(yobs)
    xi <- as.matrix(xi)
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
    for(i in 1:nadd)
    {
        gpobj <- if(mtype=="zmean") gpsepms(likratio,xi,d,g) else gpseplmms(likratio,xi,mtype,d,g)
        py <- predict(gpobj,feasible)
        delete(gpobj)
        info <- mininfo(py,wmin)
        newidx <- which.max(info)
        newx <- feasible[newidx,]
        newy <- fearesp[,newidx]
        xi <- rbind(xi,newx)
        yi <- cbind(yi,newy)
        feasible <- feasible[-newidx,]
        fearesp <- fearesp[,-newidx]
        newdelta <- newy-yobs
        newuclik <- evallik(newdelta,tinput,mtype,d,gl)
        newconlik <- -0.5*tlen*log(mean(newdelta^2))
        newlikratio <- -2*(newconlik-newuclik)
        wmin <- min(wmin,newlikratio)
        likratio <- c(likratio,newlikratio)
    }
    xopt <- lrinv(xi,yi,yobs,timepoints,grid,mtype,d=d,g=g,gl=gl)
    ret <- list(xx=xi,yy=yi,xopt=xopt)
    return(ret)
}
histMatchset <- function(xi,yi,yobs,timepoints,cutoff,budget,nfea,
                         feasible, fearesp, mtype=c("zmean","cmean","lmean"),
                         d=NULL,g=0.001,nthread=4)
{
    mtype <- match.arg(mtype)
    featype <- match.arg(featype)
    lhsname <- paste(featype,"LHS",sep="")
    lhsgen <- get(lhsname)
    xi <- as.matrix(xi)
    din <- ncol(xi)
    ndes <- nrow(xi)
    sig2hat <- ssanova(yobs~timepoints,type="cubic",method="v")$varht
    nrem <- budget-ndes
    iteradd <- NULL
    while(nrem > 0)
    {
        impval <- maximp(yi,xi,yobs,feasible,mtype,sig2hat,d,g,nthread)
        nnonimp <- sum(impval<cutoff)
        if(nnonimp <= 0) break
        nadd <- min(nnonimp,nrem)
        sidx <- order(impval)[1:nadd]
        newx <- feasible[sidx,,drop=FALSE]
        newy <- fearesp[,sidx,drop=FALSE]
        feasible <- feasible[-sidx,,drop=FALSE]
        fearesp <- fearesp[,-sidx,drop=FALSE]
        iteradd <- c(iteradd,nadd)
        xi <- rbind(xi,newx)
        yi <- cbind(yi,newy)
        nrem <- nrem-nadd
    }
    devs <- sqrt(apply((yobs-yi)^2,2,sum))
    optidx <- which.min(devs)
    xopt <- xi[optidx,]
    yopt <- yi[,optidx]
    ret <- list(xx=xi,yy=yi,iteradd=iteradd,xopt=xopt,yopt=yopt)
    return(ret)
}
