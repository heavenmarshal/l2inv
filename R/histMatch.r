## scalar valued implausible note: observation
## value is writen in the head of vector resp
implausible <- function(resp,design,feasible,mtype,sig2,d,g)
{
    obsval <- resp[1]
    resp <- resp[-1]
    gpobj <- if(mtype=="zmean") gpsepms(resp,design,d,g) else gpseplmms(resp,design,mtype,d,g)
    py <- predict(gpobj,feasible)
    delete(gpobj)
    impl <- abs(py$mean-obsval)/sqrt(py$s2+sig2)
    return(impl)
}
## maximum implausible criterion for the time series
maximp <- function(resp,design,yobs,feasible,mtype,sig2,
                   d,g,nthread)
{
    argresp <- cbind(yobs,resp)
    cl <- parallel::makeCluster(nthread)
    impl <- tryCatch(parallel::parApply(cl,argresp,1,implausible,
                                        design,feasible,mtype,sig2,d,g),
                     finally=parallel::stopCluster(cl))
    maximpl <- apply(impl,1,max)
    return(maximpl)
}
histMatch <- function(xi,yi,yobs,timepoints,cutoff,budget,nfea,
                      mtype=c("zmean","cmean","lmean"),
                      featype=c("maximin","random"),
                      func,...,d=NULL,g=0.001,nthread=4)
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
        feasible <- lhsgen(nfea,din)
        impval <- maximp(yi,xi,yobs,feasible,mtype,sig2hat,d,g,nthread)
        nnonimp <- sum(impval<cutoff)
        if(nnonimp <= 0) break
        nadd <- min(nnonimp,nrem)
        sidx <- order(impval)[1:nadd]
        newx <- feasible[sidx,,drop=FALSE]
        newy <- apply(newx,1,func,...)
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
