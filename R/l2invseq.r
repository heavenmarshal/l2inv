naiveinfo <- function(py,alpha,cht,valist)
{
    infomat <- alpha*py$coeffs2-(py$coeff-cht)^2
    info <- infomat*py$d2
    info <- apply(info,2,sum)
}
norminfo <- function(py,alpha,cht,valist)
{
    varp <- apply(py$coeffs2*py$d2,2,sum)
    devp <- apply(py$d2*(py$coeff-cht)^2,2,sum)
    stdvarp <- (varp-mean(varp))/sd(varp)
    stddevp <- (devp-mean(devp))/sd(devp)
    info <- alpha*stdvarp - stddevp
}
mvconinfo <- function(py,alpha,cht,valist)
{
    barrier <- alpha * apply(py$coeffs2*py$d2,2,sum)
    numbas <- nrow(py$coeff)
    nfea <- ncol(py$coeff)
    out <- .C("mvcon_R", as.double(py$coeff), as.double(py$coeffs2),
              as.double(barrier), as.double(cht), as.double(py$d2),
              as.double(nfea), as.double(valist$nmc), as.integer(numbas),
              ans=double(nfea), PACKAGE="l2inv")
    info <- out$ans
}
coneinfo <- function(py,alpha,cht,valist)
{
    mse <- py$coeffs2
    se <- sqrt(mse)
    norm <- (cht-py$coeff)/se
    u2 <- norm+alpha
    u1 <- norm-alpha
    du2 <- dnorm(u2)
    du1 <- dnorm(u1)
    part1 <- mse*(alpha*alpha-norm*norm-1)*(pnorm(u2)-pnorm(u1))
    part2 <- mse*(u2*du2-u1*du1)
    part3 <- -2*norm*mse*(du2-du1)
    info <- part1+part2+part3
    info <- ifelse(is.na(info) | is.infinite(info),0,info)
    info <- apply(info*py$d2,2,sum)
}
l2invseq <- function(xi,yi,yobs,nadd,feasible,grid,alpha,
                     func,...,type=c("naive","norm","mvcon","cone"),
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
        newy <- func(newx,...)
        xi <- rbind(xi,newx)
        yi <- cbind(yi,newy)
        feasible <- feasible[-newidx,]
    }
    xopt <- l2inv(xi,yi,yobs,grid,frac,d=d,g=g)
    ret <- list(xx=xi,yy=yi,xopt=xopt)
    return(ret)
}
mininfo <- function(py,cmin)
{
    ave <- py$mean
    se <- sqrt(py$s2)
    norm <- (cmin-ave)/se
    info <- ifelse(is.infinite(norm),0,
                   se*(norm*pnorm(norm)+dnorm(norm)))
}
ojsinvseq <- function(xi,yi,yobs,nadd,feasible,grid,mtype=c("zmean","cmean","lmean"),
                      func,...,d=NULL,g=0.001)
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
        newy <- func(newx,...)
        xi <- rbind(xi,newx)
        yi <- cbind(yi,newy)
        feasible <- feasible[-newidx,]
        newlw <- sqrt(sum((newy-yobs)^2))
        wmin <- min(wmin,newlw)
        lw <- c(lw,newlw)
    }
    xopt <- ojsinv(xi,yi,yobs,grid,mtype,d=d,g=g)
    ret <- list(xx=xi,yy=yi,xopt=xopt)
    return(ret)
}
lrinvseq <- function(xi,yi,yobs,timepoints,nadd,feasible,grid,
                     mtype=c("zmean","cmean","lmean"),
                     func,...,d=NULL,g=0.001,gl=0.1,nthread=4)
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
        newy <- func(newx,...)
        xi <- rbind(xi,newx)
        yi <- cbind(yi,newy)
        feasible <- feasible[-newidx,]
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
