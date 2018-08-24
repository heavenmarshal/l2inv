esl2dinvoptR <- function(design,resp,yobs,frac=.95,
                        mtype=c("zmean","cmean","lmean"),
                        critype=c("esl2d","naive"),
                        d=NULL,g=0.001)
{
    ndim <- ncol(design)
    mtype <- match.arg(mtype)
    critype <- match.arg(critype)
    design <- as.matrix(design)
    lbasis <- buildBasis(resp,frac)
    cht <- drop(t(lbasis$basis)%*%yobs/lbasis$redd^2)
    numbas <- lbasis$numbas
    critfun <- get(critype)
    coeff <- lbasis$coeff
    gpobjs <- vector("list",numbas)
    reddsq <- lbasis$redd^2
    for(i in 1:numbas)
    {
        ccoeff <- coeff[,i]
        gpobjs[[i]] <- if(mtype=="zmean") gpsepms(ccoeff,design,d,g)
                       else gpseplmms(ccoeff,design,mtype,d,g)
    }
    opt <- genoud(critfun,ndim,Domains=cbind(rep(0,ndim),rep(1,ndim)),
                  boundary.enforcement=2,print.level=0,
                  gpobjs=gpobjs, nbas=numbas, reddsq=reddsq,cht=cht)
    xopt <- opt$par
    for(i in 1:numbas) delete(gpobjs[[i]])
    return(xopt)
}

esl2dinvopt <- function(design,resp,yobs,frac=.95,
                        mtype=c("zmean","cmean","lmean"),
                        d=NULL,g=0.001,gactl=list())
{
    ndim <- ncol(design)
    mtype <- match.arg(mtype)
    design <- as.matrix(design)
    lbasis <- buildBasis(resp,frac)
    cht <- drop(t(lbasis$basis)%*%yobs/lbasis$redd^2)
    numbas <- lbasis$numbas
    critfun <- get(paste("esl2d",mtype,sep=""))
    coeff <- lbasis$coeff
    gpobjs <- vector("list",numbas)
    reddsq <- lbasis$redd^2
    gactl.default <- list(pop.size=1000, max.generations = 100,
                          wait.generations = 10)
    remnames <- setdiff(names(gactl.default),names(gactl))
    gactl <- c(gactl,gactl.default[remnames])
    for(i in 1:numbas)
    {
        ccoeff <- coeff[,i]
        gpobjs[[i]] <- if(mtype=="zmean") gpsepms(ccoeff,design,d,g)
                       else gpseplmms(ccoeff,design,mtype,d,g)
    }
    opt <- genoud(critfun,ndim,pop.size=gactl$pop.size,
                  max.generations=gactl$max.generations,
                  wait.generations=gactl$wait.generations,
                  Domains=cbind(rep(0,ndim),rep(1,ndim)),
                  boundary.enforcement=2,print.level=0,
                  gpobjs=gpobjs, nbas=numbas, reddsq=reddsq,cht=cht)
    xopt <- opt$par
    for(i in 1:numbas) delete(gpobjs[[i]])
    return(xopt)
}

sl2invopt <- function(design,resp,yobs,mtype=c("zmean","cmean","lmean"),
                      q=.5,d=NULL,g=0.001,gactl=list())
{
    ndim <- ncol(design)
    mtype <- match.arg(mtype)
    design <- as.matrix(design)
    isnmis <- !is.na(yobs)
    yobsred <- yobs[isnmis]
    respred <- resp[isnmis,]
    dev <- sqrt(apply((respred-yobsred)^2,2,sum))
    gactl.default <- list(pop.size=1000, max.generations = 100,
                          wait.generations = 10)
    remnames <- setdiff(names(gactl.default),names(gactl))
    gactl <- c(gactl,gactl.default[remnames])
    gpobj <- if(mtype=="zmean") gpsepms(dev,design,d,g) else gpseplmms(dev,design,mtype,d,g)
    opt <- genoud(sl2dev,ndim,pop.size=gactl$pop.size,
                  max.generations=gactl$max.generations,
                  wait.generations=gactl$wait.generations,
                  Domains=cbind(rep(0,ndim),rep(1,ndim)),
                  boundary.enforcement=2,print.level=0,
                  gpobj=gpobj, q=q)
    xopt <- opt$par
    delete(gpobj)
    return(xopt)
}
