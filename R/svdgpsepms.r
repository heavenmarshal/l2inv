fitgps_zmean <- function(resp,design,test,nstarts=5,d=NULL,g=0.001,type="zmean")
{
    din <- ncol(design)
    d <- darg(d,design)
    g <- garg(g,resp)
    parb <- lhs::maximinLHS(nstarts,din+1)
    ranb <- parb[,1:din,drop=FALSE]
    nugb <- parb[,din+1]
    ldmin <- log(d$min)
    ldmax <- log(d$max)
    ranb <- ldmin + ranb*(ldmax-ldmin)
    ranb <- exp(ranb)
    lgmin <- log(g$min)
    lgmax <- log(g$max)
    nugb <- lgmin + nugb*(lgmax-lgmin)
    nugb <- exp(nugb)
    gpis <- llik <- rep(NA,nstarts)
    for(i in 1:nstarts)
    {
        gpis[i] <- newGPsep(design,resp,ranb[i,],nugb[i],TRUE)
        mle <- jmleGPsep(gpis[i],drange=c(d$min,d$max),
                               grange=c(g$min,g$max),
                               dab=d$ab,gab=g$ab)
        llik[i] <- llikGPsep(gpis[i],dab=d$ab,gab=g$ab)
    }
    optidx <- which.max(llik)
    pred <- predGPsep(gpis[optidx],test,TRUE)
    for(i in 1:nstarts) deleteGPsep(gpis[i])
    return(pred)
}
fitgps_lmean <- function(resp,design,test,nstarts=5,d=NULL,g=0.001,type=c("cmean","lmean"))
{
    din <- ncol(design)
    d <- darg(d,design)
    g <- garg(g,resp)
    parb <- lhs::maximinLHS(nstarts,din+1)
    ranb <- parb[,1:din,drop=FALSE]
    nugb <- parb[,din+1]
    ldmin <- log(d$min)
    ldmax <- log(d$max)
    ranb <- ldmin + ranb*(ldmax-ldmin)
    ranb <- exp(ranb)
    lgmin <- log(g$min)
    lgmax <- log(g$max)
    nugb <- lgmin + nugb*(lgmax-lgmin)
    nugb <- exp(nugb)
    gpis <- llik <- rep(NA,nstarts)
    for(i in 1:nstarts)
    {
        gpis[i] <- newGPsepLm(design,resp,ranb[i,],nugb[i],TRUE,type)
        mle <- jmleGPsepLm(gpis[i],drange=c(d$min,d$max),
                           grange=c(g$min,g$max),
                           dab=d$ab,gab=g$ab)
        llik[i] <- llikGPsepLm(gpis[i],dab=d$ab,gab=g$ab)
    }
    optidx <- which.max(llik)
    pred <- predGPsepLm(gpis[optidx],test,type)
    for(i in 1:nstarts) deleteGPsepLm(gpis[i])
    return(pred)

}
svdgpsepms <- function(X0,design,resp,frac=.95,nstart=5,mtype=c("zmean","cmean","lmean"),
                       d=NULL,g=0.001,nthread=4)
{
    mtype <- match.arg(mtype)
    fitname <- if(mtype=="cmean") "lmean" else mtype
    fitgps <- get(paste("fitgps",fitname,sep="_"))
    lenresp <- length(resp)
    lbasis <- buildBasis(resp,frac)
    basis <- lbasis$basis
    numbas <- lbasis$numbas
    coeff <- lbasis$coeff
    resid <- resp-basis%*%t(coeff)
    varres <- drop(crossprod(as.vector(resid)))
    varres <- varres/(lenresp+2)
    cl <- parallel::makeCluster(nthread)
    ret <- tryCatch(parallel::parApply(cl,coeff,2,fitgps,
                                       design,X0,nstart,d,g,mtype),
                    finally=parallel::stopCluster(cl))
    vmean <- matrix(unlist(sapply(ret,`[`,"mean")),nrow=numbas,byrow=TRUE)
    vsigma2 <- matrix(unlist(sapply(ret,`[`,"s2")),nrow=numbas,byrow=TRUE)
    pmean <- basis%*%vmean
    psd <- sqrt(basis^2%*%vsigma2+varres)
    ret <- list(mean=pmean,sd=psd,coeff=vmean,coeffs2=vsigma2,d2=lbasis$redd^2)
    return(ret)
}
