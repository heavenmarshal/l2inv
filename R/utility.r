buildBasis <- function(response,percent=.95,numbas=NULL)
{
    svdm <- svd(response)
    if(is.null(numbas))
    {
        cumd <- cumsum(svdm$d)/sum(svdm$d)
        numbas <- min(which(cumd>=percent))
    }
    basis <- svdm$u[,1:numbas]          #basis matrix TT by numbas
    redd <- svdm$d[1:numbas]            #reduced d vector length numbas
    redv <- svdm$v[,1:numbas]           #reduced v matrix nn by mumbas
    basis <- t(t(basis)*redd)            #coefficient matrix nn by numbas
    coeff <- redv
    if(!is.matrix(coeff))
        coeff <- matrix(coeff,ncol=numbas)
    ret <- list(basis=basis,redd=redd,coeff=coeff,numbas=numbas)
}

imputels <- function(yobs,lbasis,maxiter,tol)
{
    yimp <- yobs
    msidx <- which(is.na(yobs))
    yimp[msidx] <- mean(yobs,na.rm=TRUE)
    basis <- lbasis$basis
    reddsq <- lbasis$redd^2
    pcht <- rep(0,lbasis$numbas)
    for(i in 1:maxiter)
    {
        cht <- drop(t(basis)%*%yimp/reddsq)
        err <- sqrt(sum((cht-pcht)^2))
        if(err<tol)
            break
        pcht <- cht
        yimp[msidx] <- basis[msidx,]%*%cht
    }
    converge <- (i<maxiter)
    ysmooth <- drop(basis%*%cht)
    ret <- list(converge=converge,ysmooth=ysmooth,cht=cht)
    return(ret)
}
