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
getcols <- function(idx,mat) mat[,idx]
rename <- function(lst,name) {names(lst) <- name; return(lst)}

## a queue for diff stopping criterion

initQueue <- function(capacity)
{
    if(capacity<=0) error("capacity of a queue must be positive!")
    qv <- rep(NA,capacity)
    queue <- list(qv=qv,capacity=capacity,length=0)
    return(queue)
}
isFull <- function(queue)
{
    return(queue$capacity==queue$length)
}
## if is full automatically discard the end item
enQueue <- function(queue,x)
{
    if(queue$length==0 || queue$capacity==1)
    {
        queue$qv[1]=x
        queue$length=1
        return(queue)
    }
    qv <- queue$qv
    if(isFull(queue))
    {
        shiftidx <- 1:(queue$capacity-1)
        qv[shiftidx+1] <- qv[shiftidx]
        qv[1] <- x
        queue$qv <- qv
        return(queue)
    }
    len <- queue$length
    qv[2:(len+1)] <- qv[1:len]
    qv[1] <- x
    queue$qv <- qv
    queue$length <- len+1
    return(queue)
}
evalRatio <- function(queue,newx)
{
    if(queue$length==0) return(Inf)
    qv <- queue$qv
    denominator <- sum(abs(qv),na.rm=TRUE)
    extqv <- c(newx,qv)
    difeqv <- abs(diff(extqv))
    numerator <- sum(difeqv,na.rm=TRUE)
    return(numerator/denominator)
}
