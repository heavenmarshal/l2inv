esl2d <- function(x,gpobjs,nbas,reddsq,cht)
{
    x <- matrix(x,nrow=1)
    pmean <- ps2 <- vector("double",nbas)
    for(i in 1:nbas)
    {
        pred <- predict(gpobjs[[i]],x)
        pmean[i] <- pred$mean
        ps2[i] <- pred$s2
    }
    crit <- sum(reddsq*(ps2+(pmean-cht)^2))
}
naive <- function(x,gpobjs,nbas,reddsq,cht)
{
    x <- matrix(x,nrow=1)
    pmean <- ps2 <- vector("double",nbas)
    for(i in 1:nbas)
    {
        pred <- predict(gpobjs[[i]],x)
        pmean[i] <- pred$mean
        ps2[i] <- pred$s2
    }
    crit <- sum(reddsq*(pmean-cht)^2)
}
ojsdev <- function(x,gpobj,q)
{
    x <- matrix(x,nrow=1)
    py <- predict(gpobj,x)
    crit <- py$mean + qnorm(q)*sqrt(py$s2)
}

esl2dzmean <- function(x,gpobjs,nbas,reddsq,cht)
{
    gpidces <- sapply(gpobjs,`[`,"optgpi")
    out <- .C("esl2dZmean", as.integer(nbas), as.integer(gpidces),
              as.double(x), as.double(reddsq), as.double(cht),
              esl2d=double(1))
    return(out$esl2d)
}

esl2dcmean <- function(x,gpobjs,nbas,reddsq,cht)
{
    gpidces <- sapply(gpobjs,`[`,"optgpi")
    out <- .C("esl2dCmean", as.integer(nbas), as.integer(gpidces),
              as.double(x), as.double(reddsq), as.double(cht),
              esl2d=double(1))
    return(out$esl2d)
}
esl2dlmean <- function(x,gpobjs,nbas,reddsq,cht)
{
    ndim <- length(x)
    gpidces <- sapply(gpobjs,`[`,"optgpi")
    out <- .C("esl2dLmean", as.integer(nbas), as.integer(ndim),
              as.integer(gpidces), as.double(x), as.double(reddsq),
              as.double(cht), esl2d=double(1))
    return(out$esl2d)
}
