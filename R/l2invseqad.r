l2invseqad <- function(xi,yi,yobs,nadd,feasible,grid,alphafrom,alphato,
                       func,...,type=c("naive","cone"),
                       frac=.95,d=NULL,g=0.001,nthread=4)
{
    xi <- as.matrix(xi)
    type <- match.arg(type)
    infoname <- paste(type,"info",sep="")
    infofun <- get(infoname)
    alpha <- seq(alphafrom,alphato,len=nadd)
    for(i in 1:nadd)
    {
        lbasis <- buildBasis(yi,frac)
        cht <- drop(t(lbasis$basis)%*%yobs/lbasis$redd^2)
        py <- svdgpsepms(feasible,xi,yi,frac,nthread=nthread)
        info <- infofun(py,alpha[i],cht)
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
