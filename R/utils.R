invlogit <- function(x) 1/(1+exp(-x))
rbern <- function(p) drop0((p > Matrix(runif(length(p)),nrow(p),ncol(p)))*1 )

testMatrixRowwiseConversion <- function() {
  T <- 1000
  m <- Matrix(0,T,T)
  m[sample(1:(T^2),T)] <- 1
  b <- Matrix2rowwiselist(m)
  a <- rowwiselist2Matrix(b$rows,b$cols,dim(m))
  stopifnot(all(m==a))
}

rowwiselist2Matrix <- function(rows,cols,dims) {
  # Use dims=c(1,1) if want matrix to fit exactly to max of rows and columns
  tmp <- list()
  for (i in 1:length(rows)) {
    tmp[[i]] <- cbind(rows[i],cols[[i]])
  }
  subs <- do.call(rbind,tmp)
  T <- max(max(subs[,1]),dims[1])
  N <- max(max(subs[,2]),dims[2])
  x <- Matrix(0,T,N)
  x[sub2ind(subs,c(T,N))] <- 1
  return(x)
}

Matrix2rowwiselist <- function(m) {
  tmp <- as(m,"dgTMatrix")
  urows <- unique(tmp@i)
  cols <- list()
  for (i in 1:length(urows)) {
    ix <- which(tmp@i == urows[i])
    cols[[i]] <- tmp@j[ix] + 1 # tmp$j is 0 based
  }
  urows <- urows + 1
  o <- order(urows)
  return(list(rows=urows[o],cols=cols[o]))
}
empirical.dist.of.sets <- function(y) {
  yr <- Matrix2rowwiselist(y)$cols  
  ys <- unique(yr)
  id <- match(yr, ys)
  tb <- table(id)
  N <- max(sapply(yr,max))
  events <- rowwiselist2Matrix(1:length(ys),ys,c(length(ys),N))
  o <- order(as.vector(tb),decreasing=TRUE)
  events <- events[o,]
  tb <- tb[o]/sum(tb)
  return(list(x=events,prob=tb))
}

rownormalize <- function(x) {
  if (is.null(dim(x)))
    return(x/sum(x))
  rs <- rowSums(x)
  tmp <- matrix(rs,nrow(x),ncol(x))
  x / tmp
}
rownormalizelogs <- function(x) {
  I <- nrow(x)
  J <- ncol(x)
  ix <- sub2ind(cbind(1:I,max.col(x)),c(I,J))
  rowmaxs <- matrix(x[ix],I,J)
  x <- rownormalize(exp(x - rowmaxs))  # for numerical stability
  return(x)
}
sub2ind <- function(subs, dims) {
  # matrix where each row for a different element,
  # each column is for a particular dimension
  D <- length(dims)
  inds <- subs[,1]
  for (d in (1:D)[-1]) {
    inds <- inds + (subs[,d]-1) * prod(dims[-(d:D)])
  }
  inds
}
ind2sub <- function(ind,dims) {
  # Only supports matrices right now
#dims <- c(6,4)
#all(sub2ind(ind2sub(1:30,dims),dims)==1:30)
  D <- length(dims)
  subs <- matrix(0,length(ind),D)
  subs[,1] <- ((ind-1) %% dims[1])+1
  subs[,2] <- ((ind-1) %/% dims[1])+1
  return(subs)
}
rbern <- function(n,p)  (runif(n) < p) * 1
rdiscrete <- function(n,p) {
  sample(1:length(p),size=n,prob=p)
}

logOfSum <- function(logValues) {
  a0 <- max(logValues)
  index <- which(logValues == a0)
  logValues <- logValues[-index]
  accum <- 1
  for (i in 1:length(logValues)) {
     accum <- accum + exp(logValues[i] - a0)
  }
  a0 + log(accum)
}
samplelp <- function(pz) {
  # Given a vector of log values, sample an index using them as probability weights.
  pz <- exp(pz - max(pz))  # for numerical stability
  pz <- pz/sum(pz)
  rdiscrete(1,pz)
}
samplelpmat <- function(pz) {
  pz <- rownormalizelogs(pz)
  x <- (runif(nrow(pz)) < pz[,2])*1
  return(x)
}


llk <- function(y,dims,params) {
  lk <- lpy(y,dims,params,wsum="all")
  z.i <- rowSums(params$z)
  w.i <- rowSums(params$w)
  lk <- lk +
    sum(z.i * log(params$phi)) +
    sum((dims$K-z.i) * log(1-params$phi)) + 
    sum(w.i * log(params$omega)) +
    sum((dims$K-w.i) * log(1-params$omega))
  return(lk)
}

plotProgress <- function(fit) {
# TODO: Need to fix scaling on jjplot tile().  Right now just sets one of the elements to 1.

  y <- fit$y
  if (any(dim(y)>100)) error("Progress plots on a problem of this size not recommended")
  params <- fit$params
  T <- nrow(y)
  N <- ncol(y)
  theta <- params$theta
  colnames(theta) <- colnames(z)
  rownames(theta) <- rownames(z)
  z <- params$z
  rownames(z) <- colnames(y)
  colnames(z) <- paste("set",0:(K-1))
  w <- params$w
  rownames(w) <- rownames(y)
  colnames(w) <- paste("set",0:(K-1))
  subs <- ind2sub(1:(T*N),c(T,N))
  yhat <- matrix(predict.lsm(fit=fit,subs=subs),T,N)
  colnames(yhat) <- colnames(y)
  rownames(yhat) <- rownames(y)
  yhat <- Matrix(yhat)
  
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(2, 3)))
  pushViewport(viewport(layout.pos.col=1,  layout.pos.row=1))
  disp(z)
  popViewport()
  pushViewport(viewport(layout.pos.col=2,  layout.pos.row=1))
  disp(theta)
  popViewport()
  pushViewport(viewport(layout.pos.col=3,  layout.pos.row=1))
  disp(z*theta)
  popViewport()
  pushViewport(viewport(layout.pos.col=1,  layout.pos.row=2))
  disp(w)
  popViewport()
  pushViewport(viewport(layout.pos.col=2,  layout.pos.row=2))
  disp(t(y))
  popViewport()
  pushViewport(viewport(layout.pos.col=3,  layout.pos.row=2))
#  yhat[1,1] <- 0  # Temporary hack
#  yhat[1,2] <- 1
  disp(t(yhat))
  popViewport()

}

disp <- function (z, rowOrderNames = NULL, colOrderNames = NULL, limits = c(0,1)) 
{

    require(jjplot)
    require(reshape)
    
    if (class(z) != "matrix") {
        cat("converting to matrix...\n")
        z <- as.matrix(t(z))
    }
    df <- melt(z)
    
    if (is.null(rowOrderNames) & !is.null(rownames(z))) {
        rowOrderNames <- rownames(z)
    }
    if (is.null(colOrderNames) & !is.null(colnames(z))) {
        colOrderNames <- colnames(z)
    }

    if (!is.null(rowOrderNames)) {
        df$X1 <- reorder(df$X1, match(as.character(df$X1), rowOrderNames))
    }
    if (!is.null(colOrderNames)) {
        df$X2 <- reorder(df$X2, match(as.character(df$X2), colOrderNames))
    }
    
    jjplot(X2 ~ tile(border=NA) : color(value) + X1,  data = df, newpage=FALSE,theme=jjplot.theme("bw"),xlab="",ylab="")
}
