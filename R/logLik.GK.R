logLik.GK <- function(z) {
    w <- z$w/sum(z$w)
    uv <- expand.grid(z$u, z$v)
    B <- as.matrix(cbind(uv, 1))
    Xp <- z$X[z$y==1,]
    Xn <- z$X[z$y==0,]
    pos <- apply(Xp %*% t(B), 1, function(a) sum(w[a > 0]))
    neg <- 1 - apply(Xn %*% t(B), 1, function(a) sum(w[a > 0]))
    sum(log(c(pos, neg)))
}

