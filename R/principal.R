### CLOSED PRINCIPAL CURVE

## Instead of the angle in PC1/PC2, fit a closed curve through the
## cells in the first k PCs (Hastie & Stuetzle 1989 iterations, with a
## Fourier series in each coordinate as the periodic smoother) and use
## the arc length along the curve as the phase.

## Fourier design matrix of order K
fourier <- function(lambda, K) {
    m <- outer(lambda, seq_len(K))
    cbind(1, cos(m), sin(m))
}

#' Phase along a closed principal curve.
#'
#' Fits a closed curve through the cells in the first \code{k} PCs,
#' each PC weighted by its standard deviation, and projects each cell
#' to the nearest point of the curve. Starting from the angle in
#' PC1/PC2, it alternates fitting a Fourier series of order \code{K} in
#' each coordinate against the cell phase, and re-projecting the cells
#' onto the curve, with the phase as arc length. The phase
#' \code{lambda} is then aligned with the current \code{theta} (after
#' \link{revert} and \link{shift}): same direction, mean offset 0.
#' Experimental: on yeast single-cell data the curve ordered cells less
#' reproducibly than the angle (ROADMAP.md); if the cells collapse onto
#' a few points of the curve, the iteration stops with a warning and
#' keeps the previous phase.
#' @inheritParams revert
#' @param k number of PCs.
#' @param K order of the Fourier series.
#' @param iter maximal number of iterations.
#' @param tol stop when the mean absolute phase change is smaller.
#' @param ngrid number of points on the curve for the projection.
#' @return the phases object with the cell phase \code{lambda} (in
#'     -pi:pi) and the distance to the curve \code{curve.dist} in
#'     \code{rotation.phase}; the curve in \code{phases$curve}
#'     (\code{grid}: phase before alignment, \code{lambda} is
#'     \code{sign*grid$lambda + offset}, and coordinates; \code{coef};
#'     \code{iter}).
#' @export
principal_phase <- function(phases, k=3, K=4, iter=20, tol=1e-3,
                            ngrid=500) {

    k <- min(k, ncol(phases$rotation))
    Y <- t(t(phases$rotation[, 1:k, drop=FALSE])*phases$sdev[1:k])

    ## start: angle in PC1/PC2 of the original PCA
    lambda <- atan2(Y[,2], Y[,1])

    ## least squares Fourier fit, slightly regularized
    fitc <- function(l) {
        X <- fourier(l, K)
        solve(crossprod(X) + diag(1e-8*nrow(X), ncol(X)), crossprod(X, Y))
    }

    g <- seq(-pi, pi, length.out=ngrid+1)[-1]
    nn <- NULL
    for ( i in seq_len(iter) ) {

        ## fit curve and evaluate on grid
        cf1 <- fitc(lambda)
        crv1 <- fourier(g, K) %*% cf1

        ## arc length parametrization of the grid, in -pi:pi
        seg <- sqrt(rowSums((crv1 - crv1[c(ngrid, 1:(ngrid-1)),])^2))
        s1 <- cumsum(seg)
        s1 <- 2*pi*s1/s1[ngrid] - pi

        ## project cells to the nearest grid point; stop if the cells
        ## collapse onto a few points of the curve
        nn1 <- nearest(Y, crv1)
        if ( length(unique(nn1$index)) < 4*(2*K+1) ) {
            warning("principal curve collapsed at iteration ", i,
                    ", keeping the previous phase")
            break
        }
        nn <- nn1
        cf <- cf1
        s <- s1
        new <- s[nn$index]

        dl <- abs(center_phase(new - lambda))
        lambda <- new
        if ( mean(dl) < tol ) break
    }
    if ( is.null(nn) ) stop("principal curve collapsed at the start")

    ## align with the current theta: direction and mean offset
    th <- phases$rotation.phase$theta
    rs <- sapply(c(1, -1), function(sg) {
        d <- th - sg*lambda
        c(R=abs(mean(exp(1i*d))), off=Arg(mean(exp(1i*d))))
    })
    sg <- c(1, -1)[which.max(rs['R',])]
    off <- unname(rs['off', which.max(rs['R',])])
    lambda <- center_phase(sg*lambda + off)

    phases$rotation.phase$lambda <- lambda
    phases$rotation.phase$curve.dist <- nn$dist
    phases$curve <- list(grid=data.frame(lambda=s, fourier(g, K) %*% cf),
                         coef=cf,
                         k=k, K=K, iter=i, sign=sg, offset=off)
    phases$processing <- c(phases$processing, "principal curve")
    phases
}

## nearest row of B for each row of A, in chunks
nearest <- function(A, B, chunk=2000) {
    idx <- dst <- numeric(nrow(A))
    bb <- rowSums(B^2)
    for ( st in seq(1, nrow(A), by=chunk) ) {
        j <- st:min(nrow(A), st+chunk-1)
        d2 <- outer(rowSums(A[j,,drop=FALSE]^2), bb, '+') -
            2*A[j,,drop=FALSE] %*% t(B)
        idx[j] <- max.col(-d2, ties.method='first')
        dst[j] <- sqrt(pmax(d2[cbind(seq_along(j), idx[j])], 0))
    }
    list(index=idx, dist=dst)
}
