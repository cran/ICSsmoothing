library(polynom)
library(ggplot2)

.tridiagonal_determinant <- function(n, a, b) {
  result <- array(0, n + 1)
  result[1] = 1
  result[2] = b
  if(n > 1){
    for (i in 3:(n + 1)) {
      result[i] <- (b * result[i - 1]) - (a * a * result[i - 2])
    }
  }
  return(result)
}

.tridiagonal_inverse <- function(n, a, b, determinant) {
  Inv <- matrix(0, nrow = n, ncol = n)
  if (n == 1) {
    Inv[1, 1] <- (1 / b)
    return(Inv)
  }

  Inv[1, n] <-
    ((-1) ^ (1 + n)) * (a ^ (n - 1)) / determinant[n + 1]
  Inv[n, 1] <- Inv[1, n]
  for (j in 1:(n - 1)) {
    Inv[n - j + 1, n] <- ((-1) ^ (n + j)) * Inv[1, n] *
      (a ^ (j - n)) * determinant[n - j + 1]
    Inv[1, j] <- Inv[n - j + 1, n]
    Inv[n, n - j + 1] <- Inv[n - j + 1, n]
    Inv[j, 1] <- Inv[n - j + 1, n]
  }


  s <- 0
  for (j in (n - 1):((n + (n %% 2)) / 2)) {
    s <- s + 1
    for (i in (n + 1 - j):(n - s)) {
      if (s + 1 - j <= 0) {
        Inv[i, j] <- (Inv[i + 1, j + 1] + Inv[i + j -
                                                n, n])
        Inv[j, i] <- Inv[i, j]
        Inv[n - i + 1, n - j + 1] <- Inv[i, j]
        Inv[n - j + 1, n - i + 1] <- Inv[i, j]
      }
    }
  }
  return(Inv)
}

#' Construct inverse of a tridiagonal matrix T_n(a,b,a).
#'
#' \code{tridiag_inv_unif_by_sums} constructs inverse of a regular tridiagonal matrix \code{T}_{\code{n}}(\code{a},\code{b},\code{a})
#' with constant entries by a special algorithm using sums of matrix elements.
#'
#' @param n an order of given tridiagonal matrix.
#' @param a a value of tridiagonal matrix elements that are off-diagonal.
#' @param b a value of tridiagonal matrix diagonal elements.
#' @return The inverse of matrix \code{T}_{\code{n}}(\code{a},\code{b},\code{a}).
#' @examples
#' tridiag_inv_unif_by_sums(5, 1, 4)
#' tridiag_inv_unif_by_sums(9, 10, -1)
#' @export
tridiag_inv_unif_by_sums <- function(n, a, b) {
  if (n %% 1 != 0 || n <= 0) {
    stop("n is not a positive integer!")
  }
  if (!is.numeric(a) || !is.numeric(b)) {
    stop("a or b are not numeric!")
  }
  if (n == 1 && b == 0) {
    stop("The matrix is not invertible.")
  }

  det <- .tridiagonal_determinant(n, a, b)

  return (.tridiagonal_inverse(n, a, b, det))
}

#' Construct inverse of a general tridiagonal matrix.
#'
#' \code{tridiag_inv_general} constructs inverse of a general tridiagonal matrix \code{T} of order \code{n},
#' using Usmani's theorem.
#'
#' @param T a tridiagonal matrix.
#' @param n an order of given tridiagonal matrix.
#' @return The inverse of matrix T.
#' @examples
#' tridiag_inv_general(matrix(c(1, 4, 0, -9), 2, 2), 2)
#' tridiag_inv_general(matrix(c(1, 3, 5, -2, 0, 8, 7, 6, 6), 3, 3), 3)
#' @export
tridiag_inv_general <- function(T, n) {
  if (!is.matrix(T)) {
    stop("T is not a matrix!")
  }
  if (n == 1) {
    return(matrix((1 / T[1, 1]), 1, 1))
  }
  if (n %% 1 != 0 || n <= 0) {
    stop("n is not a positive integer!")
  }
  if (det(T) == 0) {
    stop("Matrix T is not invertible!")
  }

  a <- rep(0, n - 1)
  b <- rep(0, n)
  c <- rep(0, n - 1)
  for (i in 1:(n - 1)) {
    b[i] <- T[i, i]
    c[i] <- T[i, i + 1]
    a[i] <- T[i + 1, i]
  }
  b[n] <- T[n, n]

  theta <- rep(0, n + 1)
  phi <- rep(0, n + 1)
  theta[1] <- 1
  theta[2] <- b[1]
  for (i in 3:(n + 1)) {
    theta[i] <- b[i - 1] * theta[i - 1] - a[i - 2] * c[i -
                                                         2] * theta[i - 2]
  }
  phi[n + 1] <- 1
  phi[n] <- b[n]
  for (i in (n - 1):1) {
    phi[i] <- b[i] * phi[i + 1] - a[i] * c[i] * phi[i +
                                                      2]
  }
  Tau <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      if (i < j) {
        Tau[i, j] <-
          ((-1) ^ (i + j)) * prod(c[i:(j - 1)]) * theta[i] * phi[j + 1] / theta[n + 1]
      } else if (i > j) {
        Tau[i, j] <-
          ((-1) ^ (i + j)) * prod(a[j:(i - 1)]) * theta[j] * phi[i + 1] / theta[n + 1]
      } else {
        Tau[i, j] <- theta[i] * phi[i + 1] / theta[n + 1]
      }
    }
  }

  return(Tau)
}

#' Construct 4 Hermite basis functions.
#'
#' \code{hermite_bf_matrix} constructs matrix of Hermite basis functions' coefficients on
#' [\code{u},\code{v}], that is the matrix of 4 cubic polynomials' coefficients of
#' one-component Hermite cubic spline.
#'
#' @param u a left border of interval [\code{u},\code{v}].
#' @param v a right border of interval [\code{u},\code{v}], \code{u}\eqn{\le}\code{v}.
#' @return The matrix of 4 Hermite basis functions' coefficients.
#' @examples
#' hermite_bf_matrix(0,1)
#' hermite_bf_matrix(-2,3)
#' @export
hermite_bf_matrix <- function(u, v) {
  if (!is.numeric(u) || !is.numeric(v)) {
    stop("u or v are not numeric!")
  }
  if (u >= v) {
    stop("Badly chosen nodes!")
  }
  bf_matrix <- matrix(0, 4, 4)
  h <- v - u
  hh <- v + u
  bf_matrix[1, 1] <- -3 * u * v ^ 2 + v ^ 3
  bf_matrix[1, 2] <- 3 * v * u ^ 2 - u ^ 3
  bf_matrix[1, 3] <- -u * v ^ 2 * h
  bf_matrix[1, 4] <- -u ^ 2 * v * h
  bf_matrix[2, 1] <- 6 * u * v
  bf_matrix[2, 2] <- -6 * u * v
  bf_matrix[2, 3] <- (2 * u * v + v ^ 2) * h
  bf_matrix[2, 4] <- (u ^ 2 + 2 * u * v) * h
  bf_matrix[3, 1] <- -3 * (hh)
  bf_matrix[3, 2] <- 3 * (hh)
  bf_matrix[3, 3] <- (-2 * v - u) * h
  bf_matrix[3, 4] <- (-2 * u - v) * h
  bf_matrix[4, 1] <- 2
  bf_matrix[4, 2] <- -2
  bf_matrix[4, 3] <- h
  bf_matrix[4, 4] <- h
  return(bf_matrix / (h ^ 3))
}

.uniform_auxiliary_matrix_a <- function(n, h, Tau) {
  A <- matrix(0, n, (n + 4))
  if (n == 1) {
    A[1, 1] <- -3 / (4 * h)
    A[1, 2] <- 0
    A[1, 3] <- 3 / (4 * h)
    A[1, 4] <- -1 / 4
    A[1, 5] <- -1 / 4
  } else {
    for (i in 1:n) {
      for (j in 1:(n + 4)) {
        if (j == 1 || j == 2) {
          A[i, j] <- -Tau[i, j] * (3 / h)
        } else if (j > 2 && j < (n + 1)) {
          A[i, j] <- (Tau[i, j - 2] - Tau[i, j]) * (3 / h)
        } else if (j == n + 1 || j == n + 2) {
          A[i, j] <- Tau[i, j - 2] * (3 / h)
        } else if (j == n + 3) {
          A[i, j] <- -Tau[i, 1]
        } else {
          A[i, j] <- -Tau[i, n]
        }
      }
    }
  }
  return (A)
}

.basis_matrices <- function(n, uu) {
  vec_mat <- array(0, c(4, 4, (n + 1)))
  for (i in 1:(n + 1)) {
    m <- hermite_bf_matrix(uu[i], uu[i + 1])
    vec_mat[, , i] <- m
  }
  return (vec_mat)
}

.basis_functions <- function(n, uu) {
  vec_mat <- .basis_matrices(n, uu)
  BF <- array(dim = c(n + 1, 4, 4))
  for (i in 1:(n + 1)) {
    BF[i, 1,] <- vec_mat[, , i][, 1]
    BF[i, 2,] <- vec_mat[, , i][, 2]
    BF[i, 3,] <- vec_mat[, , i][, 3]
    BF[i, 4,] <- vec_mat[, , i][, 4]
  }
  return (BF)
}

.uniform_auxiliary_matrix_b <- function(n, A, BF) {
  B <- array(dim = c(n + 1, n + 4, 4))
  B[1, 1,] <- BF[1, 1,] + (A[1, 1] * BF[1, 4,])
  B[1, 2,] <- BF[1, 2,] + (A[1, 2] * BF[1, 4,])
  for (i in 3:(n + 2)) {
    B[1, i,] <- A[1, i] * BF[1, 4,]
  }
  B[1, n + 3,] <- BF[1, 3,] + (A[1, n + 3] * BF[1, 4,])
  B[1, n + 4,] <- A[1, n + 4] * BF[1, 4,]
  if (n >= 2) {
    for (i in 2:n) {
      for (j in 1:(n + 4)) {
        if (i == j) {
          B[i, j,] <- BF[i, 1,] + A[i - 1, j] * BF[i,
                                                   3,] + A[i, j] * BF[i, 4,]
        } else if (i == j - 1) {
          B[i, j,] <- BF[i, 2,] + A[i - 1, j] * BF[i,
                                                   3,] + A[i, j] * BF[i, 4,]
        } else {
          B[i, j,] <- A[i - 1, j] * BF[i, 3,] + A[i,
                                                  j] * BF[i, 4,]
        }
      }
    }
  }
  for (i in 1:n) {
    B[n + 1, i,] <- A[n, i] * BF[n + 1, 3,]
  }
  B[n + 1, n + 1,] <- BF[n + 1, 1,] + A[n, n + 1] * BF[n +
                                                         1, 3,]
  B[n + 1, n + 2,] <- BF[n + 1, 2,] + A[n, n + 2] * BF[n +
                                                         1, 3,]
  B[n + 1, n + 3,] <- A[n, n + 3] * BF[n + 1, 3,]
  B[n + 1, n + 4,] <- BF[n + 1, 4,] + A[n, n + 4] * BF[n +
                                                         1, 3,]
  return(B)
}

.explicit_spline <- function(n, B, gamma) {
  expl_spline <- array(0, c(n + 1, 1, 4))
  for (i in 1:(n + 1)) {
    pom <- array(0, c(1, 1, 4))
    for (j in 1:(n + 4)) {
      pom <- pom + gamma[j] * B[i, j,]
    }
    expl_spline[i, 1,] <- pom
  }
  return (expl_spline)
}

.create_polynomial <- function(n, uu, expl_spline, clrs, pp) {
  pL <- list(polynomial(1))
  for (j in 2:(n + 1))
    pL[[j]] <- 1 + pL[[j - 1]]
  class(pL) <- "polylist"

  for (i in 1:(n + 1)) {
    pL[[i]] <- polynomial(coef = c(expl_spline[i, 1,]))
    frb = ifelse(i %% 2 == 0, clrs[2], clrs[1])

    fnc=as.function(pL[[i]])
    pp = pp + stat_function(fun = fnc, color=frb, xlim = c(uu[i], uu[i + 1]))
    #lines(pL[[i]], xlim = c(uu[i], uu[i + 1]), col = frb)
  }
  return (list(pL = pL, pp = pp))
}

#' Construct the explicit form of uniform clamped interpolating cubic spline (UcICS).
#'
#' \code{cics_unif_explicit} constructs the explicit form of uniform clamped interpolating cubic spline
#' (via Hermite cubic spline) for nodes \code{uu}, function values \code{yy} and exterior-node derivatives
#' \code{d}.
#'
#' @param uumin a starting node.
#' @param uumax an ending node.
#' @param yy a vector of function values pertaining to nodes in \code{uu}.
#' @param d a vector of two values of derivative, in the first and the last node of \code{uu}.
#' @param clrs a vector (optional parameter) of colours that are used alternately to plot the graph of spline's components.
#' @param xlab a title (optional parameter) for the \code{x} axis.
#' @param ylab a title (optional parameter) for the \code{y} axis.
#' @param title a title (optional parameter) for the plot.
#' @return A list of spline components
#' \item{spline_coeffs}{matrix, whose \code{i}-th row contains coefficients of uniform ICS's \code{i}-th component.}
#' \item{spline_polynomials}{list of UcICS's components string representations.}
#' \item{B}{\code{4}-element array of \code{(n+1)x(n+4)} matrices, whereas element in \code{i}-th row
#' and \code{j}-th column of \code{l}-th matrix contains coefficient by \code{x^{l-1}} of cubic polynomial that is in \code{i}-th row
#' and \code{j}-th column of matrix \code{B} from spline's explicit form \deqn{S=B.\gamma.}}
#' \item{gamma}{\eqn{\gamma=} vector of spline coefficients - function values and exterior-node derivatives that
#'  takes part in the explicit form \eqn{S=B.\gamma}.}
#' \item{aux_BF}{A basis function of the spline}
#' \item{aux_tridiag_inverse}{An inverse of the tridiagonal matrix used for spline derivatives construction}
#' @examples
#' yy <- c(4, 5, 2, 1.8);
#' sp <- cics_unif_explicit(0, 6, yy, c(2, 0.9))
#' sp$spline_polynomials
#' ### <~~>
#' ### Spline components' coefficients
#' explicit_spline(sp$B, sp$gamma)
#' sp$spline_coeffs == .Last.value
#' @export
#' @importFrom grDevices colours
#' @import ggplot2
#' @import polynom
cics_unif_explicit <-
  function(uumin,
           uumax,
           yy,
           d,
           clrs = c('blue', 'red'),
           xlab = NULL,
           ylab = NULL,
           title = NULL) {
    if (uumin >= uumax) {
      stop("uumin must be smaller than uumax")
    }

    if (length(yy) <= 2) {
      stop("There are not at least 3 nodes.")
    }
    if (is.null(title)) {
      title = paste0(length(yy)-1,"-component uniform interpolation spline")
    }
    if (is.null(xlab)) {
      xlab <- "x"
    }
    if (is.null(ylab)) {
      ylab <- deparse(substitute(yy))
    }

    x <- y <- NULL
    n <- length(yy) - 2
    dist <- (uumax - uumin) / (n + 1)
    if (length(d) != 2) {
      stop("d isn't of length 2.")
    }
    if (!all(clrs %in% colours())) {
      stop("Not every string in clrs represents a colour!")
    }

    gam <- c(yy, d)
    Tau <- tridiag_inv_unif_by_sums(n, 1, 4)
    h <- dist
    uu <- c(seq(uumin, uumax, length.out = length(yy)))
    A <- .uniform_auxiliary_matrix_a(n, h, Tau)
    BF <- .basis_functions(n, uu)
    B <- .uniform_auxiliary_matrix_b(n, A, BF)
    expl_spline <- .explicit_spline(n, B, gam)

    df <- data.frame(x = uu, y = yy)#do.call(rbind, Map(data.frame, x=uu, y=yy)),
    pp = ggplot(
      df,
      aes(x, y)
    ) + xlab(xlab) + ylab(ylab) + ggtitle(title)

    # aug 31
    pp = pp + geom_point(data = df)
    #pp = pp + geom_point(shape=1, data = df)
    pol <- .create_polynomial(n, uu, expl_spline, clrs, pp)

    pL <- pol$pL
    pp <- pol$pp
    print(pp)

    list(
      spline_coeffs = expl_spline[, ,],
      spline_polynomials = pL,
      B = B,
      gamma = gam,
      aux_BF = BF,
      aux_tridiag_inverse = Tau
    )
  }

.coordinate_indices <- function(lngth, k, uu, xx) {
  xx_idx <- c(rep(0, lngth))

  for (i in 1:lngth) {
    for (j in 1:(k - 1)) {
      if (uu[j] <= xx[i] && xx[i] < uu[j + 1]) {
        xx_idx[i] <- j
      }
      if (uu[k] <= xx[i] && xx[i] <= uu[k + 1]) {
        xx_idx[i] <- k
      }
    }
  }
  return (xx_idx)
}

.auxiliary_estimation_matrix <- function(lngth, k, xx, B, xx_idx) {
  M <- matrix(0, lngth, k + 3)
  for (i in 1:lngth) {
    for (j in 1:(k + 3)) {
      sup_poly <- B[xx_idx[i], j,]
      M[i, j] <- sup_poly[1] + sup_poly[2] * xx[i] + sup_poly[3] *
        (xx[i]) ^ 2 + sup_poly[4] * (xx[i]) ^ 3
    }
  }
  return (M)
}

.estimate_gamma <- function(k, d, yy, M) {
  if (missing(d)) {
    MT <- t(M)
    XX <- MT %*% M
    XY <- MT %*% yy
    ##est_gam <- t(solve(XX)) %*% XY
    #est_gam <- t(solve(XX, XY))
    est_gam <- qr.solve(M, yy)
  } else {
    if (length(d) != 2) {
      stop("d isn't a numeric vector of length 2.")
    }
    M1 <- M[, 1:(k + 1)]
    M2 <- M[, (k + 2):(k + 3)]
    yy_new <- yy - (M2 %*% d)
    M1T <- t(M1)
    XX1 <- M1T %*% M1
    XY1 <- M1T %*% yy_new
    ##est_gam_fv <- t(solve(XX1)) %*% XY1
    #est_gam_fv <- t(solve(XX1, XY1))
    est_gam_fv <- qr.solve(M1, yy_new)
    est_gam <- c(est_gam_fv, d)
  }
  return(est_gam)
}

#' Smooth given data set by k-component uniform clamped interpolating spline (UcICS).
#'
#' \code{cics_unif_explicit_smooth} constructs the uniform clamped interpolating spline with \code{k} components that smoothes
#'  given data set \code{{(xx[i],yy[i]), i=1..length(xx)}}.
#'
#' @param xx a vector of data set's \code{x}-coordinates (that are in increasing order).
#' @param yy a vector of datanvidi set's \code{y}-coordinates.
#' @param k a chosen number of components of smoothing UcICS (integer \eqn{\ge 2}).
#' @param clrs a vector of colours that are used alternately to plot the graph of spline's components.
#' @param d a vector (optional parameter) that contains two values of derivative, in the first and the last
#'  computed node. If missing, values of derivative are estimated by given linear regression model. If present, their
#'  contribution is removed from linear model and only function values are estimated.
#' @param xlab a title (optional parameter) for the \code{x} axis.
#' @param ylab a title (optional parameter) for the \code{y} axis.
#' @param title a title (optional parameter) for the plot.
#' @param plotTF a boolean value (optional parameter), if TRUE then plot.
#' @return a list with components
#' \item{nodes}{vector of equidistant nodes, based on which we construct the smoothing spline.}
#' \item{est_spline_coeffs}{\code{4}-element array of \code{(k)x(k+3)} matrices, whereas element in \code{i}-th row and \code{j}-th
#'  of \code{l}-th matrix contains coefficient by \code{x^{l-1}} of cubic polynomial, which is in \code{i}-th row and \code{j}-th column  of matrix
#' \code{B} from smoothing spline's explicit form \deqn{S=B.\gamma.}}
#' \item{est_spline_polynomials}{list of string representations of smoothing UcICS.}
#' \item{est_gamma}{vector of estimated smoothing spline's coefficients (function values and exterior-node derivatives).}
#' \item{aux_BF}{A basis function of the spline}
#' \item{aux_tridiag_inverse}{An inverse of the tridiagonal matrix used for spline derivatives construction}
#' \item{aux_M}{An estimation matrix used to compute \code{est_gamma}}
#' @examples
#'
#' cp <- cics_unif_explicit_smooth(
#' xx = CERN$x,
#'   yy = CERN$y,
#'   k = 19, #23,
#'   d = c(1, 0),
#'   xlab = "X axis",
#'   ylab = "Y axis"
#' )
#' @export
#' @importFrom grDevices colours
#' @import ggplot2
#' @import polynom
cics_unif_explicit_smooth <-
  function(xx,
           yy,
           k,
           clrs = c('blue', 'red'),
           d,
           xlab = NULL,
           ylab = NULL,
           title = NULL,
           plotTF=TRUE) {
    if (is.unsorted(xx)) {
      stop("x-coordinates of measurements are not in increasing order!")
    }
    if (length(xx) != length(yy)) {
      stop("Lengths of x-coordinate and y-coordinate sequences differ!")
    }
    if (!all(clrs %in% colours())) {
      stop("Not every string in clrs represents a colour!")
    }
    if (k < 2) {
      stop("k is either not greater or equal to 2.")
    }
    if (is.null(title)) {
      title = paste0("Smoothing with a ", k,"-component uniform interpolation spline")
    }
    if (is.null(xlab)) {
      xlab <- deparse(substitute(xx))
    }
    if (is.null(ylab)) {
      ylab <- deparse(substitute(yy))
    }

    x <- y <- NULL
    lngth <- length(xx)
    n <- k - 1
    uu <- c(seq(xx[1], xx[length(xx)], length.out = k + 1))
    h <- (xx[length(xx)] - xx[1]) / k
    Tau <- tridiag_inv_unif_by_sums(n, 1, 4)

    A <- .uniform_auxiliary_matrix_a(n, h, Tau)
    BF <- .basis_functions(n, uu)
    B <- .uniform_auxiliary_matrix_b(n, A, BF)

    xx_idx <- .coordinate_indices(lngth, k, uu, xx)

    M <- .auxiliary_estimation_matrix(lngth, k, xx, B, xx_idx)
    est_gam <- .estimate_gamma(k, d, yy, M)
    expl_spline <- .explicit_spline(n, B, est_gam)

    df <- data.frame(x = xx, y = yy)#do.call(rbind, Map(data.frame, x=xx, y=yy)),

    pp = ggplot(
      df,
      aes(x, y)
    ) + xlab(xlab) + ylab(ylab) + ggtitle(title)

    # aug 31
    #pp = pp + geom_point(data = df, size=1, alpha = 0.2)
    pp = pp + geom_point(data = df, size=1, alpha = 0.2, shape=1)
    #pp = pp + geom_point(data = df, size=1, alpha = 0.2, color='red')
    #pp = pp + geom_point(shape = 0,color='gray',data = df, size=1, alpha = 1)
    pol <- .create_polynomial(n, uu, expl_spline, clrs, pp)

    pL <- pol$pL
    pp <- pol$pp

    if (plotTF) {
      print(pp)
    }

    list(
      pp = pp,
      nodes = uu,
      est_spline_coeffs = expl_spline[, ,],
      est_spline_polynomials = pL,
      est_gamma = est_gam,
      B = B,
      aux_BF = BF,
      aux_tridiag_inverse = Tau,
      aux_M = M
    )
  } # END cics_unif_explicit_smooth

.nonuniform_tridiagonal_inverse_matrix <- function(n, hh) {
  i <- NULL
  T <- matrix(0, n, n)
  if (n > 1) {
    for (i in 1:(n - 1)) {
      T[i, i] <- (2 / hh[i]) + (2 / hh[i + 1])
      T[i, i + 1] <- (1 / hh[i + 1])
      T[i + 1, i] <- (1 / hh[i + 1])
    }
  }
  T[n, n] <- (2 / hh[n]) + (2 / hh[n + 1])
  Tau <- tridiag_inv_general(T, n)
  return (Tau)
}

.nonuniform_auxiliary_matrix_a <- function(n, hh, Tau) {
  i <- j <- NULL
  A <- matrix(0, n, n + 4)
  if (n == 1) {
    A[1, 1] <- -3 * hh[2] / (2 * hh[1] * (hh[1] + hh[2]))
    A[1, 2] <- 3 * (hh[2] - hh[1]) / (2 * hh[1] * hh[2])
    A[1, 3] <- 3 * hh[1] / (2 * hh[2] * (hh[1] + hh[2]))
    A[1, 4] <- -hh[2] / (2 * (hh[1] + hh[2]))
    A[1, 5] <- -hh[1] / (2 * (hh[1] + hh[2]))
  } else {
    for (i in 1:n) {
      for (j in 1:(n + 4)) {
        if (j == 1) {
          A[i, j] <- (-3 * Tau[i, 1] / (hh[1] ^ 2))
        } else if (j == 2) {
          A[i, j] <- (3 * Tau[i, 1] / hh[1] ^ 2) - (3 * Tau[i,
                                                            1] / hh[2] ^ 2) - (3 * Tau[i, 2] / hh[2] ^ 2)
        } else if (j > 2 && j < n + 1) {
          A[i, j] <- (3 * Tau[i, j - 2] / hh[j - 1] ^ 2) +
            (3 * Tau[i, j - 1] / hh[j - 1] ^ 2) - (3 * Tau[i,
                                                           j - 1] / hh[j] ^ 2) - (3 * Tau[i, j] / hh[j] ^ 2)
        } else if (j == n + 1) {
          A[i, j] <- (3 * Tau[i, n - 1] / hh[n] ^ 2) + (3 *
                                                          Tau[i, n] / hh[n] ^ 2) - (3 * Tau[i, n] / hh[n +
                                                                                                         1] ^ 2)
        } else if (j == n + 2) {
          A[i, j] <- (3 * Tau[i, n] / hh[n + 1] ^ 2)
        } else if (j == n + 3) {
          A[i, j] <- (-Tau[i, 1] / hh[1])
        } else {
          A[i, j] <- (-Tau[i, n] / hh[n + 1])
        }
      }
    }
  }
  return (A)
}

#' Construct the explicit form of non-uniform clamped interpolating cubic spline (NcICS).
#'
#' \code{cics_explicit} constructs the explicit form of non-uniform clamped interpolating cubic spline
#' (via Hermite cubic spline) for nodes \code{uu}, function values \code{yy} and exterior-node derivatives
#' \code{d}.
#'
#' @param uu a vector of arbitrary nodes (ordered ascendingly), with magnitude \code{n+2}, \code{n}\eqn{\ge}\code{1}.
#' @param yy a vector of function values pertaining to nodes in \code{uu}.
#' @param d a vector of two values of derivative, in the first and the last node of \code{uu}.
#' @param clrs a vector of colours that are used alternately to plot the graph of spline's components.
#' @param xlab a title (optional parameter) for the \code{x} axis.
#' @param ylab a title (optional parameter) for the \code{y} axis.
#' @param title a title (optional parameter) for the plot.
#' @return a list with components
#' \item{spline_coeffs}{matrix, whose \code{i}-th row contains coefficients of non-uniform ICS's \code{i}-th component.}
#' \item{spline_polynomials}{list of NcICS's components string representations.}
#' \item{B}{\code{4}-element array of \code{(n+1)x(n+4)} matrices, whereas element in \code{i}-th row
#' and \code{j}-th column of \code{l}-th matrix contains coefficient by \code{x^{l-1}} of cubic polynomial that is in \code{i}-th row
#' and \code{j}-th column of matrix \code{B} from spline's explicit form \deqn{S=B.\gamma.}}
#' \item{gamma}{\eqn{\gamma=} vector of spline coefficients - function values and exterior-node derivatives that
#'  takes part in the explicit form \eqn{S=B.\gamma}.}
#' \item{aux_BF}{A basis function of the spline}
#' \item{aux_tridiag_inverse}{An inverse of the tridiagonal matrix used for spline derivatives construction}
#' @examples
#' cics_explicit(
#'  uu = c(1, 2.2, 3, 3.8, 7),
#'  CERN$y[1:5],
#'  d=c(0,-2),
#'  xlab="X axis",
#'  ylab="Y axis"
#' )
#'
#' uu <- c(0, 1, 4, 6);
#' yy <- c(4, 5, 2, 1.8);
#' sp <- cics_explicit(uu, yy, c(1,0))
#' sp$spline_polynomials
#' ### <~~>
#' ### Spline components' coefficients
#' explicit_spline(sp$B, sp$gamma)
#' sp$spline_coeffs == .Last.value
#' @export
#' @importFrom grDevices colours
#' @import ggplot2
#' @import polynom
cics_explicit <-
  function(uu,
           yy,
           d,
           clrs = c('blue', 'red'),
           xlab = NULL,
           ylab = NULL,
           title = NULL) {
    if (is.unsorted(uu)) {
      stop("nodes are not in increasing order!")
    }
    if (length(uu) != length(yy)) {
      stop("Lengths of node sequence and
                  sequence of function values differ!")
    }
    if (length(uu) <= 2) {
      stop("There should be at least 3 nodes.")
    }
    if (length(d) != 2) {
      stop("d isn't of length 2.")
    }
    if (!all(clrs %in% colours())) {
      stop("Not every string in clrs represents a colour!")
    }
    if (is.null(title)) {
      title = paste0(length(uu)-1,"-component non-uniform interpolation spline")
    }
    if (is.null(xlab)) {
      xlab <- deparse(substitute(uu))
    }
    if (is.null(ylab)) {
      ylab <- deparse(substitute(yy))
    }

    x <- y <- NULL
    n <- length(uu) - 2
    hh <- c(uu[2:length(uu)] - uu[1:(length(uu) - 1)])
    gam <- c(yy, d)
    Tau <- .nonuniform_tridiagonal_inverse_matrix(n, hh)

    A <- .nonuniform_auxiliary_matrix_a(n, hh, Tau)
    vec_mat <- .basis_matrices(n, uu)

    BF <- .basis_functions(n, uu)
    B <- .uniform_auxiliary_matrix_b(n, A, BF)
    expl_spline <- .explicit_spline(n, B, gam)

    df <- data.frame(x = uu, y = yy)#do.call(rbind, Map(data.frame, x=uu, y=yy)),
    pp = ggplot(
      df,
      aes(x, y)
    ) + xlab(xlab) + ylab(ylab) + ggtitle(title)


    pp = pp + geom_point(data = df)
    pol <- .create_polynomial(n, uu, expl_spline, clrs, pp)

    pL <- pol$pL
    pp <- pol$pp
    print(pp)

    list(
      spline_coeffs = expl_spline[, ,],
      spline_polynomials = pL,
      B = B,
      gamma = gam,
      aux_BF = BF,
      aux_tridiag_inverse = Tau
    )
  }

#' Smooth given data set by k-component non-uniform clamped interpolating spline (NcICS).
#'
#' \code{cics_explicit_smooth} constructs the non-uniform clamped interpolating spline with \code{k}
#'  components that smoothes
#'  given data set \code{{(xx[i],yy[i]), i=1..length(xx)}}.
#'
#' @param xx a vector of data set's \code{x}-coordinates (that are in increasing order).
#' @param yy a vector of data set's \code{y}-coordinates.
#' @param uu a vector of arbitrary nodes, based on which we construct the smoothing spline. uu[1] and uu[length(uu)] must be equal to xx[1] and xx[length(xx)], respectively.
#' @param clrs a vector of colours that are used alternately to plot the graph of spline's components.
#' @param d a vector (optional parameter) that contains two values of derivative, in the first and the last
#'  node from \code{uu}. If missing, values of derivative are estimated by given linear regression model. If present, their
#'  contribution is removed from linear model and only function values are estimated.
#' @param xlab a title (optional parameter) for the \code{x} axis.
#' @param ylab a title (optional parameter) for the \code{y} axis.
#' @param title a title (optional parameter) for the plot.
#' @return a list with components
#' \item{est_spline_coeffs}{\code{4}-element array of \code{(k)x(k+3)} matrices, whereas element in \code{i}-th row and \code{j}-th
#'  of \code{l}-th matrix contains coefficient by \code{x^{l-1}} of cubic polynomial, which is in \code{i}-th row and \code{j}-th column  of matrix
#' \code{B} from smoothing spline's explicit form \deqn{S=B.\gamma.}}
#' \item{est_spline_polynomials}{list of string representations of smoothing NcICS.}
#' \item{est_gamma}{vector of estimated smoothing spline's coefficients (function values and exterior-node derivatives).}
#' \item{aux_BF}{A basis function of the spline}
#' \item{aux_tridiag_inverse}{An inverse of the tridiagonal matrix used for spline derivatives construction}
#' \item{aux_M}{An estimation matrix used to compute \code{est_gamma}}
#' @examples
#'
#' cics_explicit_smooth(
#' xx = CERN$x,
#' yy = CERN$y,
#' d = c(0, 1),
#' uu = c(1, sort(runif(20,1,277)), 277),
#' xlab = "X axis",
#' ylab = "Y axis"
#' )
#'
#' yy <- c(1, 2, 3, 4, 3, 2, 2, 3, 5, 6, 7, 6, 5, 5, 4, 3, 2, 1, 0)
#' xx <- c(1:length(yy))
#' uu <- c(1,7,10,19)
#' sp <- cics_explicit_smooth(xx,yy,uu)
#' ### We can change the derivatives at the end nodes:
#' sp <- cics_explicit_smooth(xx,yy, uu, d=c(3,-7/10))
#'
#' ### CERN:
#' uu <- c(1, 15, 26, 63, 73, 88, 103, 117, 132, 200, 203, 219, 258, 277)
#' sp <- cics_explicit_smooth(
#'   xx = CERN$x,
#'   yy = CERN$y,
#'   d = c(1, 0),
#'   uu
#' )
#' @export
#' @importFrom grDevices colours
#' @import ggplot2
#' @import polynom
cics_explicit_smooth <-
  function(xx,
           yy,
           uu,
           clrs = c('blue', 'red'),
           d,
           xlab = NULL,
           ylab = NULL,
           title = NULL) {
    k <- length(uu) - 1
    if (is.unsorted(xx)) {
      stop("x-coordinates of measurements are not in increasing order!")
    }
    if (is.unsorted(uu)) {
      stop("nodes uu are not in increasing order!")
    }
    if (length(xx) != length(yy)) {
      stop("Lengths of x-coordinate and y-coordinate sequences differ!")
    }
    if (!all(clrs %in% colours())) {
      stop("Not every string in clrs represents a colour!")
    }
    if (k < 2) {
      stop("k is not greater or equal to 2.")
    }
    if (is.null(title)) {
      title = paste0("Smoothing with a ", k,"-component non-uniform interpolation spline")
    }
    if (is.null(xlab)) {
      xlab <- deparse(substitute(xx))
    }
    if (is.null(ylab)) {
      ylab <- deparse(substitute(yy))
    }

    x <- y <- NULL
    lngth <- length(xx)
    n <- length(uu) - 2
    hh <- c(uu[2:length(uu)] - uu[1:(length(uu) - 1)])
    Tau <- .nonuniform_tridiagonal_inverse_matrix(n, hh)

    A <- .nonuniform_auxiliary_matrix_a(n, hh, Tau)

    vec_mat <- .basis_matrices(n, uu)
    BF <- .basis_functions(n, uu)
    B <- .uniform_auxiliary_matrix_b(n, A, BF)

    xx_idx <- .coordinate_indices(lngth, k, uu, xx)

    M <- .auxiliary_estimation_matrix(lngth, k, xx, B, xx_idx)
    est_gam <- .estimate_gamma(k, d, yy, M)
    expl_spline <- .explicit_spline(n, B, est_gam)

    df <- data.frame(x = xx, y = yy)#do.call(rbind, Map(data.frame, x=xx, y=yy)),

    pp = ggplot(
      df,
      aes(x, y)
    ) + xlab(xlab) + ylab(ylab) + ggtitle(title)

    # aug 31
    #pp = pp + geom_point(data = df, size=1, alpha = 0.2)
    pp = pp + geom_point(data = df, size=1, alpha = 0.2, shape=1)
    pol <- .create_polynomial(n, uu, expl_spline, clrs, pp)

    pL <- pol$pL
    pp <- pol$pp

    print(pp)

    list(
      est_spline_coeffs = expl_spline[, ,],
      est_spline_polynomials = pL,
      est_gamma = est_gam,
      aux_BF = BF,
      aux_A = A,
      aux_tridiag_inverse = Tau,
      aux_M = M
    )
  }

#' The function computes the coefficients of the cubic polynomials as spline components of the clamped interpolating cubic spline of class \code{C^2} in its explicit form \code{S=B * gamma}.
#'
#' @param B a \code{4}-element array of \code{(n+1)x(n+4)} matrices, whereas element in \code{i}-th row and \code{j}-th column of \code{l}-th matrix contains coefficient by \code{x^{l-1}} of cubic polynomial that is in \code{i}-th row and \code{j}-th column of matrix \code{B} from spline's explicit form \code{S=B.gamma}.
#' @param gamma a vector of spline coefficients - function values and exterior-node derivatives that takes part in the explicit form \eqn{S=B.\gamma}.
#' @return a matrix with four columns, whose \code{i}-th row contains the coefficients of the splines's \code{i}-th component.
#' @examples
#' # See functions cics_explicit, cics_unif_explicit and the vignette.
#' @export
explicit_spline <- function(B, gamma) {
  n=length(gamma)-4
  expl_spline <- array(0, c(n + 1, 4))
  for (i in 1:(n + 1)) {
    pom <- array(0, c(1, 4))
    for (j in 1:(n + 4)) {
      pom <- pom + gamma[j] * B[i, j,];
    }
    expl_spline[i,] <- pom
  }
  return (expl_spline)
}
