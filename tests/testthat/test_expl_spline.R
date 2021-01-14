library(testthat)
library(ICSsmoothing)

test_that("cics_unif_explicit works correctly.", {
  n <- as.integer(runif(1, 3, 100))
  a <- -10
  b <- 10
  uu <- c(seq(a, b, length.out = n + 2))
  yy <- c(runif(n + 2, -20, 20))
  d <- c(runif(2, -20, 20))
  clrs <- c(sample(colours(), 2))
  exp_sp <- cics_unif_explicit(a, b, yy, d, clrs)
  fl <- rep(0, n)
  fr <- rep(0, n)
  dl <- rep(0, n)
  dr <- rep(0, n)
  ddl <- rep(0, n)
  ddr <- rep(0, n)
  for (i in 1:(length(uu) - 2)) {
    fv_l <- as.function(exp_sp$spline_polynomial[[i]])
    fv_r <- as.function(exp_sp$spline_polynomial[[i + 1]])
    fl[i] <- fv_l(uu[i + 1])
    fr[i] <- fv_r(uu[i + 1])
    poly_der_l <- deriv(exp_sp$spline_polynomial[[i]])
    der_l <- as.function(poly_der_l)
    poly_der_r <- deriv(exp_sp$spline_polynomial[[i + 1]])
    der_r <- as.function(poly_der_r)
    dl[i] <- der_l(uu[i + 1])
    dr[i] <- der_r(uu[i + 1])
    der2l <- as.function(deriv(poly_der_l))
    der2r <- as.function(deriv(poly_der_r))
    ddl[i] <- der2l(uu[i + 1])
    ddr[i] <- der2r(uu[i + 1])
  }
  y_l <- as.function(exp_sp$spline_polynomial[[1]])(uu[1])
  y_r <- as.function(exp_sp$spline_polynomial[[n + 1]])(uu[n +
    2])
  ext_d <- c(as.function(deriv(exp_sp$spline_polynomial[[1]]))(uu[1]),
    as.function(deriv(exp_sp$spline_polynomial[[n + 1]]))(uu[length(uu)]))

  tf1 <- all.equal(fl, fr,tolerance = 1e-3) && all.equal(fl, yy[2:(length(yy) -
    1)],tolerance = 1e-3)
  tf2 <- all.equal(y_l, yy[1],tolerance = 1e-3) && all.equal(y_r, yy[n + 2],tolerance = 1e-3)
  tf3 <- all.equal(dl, dr,tolerance = 1e-3) && all.equal(ddl, ddr,tolerance = 1e-3) && all.equal(ext_d,
    d,tolerance = 1e-3)
  is_spline <- tf1 && tf2 && tf3

  expect_equal(is_spline, TRUE)
})

test_that("cics_explicit works correctly.", {
  n <- as.integer(runif(1, 3, 100))
  a <- -10
  b <- 10
  uu <- c(sort(runif(n + 2, a, b)))
  yy <- c(runif(n + 2, -20, 20))
  d <- c(runif(2, -20, 20))
  clrs <- c(sample(colours(), 2))
  exp_sp <- cics_explicit(uu, yy, d, clrs)
  fl <- rep(0, n)
  fr <- rep(0, n)
  dl <- rep(0, n)
  dr <- rep(0, n)
  ddl <- rep(0, n)
  ddr <- rep(0, n)
  for (i in 1:(length(uu) - 2)) {
    fv_l <- as.function(exp_sp$spline_polynomial[[i]])
    fv_r <- as.function(exp_sp$spline_polynomial[[i + 1]])
    fl[i] <- fv_l(uu[i + 1])
    fr[i] <- fv_r(uu[i + 1])
    poly_der_l <- deriv(exp_sp$spline_polynomial[[i]])
    der_l <- as.function(poly_der_l)
    poly_der_r <- deriv(exp_sp$spline_polynomial[[i + 1]])
    der_r <- as.function(poly_der_r)
    dl[i] <- der_l(uu[i + 1])
    dr[i] <- der_r(uu[i + 1])
    der2l <- as.function(deriv(poly_der_l))
    der2r <- as.function(deriv(poly_der_r))
    ddl[i] <- der2l(uu[i + 1])
    ddr[i] <- der2r(uu[i + 1])
  }
  y_l <- as.function(exp_sp$spline_polynomial[[1]])(uu[1])
  y_r <- as.function(exp_sp$spline_polynomial[[n + 1]])(uu[n +
    2])
  ext_d <- c(as.function(deriv(exp_sp$spline_polynomial[[1]]))(uu[1]),
    as.function(deriv(exp_sp$spline_polynomial[[n + 1]]))(uu[length(uu)]))

  tf1 <- all.equal(fl, fr,tolerance = 1e-2) && all.equal(fl, yy[2:(length(yy) -
    1)],tolerance = 1e-3)
  tf2 <- all.equal(y_l, yy[1],tolerance = 1e-3) && all.equal(y_r, yy[n + 2],tolerance = 1e-3)
  tf3 <- all.equal(dl, dr,tolerance = 1e-3) && all.equal(ddl, ddr,tolerance = 1e-3) && all.equal(ext_d,
    d,tolerance = 1e-3)
  is_spline <- tf1 && tf2 && tf3

  expect_equal(is_spline, TRUE)
})
