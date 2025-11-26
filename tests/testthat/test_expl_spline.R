library(testthat)
library(ICSsmoothing)

n_from_uu <- function(uu) {
  length(uu) - 2L
}

test_that("cics_unif_explicit constructs a C2 spline on a uniform grid", {
  n <- 5L
  a <- -2
  b <- 3
  
  uu <- seq(a, b, length.out = n + 2)
  
  yy <- c(-1, 0.5, 2, -0.5, 1.5, -2, 0.7) 
  stopifnot(length(yy) == n + 2)
  
  d  <- c(0.3, -0.8)  
  
  clrs <- c("black", "red")
  
  exp_sp <- cics_unif_explicit(a, b, yy, d, clrs)
  
  expect_length(exp_sp$spline_polynomial, n + 1L)
  
  fl  <- fr  <- numeric(n)
  dl  <- dr  <- numeric(n)
  ddl <- ddr <- numeric(n)
  
  for (i in 1:(length(uu) - 2L)) {
    fv_l <- as.function(exp_sp$spline_polynomial[[i]])
    fv_r <- as.function(exp_sp$spline_polynomial[[i + 1L]])
    
    fl[i] <- fv_l(uu[i + 1L])
    fr[i] <- fv_r(uu[i + 1L])
    
    poly_der_l <- deriv(exp_sp$spline_polynomial[[i]])
    poly_der_r <- deriv(exp_sp$spline_polynomial[[i + 1L]])
    
    der_l <- as.function(poly_der_l)
    der_r <- as.function(poly_der_r)
    
    dl[i] <- der_l(uu[i + 1L])
    dr[i] <- der_r(uu[i + 1L])
    
    der2l <- as.function(deriv(poly_der_l))
    der2r <- as.function(deriv(poly_der_r))
    
    ddl[i] <- der2l(uu[i + 1L])
    ddr[i] <- der2r(uu[i + 1L])
  }
  
  y_l <- as.function(exp_sp$spline_polynomial[[1L]])(uu[1L])
  y_r <- as.function(exp_sp$spline_polynomial[[n + 1L]])(uu[n + 2L])
  
  ext_d <- c(
    as.function(deriv(exp_sp$spline_polynomial[[1L]]))(uu[1L]),
    as.function(deriv(exp_sp$spline_polynomial[[n + 1L]]))(uu[length(uu)])
  )
  
  expect_equal(fl, fr, tolerance = 1e-3,
               info = "Left/right function values at internal knots differ (uniform grid).")
  expect_equal(fl, yy[2:(length(yy) - 1L)], tolerance = 1e-3,
               info = "Function values at internal knots do not match yy (uniform grid).")
  
  expect_equal(y_l, yy[1L],    tolerance = 1e-3,
               info = "Left boundary value does not match yy[1] (uniform grid).")
  expect_equal(y_r, yy[n + 2L], tolerance = 1e-3,
               info = "Right boundary value does not match yy[n + 2] (uniform grid).")
  
  expect_equal(dl,  dr,  tolerance = 1e-3,
               info = "First derivative is not continuous at internal knots (uniform grid).")
  expect_equal(ddl, ddr, tolerance = 1e-3,
               info = "Second derivative is not continuous at internal knots (uniform grid).")
  expect_equal(ext_d, d, tolerance = 1e-3,
               info = "Boundary derivatives do not match specified d (uniform grid).")
})


test_that("cics_explicit constructs a C2 spline on a non-uniform grid", {
  uu <- c(-2.0, -1.3, -0.4, 0.1, 1.0, 2.2, 3.0)
  n  <- n_from_uu(uu)  # očakávame n = 5L
  
  yy <- c(-0.7, 1.2, -1.0, 0.3, 2.1, -0.4, 0.8)
  stopifnot(length(yy) == length(uu))
  
  d  <- c(-0.2, 0.5)
  
  clrs <- c("blue", "green")
  
  exp_sp <- cics_explicit(uu, yy, d, clrs)
  
  expect_length(exp_sp$spline_polynomial, n + 1L)
  
  fl  <- fr  <- numeric(n)
  dl  <- dr  <- numeric(n)
  ddl <- ddr <- numeric(n)
  
  for (i in 1:(length(uu) - 2L)) {
    fv_l <- as.function(exp_sp$spline_polynomial[[i]])
    fv_r <- as.function(exp_sp$spline_polynomial[[i + 1L]])
    
    fl[i] <- fv_l(uu[i + 1L])
    fr[i] <- fv_r(uu[i + 1L])
    
    poly_der_l <- deriv(exp_sp$spline_polynomial[[i]])
    poly_der_r <- deriv(exp_sp$spline_polynomial[[i + 1L]])
    
    der_l <- as.function(poly_der_l)
    der_r <- as.function(poly_der_r)
    
    dl[i] <- der_l(uu[i + 1L])
    dr[i] <- der_r(uu[i + 1L])
    
    der2l <- as.function(deriv(poly_der_l))
    der2r <- as.function(deriv(poly_der_r))
    
    ddl[i] <- der2l(uu[i + 1L])
    ddr[i] <- der2r(uu[i + 1L])
  }
  
  y_l <- as.function(exp_sp$spline_polynomial[[1L]])(uu[1L])
  y_r <- as.function(exp_sp$spline_polynomial[[n + 1L]])(uu[n + 2L])
  
  ext_d <- c(
    as.function(deriv(exp_sp$spline_polynomial[[1L]]))(uu[1L]),
    as.function(deriv(exp_sp$spline_polynomial[[n + 1L]]))(uu[length(uu)])
  )
  
  expect_equal(fl, fr, tolerance = 1e-3,
               info = "Left/right function values at internal knots differ (non-uniform grid).")
  expect_equal(fl, yy[2:(length(yy) - 1L)], tolerance = 1e-3,
               info = "Function values at internal knots do not match yy (non-uniform grid).")
  
  expect_equal(y_l, yy[1L],     tolerance = 1e-3,
               info = "Left boundary value does not match yy[1] (non-uniform grid).")
  expect_equal(y_r, yy[n + 2L], tolerance = 1e-3,
               info = "Right boundary value does not match yy[n + 2] (non-uniform grid).")
  
  expect_equal(dl,  dr,  tolerance = 1e-3,
               info = "First derivative is not continuous at internal knots (non-uniform grid).")
  expect_equal(ddl, ddr, tolerance = 1e-3,
               info = "Second derivative is not continuous at internal knots (non-uniform grid).")
  expect_equal(ext_d, d, tolerance = 1e-3,
               info = "Boundary derivatives do not match specified d (non-uniform grid).")
})