library(testthat)
library(ICSsmoothing)

test_examples()
#test_example(path = "../../man/CERN.Rd", title = path)
#test_example(path = "../../man/hermite_bf_matrix.Rd", title = path)
##W test_example(path = "../../man/nonunif_explicit_clCS.Rd", title = path)
##W test_example(path = "../../man/smooth_by_nonunif_explicit_clCS.Rd", title = path)
##W test_example(path = "../../man/smooth_by_unif_explicit_clCS.Rd", title = path)
###F test_example(path = "../../man/tridiag_inv_general.Rd", title = path)
#test_example(path = "man/tridiag_inv_unif_by_sums.Rd", title = path)
##W test_example(path = "../../man/unif_explicit_clCS.Rd", title = path)

test_that("Algorithm for inverse of general tridiagonal matrix works correctly.", {
  n<-as.integer(runif(1,min = 1,max = 100))
  T <- matrix(0,n,n)
  diag(T) <- c(runif(n,min=-100,max=100))
  if(n>1){
  indx <- seq.int(n-1)
  T[cbind(indx+1,indx)] <- c(runif(n-1,min=-100,max=100))
  T[cbind(indx,indx+1)] <- c(runif(n-1,min=-100,max=100))
  }
  if(det(T)==0){expect_error(tridiag_inv_general(T,n),"Matrix T is not invertible!")}
  else{
  expect_equal(T %*% tridiag_inv_general(T,n),diag(n))
  expect_equal(tridiag_inv_general(T,n) %*% T,diag(n))
  }
})

test_that("Algorithm for inverse of tridiagonal matrix with constant diagonals works correctly.", {
  n<-as.integer(runif(1,min = 1,max = 100))
  T <- matrix(0,n,n)
  b<-runif(1,min=-100,max=100)
  diag(T) <- rep(b,n)
  if(n>1){
    indx <- seq.int(n-1)
    a<-runif(1,min=-100,max=100)
    T[cbind(indx+1,indx)] <- rep(a,n-1)
    T[cbind(indx,indx+1)] <- rep(a,n-1)
  }
  if(n == 1 && b == 0){expect_error(tridiag_inv_unif_by_sums(n,a,b),"The matrix is not invertible.")}
  else {
    expect_equal(T %*% tridiag_inv_unif_by_sums(n,a,b),diag(n))
    expect_equal(tridiag_inv_unif_by_sums(n,a,b) %*% T,diag(n))
  }
})

