## ----message=FALSE------------------------------------------------------------
rm(list = ls())
library(ICSsmoothing)

## ----fig.width=4, fig.height=3------------------------------------------------
yy <- c(2, 3, -2, 6)
d <- c(1, 4)
clrs <- c("blue", "red")
expl_spline <- cics_unif_explicit(-1, 2, yy, d, clrs)

## -----------------------------------------------------------------------------
expl_spline$spline_coeffs
expl_spline$spline_polynomials

## -----------------------------------------------------------------------------
expl_spline$B
expl_spline$gamma

## -----------------------------------------------------------------------------
B <- expl_spline$B
g <- expl_spline$gamma
i <- 2
k <- 3
B[,,i][k,] %*% g;

## -----------------------------------------------------------------------------
explicit_spline(expl_spline$B, expl_spline$gamma) 

## ----fig.width=6, fig.height=4------------------------------------------------
sp <- cics_unif_explicit_smooth(
  xx = CERN$x,
  yy = CERN$y,
  k = 19, # 19 or 23
  #,d = c(6,0)
  ,ylab = "pi-p"
)


## -----------------------------------------------------------------------------
sp$est_spline_polynomials[1]

## ----fig.width=6, fig.height=4------------------------------------------------
yy <- as.vector( log10(AirPassengers) )
xx <- c(1:length(yy))
k <- 24
airp_spline <- cics_unif_explicit_smooth(xx, yy, k, c("blue","red")
        ,title=paste0("Smoothing log(AirPassengers) with a ", k,"-component uniform CICS"))

## ----fig.width=6, fig.height=4------------------------------------------------
ud <- forecast_demo()

## ----fig.width=6, fig.height=4------------------------------------------------
sp <- cics_unif_explicit_smooth(CERN$x,CERN$y,19,
                                #xlab = "x",ylab = "y", 
                                plotTF =  FALSE);
sp$nodes

## ----fig.width=6, fig.height=4------------------------------------------------

uu <- c(1, 15, 26, 63, 73, 88, 103, 117, 
        132, 200, 203, 219, 258, 277)
sp <- cics_explicit_smooth(
  xx = CERN$x,
  yy = CERN$y,
  uu
  #, d = c(5.57, 0.05)
  )

