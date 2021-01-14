
 ICSsmoothing
============

Library for data smoothing using interpolating splines

Basic usage examples
--------------------

Load the library if installed.
```R
library("ICSsmoothing")
```

The library consists one testing dataframe with two numeric columns `x` and `y`.

* `CERN`, where `CERN$x = [1, 2, ..., 277]`,

and five public functions 

* `cics_unif_explicit(...)` - create and plot the uniform explicit spline,

	```R
	cics_unif_explicit(
	  uumin = CERN$x[1], # start spline at this x-coordinate
	  uumax = CERN$x[10], # finish spline at this x-coordinate
	  yy = CERN$y[1:10], # y-coordinates for spline control points
	  d = c(0,0), # derivateves at points (uumin, yy[1]) and (uumax, yy[length(yy)])
	  xlab="X axis", # x-axis label
	  ylab="Y axis" # y-axis label
	)
	```
* `cics_unif_explicit_smooth(...)` - create and plot the uniform explicit spline as smoothing curve,

	```R
	cics_unif_explicit_smooth(
	  xx = CERN$x, # x-coordinates to smooth
	  yy = CERN$y, # y-coordinates to smooth
	  k = 21, # number of components of a smoothing spline
	  d = c(0, 1), # derivateves at points (uumin, yy[1]) and (uumax, yy[length(yy)])
	  xlab = "X axis", # x-axis label
	  ylab = "Y axis"
	)
	```

* `cics_explicit(...)` - create and plot the explicit spline,

	```R
	cics_explicit(
	  uu = c(1, 2.2, 3, 3.8, 7), # x-coordinates for spline control points
	  yy = CERN$y[1:5], # y-coordinates for spline control points
	  d = c(0,0), # derivateves at points (uu[1], yy[1]) and (uu[length(uu)], yy[length(yy)])
	  xlab="X axis", # x-axis label
	  ylab="Y axis" # y-axis label
	)
	```
* `cics_explicit_smooth(...)` - create and plot the explicit spline as smoothing curve,

	```R
	cics_explicit_smooth(
	  xx = CERN$x, # x-coordinates to smooth
	  yy = CERN$y, # y-coordinates to smooth
	  uu = c(1, 4, 7, 20, 41, 57, 86, 92, 101, 121, 220, 245, 261, 277), # # x-coordinates for spline control points. uu[1] == xx[1] and uu[length(uu)] == xx[length(xx)]
	  d = c(0, 1), # derivateves at points (uumin, yy[1]) and (uumax, yy[length(yy)])
	  xlab = "X axis", # x-axis label
	  ylab = "Y axis"
	)
	```
* `forecast_demo()` - demo showing a usage of an explicit spline to forecast some data. See the commented source code of the function for more details about the functionality.

Building from source
--------------------

###  Prerequisites
* [Rtools](https://support.rstudio.com/hc/en-us/articles/200486498-Package-Development-Prerequisites)
* LaTeX

### Depends on R packages
* devtools
* roxygen2
* ggplot2

### Installation
On Windows Windows run in R console (or RStudio)
```R
install.packages("Rtools", force = TRUE)
library("Rtools")
```
On Linux (Ubuntu) install first these packages
```bash
sudo apt install libxml2-dev libssl-dev libcurl4-openssl-dev libopenblas-dev r-base r-base-dev
```

In both cases then run in R console
```R
devtools::install_github("klutometis/roxygen", force = TRUE)
library(roxygen2)
install.packages("digest")
library(digest)
install.packages("polynom")
library(polynom)
install.packages("ggplot2")
library(ggplot2)
```
Open system terminal in the parent directory of the package and type
```bash
R CMD build ICSsmoothing
R CMD install ICSsmoothing_1.1.0.tar.gz
```
**OR**
Open RStudio and install this library using top menu button *Build*->*Install and Restart*.
