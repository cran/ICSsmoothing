utils::globalVariables(c("AirPassengers","uv","yv"))

#' Forecasting demo using cics_unif_explicit_smooth.
#'
#' @return a forecast result
#' @examples
#'
#' # Plots as well as the process of computation of future derivatives and values using extrapolation.
#' ud <- forecast_demo()
#'
#' @export
#' @importFrom grDevices colours
#' @importFrom graphics plot
#' @importFrom graphics lines
#' @importFrom graphics arrows
#' @importFrom graphics par
#' @importFrom graphics points
#' @importFrom graphics title
#' @importFrom stats deriv
#' @importFrom stats lm
#' @importFrom stats predict
#' @import polynom
#' @import ggplot2
forecast_demo <- function(){


  Plot1CompHerSp = function(pp,u1,u2,n,y1,y2,d1,d2,frb){
    uuu=seq(u1,u2,by=(u2-u1)/n)
    yyy=THP_HerSp(uuu,u1,u2,y1,y2,d1,d2);
    #lines(uuu,yyy,col=frb)
    p2 <- pp + geom_line(aes(uF,yF), data.frame(uF=uuu,yF=yyy), , color=frb)
    return(p2)
    #  return(pp)
  }
  THP_HerSp = function(x, v1,v2,f1,f2,D1,D2){
    ### From THP_TaylorHermitePolynomials_LIBR_4g.mw
    f1*(x-v2)^2*(1-(2*(x-v1))/(v1-v2))/(v1-v2)^2+f2*(x-v1)^2*(1-(2*(x-v2))/(v2-v1))/(v2-v1)^2+D1*(x-v2)^2*(x-v1)/(v1-v2)^2+D2*(x-v1)^2*(x-v2)/(v2-v1)^2
  }
  forecast_demo0 = function(){
    cex0 <- 1.4

    forcastD1234 = function(uu,se,sp, plotAllTF=TRUE, colors, ...){
      nu=length(uu);
      ne=length(se);
      ns=length(sp); nu; ne; ns
      ##points(uu,ss$est[1:(length(ss$est)-2)],col='magenta')
      u1=uu[nu-1]; u2=uu[nu];
      h=u2-u1; u3=u2+h; u4=u3+h;
      u1;u2;u3;u4;
      ### We need for the last component
      ### 4 parameters y1,y2,d1,d2
      ### Three of them are given naturally by:
      #y1=se[ne-3]
      #y2=se[ne-2]
      #y1;y2; #d1
      d1=se[ne-1]
      d2_=se[ne]
      # - however d2_ is unusable as we shall see

      ### 1/4 Correction of the last derivative
      ### a) Plot derivatives of spline components:
      # plot(0, pch='',xlim=c(uu[1],u4), ylim=c(-1,1))
      # title('Spline derivatives \n and point derivatives')
      pL = deriv(sp); pL
      ### b) Plot point derivatives
      dd=rep(1,nu);
      ### "predict" computes the polynomial jaj at uu[.]
      for (i in 1:ns){ dd[i]=predict(pL[[i]],uu[i]);}
      dd[ns+1]=se[ne]
      # <=>
      #predict(pL[[ns]],uu[ns+1]);
      # ??????
      dd[ns+1]; dd[ns]; d1; ### NO d1=freg(u1,rDe)
      #d1=dd[ns];d1

      ns+1; nu; ne
      ### Plot of derivatives at uu:
      # points(uu,dd)
      # The derivative of the last spline component
      # decreases monotonously and it does not turn back.
      # This is the result of the fact that the last derivative without
      # future values reflects only tha past values and so the estimate of
      # the last derivative in the model is a bad one, it is appr. -1000.
      # Fortunatelly the derivative in our model is a well interpretable parameter.
      # It determines the direction and intensity of growth.
      # By estimating a regression line/curve from the estimated derivatives
      # we can estimate the end derivatives !!!
      # This is the benefit of interpretable parameters.

      ### c) Regression extrapolation for the last derivative
      # Leave out the last derivative
      rDe=lm(dd[1:ns]~uu[1:ns]); rDe
      #abline(rDe, col='magenta')

      de=rDe$fitted.values
      ### We at last have the corrected
      ### derivative d2 for correction
      ### of the last spline component Sk"
      ##d1=freg(u1,rDe); d1
      d1; freg(u1,rDe)
      d2=de[length(de)]; d2 #<=> de[ns] #<=> #freg(u2,rDe)
      d3=freg(u3,rDe)
      d4=freg(u4,rDe)

      if(plotAllTF){
        xy <- list(...); #print(yz)
        xlable=NULL; ylable="Derivatives";
        if(!is.null(xy$xlab)) xlable <- xy$xlab
        #if(!is.null(xy$ylab)) ylable <- xy$ylab
        #plot(xx,yy,col="gray",xlab=xlable,ylab=ylable);

        plot(0, pch='',xlim=c(uu[1],u4), ylim=c(-0.1,0.1)
             ,xlab=xlable,ylab=ylable)

        ### a) Plot derivatives of spline components:
        title('Correcting and extrapolating derivatives')
        for (i in 1:(ns)){
          frb=ifelse(i%%2==0,colors[2],colors[1]);
          lines(pL[[i]] ,xlim=c(uu[i],uu[i+1]), col=frb);
        }
        ### Plot of derivatives at uu:
        # aug 31
        points(uu,dd, col='magenta',cex=cex0)
        #abline(rDe, col='magenta')
        uuU=c(uu[1],u4)
        ddD=freg(uuU,rDe)
        lines(uuU,ddD,col='magenta')

        # The estemeted derivatives lie on the line
        #points(uu[1:ns],de[1:ns],pch=15,col='magenta')
        # The forecasted derivative also:
        points(u2,d2,pch=0,col='green',cex=cex0)
        #lines(c(uu[nu],uu[nu]),c(dd[ns+1],d2),col='green')
        arrows(uu[nu],dd[ns+1], uu[nu],d2,col='green'
               ,length=0.135,angle=12,code=2)
        points(u3,d3,pch=0,col='cyan',cex=cex0)
        points(u4,d4,pch=0,col='cyan',cex=cex0)
      }
      d14_=c(d1,d2,d3,d4)
      names(d14_)=NULL
      return ( data.frame(u14=c(u1,u2,u3,u4), d14=d14_));
    } # end forcastD1234

    forcastY1234 = function(uu,xx,yy,se,sp, dd,plotAllTF=TRUE
                            ,colors = c("red", "blue"), forecast_color="green", ...){
      d1=dd[1]; d2=dd[2]
      d3=dd[3]; d4=dd[4]
      nu=length(uu);
      ne=length(se);
      ns=length(sp); nu; ne; ns
      ##points(uu,ss$est[1:(length(ss$est)-2)],col='magenta')
      u1=uu[nu-1]; u2=uu[nu];
      h=u2-u1; u3=u2+h; u4=u3+h;
      u1;u2;u3;u4;
      ### We need for the last component
      ### 4 parameters y1,y2,d1,d2
      ### Three of them are given naturally by:
      y1=se[ne-3]
      y2=se[ne-2]

      ### 3/4 Upper and bottom Regression estimates using extrapolation
      ### a) Upper points - from j by k
      ### b) Bottom points - from j by k
      ### c) REG
      ### d) Estimates yi and di computation by extrapolation
      ###
      ### a) Upper points - from j by k
      ### a) Upper points - from j by k
      j=2; k=2
      #ye=ss$es[1:(ne-2)];ye # ~ ~
      ye=se[1:(ne-2)];ye
      yUp=c(ye[j], ye[j+1:(length(ye)-j)][1:k==k]);yUp
      uUp=c(uu[j], uu[j+1:(length(uu)-j)][1:k==k]);uUp
      #points(uUp,yUp,pch=0,col='black')
      ### b) Bottom points - from j by k
      j=1; k=2
      #ye=ss$es[1:(ne-2)];ye
      ye=se[1:(ne-2)];ye
      yBo=c(ye[j], ye[j+1:(length(ye)-j)][1:k==k]);yBo
      uBo=c(uu[j], uu[j+1:(length(uu)-j)][1:k==k]);uUp
      #points(uBo,yBo,pch=0,col='magenta')
      #  if(plotAllTF){
      #    plot(0, pch='',xlim=c(uu[1],UU[4]), ylim=c(min(yy)-0.1,max(yy)+0.1))
      #    points(xx,yy,col="gray");
      #    points(uUp,yUp,pch=0,col='magenta')
      #    points(uBo,yBo,pch=0,col='black')
      #  }
      ### c) REG
      ### Up
      rUp=lm(yUp~uUp); rUp
      #abline(rUp,col='green')
      ##yUe=freg(uUp,rUp)
      #lines(uUp,yUe,col='black')
      uu3=c(uu[1],u3)
      yy3=freg(uu3,rUp)
      #lines(uu3,yy3,pch=0,col='magenta')

      ### Bottom
      rBo=lm(yBo~uBo); rBo
      #abline(rBo,col='green')
      ##yBe=freg(uBo,rBo)
      uu4=c(uu[1],u4) ## ?? uu[1]
      yy4=freg(uu4,rBo)
      #lines(uu4,yy4,pch=0,col='green')

      ### d) Estimates yi and di by extrapolation
      y3=freg(u3,rUp); y3
      y4=freg(u4,rBo); y4
      ###print(c(y1,y2,y3,y4))
      ###
      ### 4/4 PLOT of corrected and predicted spline components
      ###
      if(plotAllTF){
        xy <- list(...); #print(yz)
        xlable=NULL; ylable=NULL;
        if(!is.null(xy$xlab)) xlable <- xy$xlab
        if(!is.null(xy$ylab)) ylable <- xy$ylab
        #plot(xx,yy,col="gray",xlab=xlable,ylab=ylable);
        plot(0, pch='',xlim=c(uu[1],u4), ylim=c(min(yy)-0.1,max(yy)+0.1)
             ,xlab=xlable, ylab=ylable)
        points(xx,yy,col="gray");
        ### a) Upper points - from j by k
        #    uu3=c(uu[1],u3)
        #    dd3=freg(uu3,rDe)
        #    points(uu3,dd3,pch=0,col='magenta')
        # aug 31
        points(uUp,yUp,col='magenta',cex=cex0)
        ### b) Bottom points - from j by k
        #???:
        points(uBo,yBo,col='black',cex=cex0)

        ### Upper
        #abline(rUp,col='magenta')
        #lines(uUp,yUe,col='magenta')
        ### Bottom
        #abline(rBo,col='black')
        #lines(uBo,yBe,col='green')
        lines(uu3,yy3,pch=0,col='magenta')
        lines(uu4,yy4,pch=0,col='black')

        ### a) Plot derivatives of spline components:
        title('Extrapolating upper and bottom function values')
        #sss=interpSpSm(xx[(1+posun):N],yy[(1+posun):N], nu-1,colors_two);
        k=length(sp)
        for (i in 1:(k)){
          frb=ifelse(i%%2==0,colors[2],colors[1]);
          lines(sp[[i]] ,xlim=c(uu[i],uu[i+1]), col=frb);
        }
        #points(u2,y2,col='green',pch=15)
        points(u3,y3,col='cyan',pch=0,cex=cex0)
        points(u4,y4,col='cyan',pch=0,cex=cex0)
      }
      y14_=c(y1,y2,y3,y4)
      names(y14_)=NULL
      return ( data.frame(u14=c(u1,u2,u3,u4), y14=y14_));
    } # end forcastY1234

    forcast2comp = function(uu,xx,yy,se,sp, dd,YY,colors = c("red", "blue"),
                            forecast_color="green", ...){
      d1=dd[1]; d2=dd[2]
      d3=dd[3]; d4=dd[4]
      nu=length(uu);
      ne=length(se);
      ns=length(sp); nu; ne; ns
      ##points(uu,ss$est[1:(length(ss$est)-2)],col='magenta')
      u1=uu[nu-1]; u2=uu[nu];
      h=u2-u1; u3=u2+h; u4=u3+h;
      u1;u2;u3;u4;
      ### We need for the last component
      ### 4 parameters y1,y2,d1,d2
      ### Three of them are given naturally by:
      y1=YY[1]; y2=YY[2]
      y3=YY[3]; y4=YY[4]
      #  if(plotAllTF){
      xy <- list(...); #print(yz)
      xlable=NULL; ylable=NULL;
      if(!is.null(xy$xlab)) xlable <- xy$xlab
      if(!is.null(xy$ylab)) ylable <- xy$ylab
      #plot(xx,yy,col="gray",xlab=xlable,ylab=ylable);
      plot(0, pch='',xlim=c(uu[1],u4), ylim=c(min(yy)-0.1,max(yy)+0.1)
           ,xlab=xlable, ylab=ylable)
      points(xx,yy,col="gray");
      ### a) Plot derivatives of spline components:
      title('Spline with forcasting')
      #sss=interpSpSm(xx[(1+posun):N],yy[(1+posun):N], nu-1,colors_two);
      k=length(sp)
      for (i in 1:(k)){
        frb=ifelse(i%%2==0,colors[2],colors[1]);
        lines(sp[[i]] ,xlim=c(uu[i],uu[i+1]), col=frb);
      }
      lines(c(uu[nu],uu[nu]),c(2.5, 2.8),col='gray')
      Plot1CompHerSp(u1,u2,25,y1,y2,d1,d2,forecast_color)
      Plot1CompHerSp(u2,u3,25,y2,y3,d2,d3,"cyan")
      Plot1CompHerSp(u3,u4,25,y3,y4,d3,d4,forecast_color)
      # szept Csaba
      return(list(uu=c(u1,u2,u3,u4), yy=yy, dd=dd))
    } # end forcast2comp

    forcast_Spline = function(xx,yy,k, plotAllTF=FALSE, colors = c("red", "blue")
                              , ...){
      posun=0;
      N <- length(xx) # szept Csaba_OK
      xy <- list(...); #print(1111111);print(xy)
      xlable=NULL; ylable=NULL;
      if(!is.null(xy$xlab)) xlable <- xy$xlab
      if(!is.null(xy$ylab)) ylable <- xy$ylab
      #plot(xx,yy,col="gray",xlab=xlable,ylab=ylable);
      #ss=interpSpSm(xx[(1+posun):N],yy[(1+posun):N],k
      #              ,colors_two,FALSE,xlab=xlable,ylab=ylable);
      ss <- cics_unif_explicit_smooth(xx[(1+posun):N],yy[(1+posun):N],k
                                      ,colors,xlab=xlable,ylab=ylable,plotTF = FALSE)
      #interpSpSm(xx[(1+posun):N],yy[(1+posun):N],k
      #             ,colors_two,FALSE,xlab=xlable,ylab=ylable);

      ###
      ### 2) Correct d2 and extrapolate d3,d4
      ###
      dOrig <- ss$est_gamma[length(ss$est_gamma)]
      uu=ss$nodes;
      se=ss$est_gamma;
      sp=ss$est_spline_polynomials;
      #?old = par()
      par(mfrow=c(2,1),mar=c(4,4,2,2))
      if(plotAllTF) par(mfrow=c(2,1),mar=c(4,4,2,2))
      uD=forcastD1234(uu,se,sp,plotAllTF,colors
                      ,xlab=xlable);
      #print(1111111);print(uD)

      #plot(uuDD$u14,uuDD$d14)

      ###
      ### 3) Extrapolate y3,y4
      ###
      DD=uD$d14; DD
      uY=forcastY1234(uu,xx,yy,se,sp, DD,plotAllTF
                      ,xlab=xlable,ylab=ylable);
      #print(222222222);print(uY)

      ###
      ### 4) Forecasting 2 components
      ###
      if(! plotAllTF) par(mfrow=c(1,1))
      YY=uY$y14; YY
      # haha
      #forcast2comp(uu,xx,yy,se,sp, DD,YY,colors=colors,forecast_color="green"
      #             ,xlab=xlable,ylab=ylable)
      #?par(old$par,old$mar)
      #return(list(uY=uY,uD=uD,dOrig=dOrig))
      return(list(UU=uY$u14,UY=uY$y14,UD=uD$d14,dOrig=dOrig))
    } # END forcast_Spline

    THP_HerSp = function(x, v1,v2,f1,f2,D1,D2){
      ### From THP_TaylorHermitePolynomials_LIBR_4g.mw
      f1*(x-v2)^2*(1-(2*(x-v1))/(v1-v2))/(v1-v2)^2+f2*(x-v1)^2*(1-(2*(x-v2))/(v2-v1))/(v2-v1)^2+D1*(x-v2)^2*(x-v1)/(v1-v2)^2+D2*(x-v1)^2*(x-v2)/(v2-v1)^2
    }

    Plot1CompHerSp = function(u1,u2,n,y1,y2,d1,d2,frb){
      uuu=seq(u1,u2,by=(u2-u1)/n)
      yyy=THP_HerSp(uuu,u1,u2,y1,y2,d1,d2);
      lines(uuu,yyy,col=frb)
    }

    freg=function(x,reg){
      reg$coefficients[1] + x*reg$coefficients[2];}

    ### End functions.

    old_par <- par(no.readonly = TRUE)
    ###
    ### yy - R dataset: AirPassengers
    ###
    yy=log10(AirPassengers); N<-length(yy);
    ###
    ### 1) xx - create corresponding x coordinates
    ###
    a <- 0; b <- 10;
    ### !!! Poucenie - derivacie sa zmensuju ak b rastie !!!
    a <- 1; b <- length(yy); # sept Csaba for cikk
    xx <- c(seq(a,b,length.out = N));

    k=2*N/12; N; k # k=24 half-years

    #uDY=forcast_Spline(xx,yy,k)
    #uDY=forcast_Spline(xx,yy,k,xlab="Years",ylab="AirPassengers")

    uDY=forcast_Spline(xx,yy,k, plotAllTF=TRUE
                       ,xlab="Years",ylab="AirPassengers")
    par(old_par)
    return(uDY)
  } # END forecast_demo

  fd <- forecast_demo0()

  { ### Call cics_unif_explicit_smooth
    yy <- as.vector( log10(AirPassengers) ); yy
    xx <- c(1:length(yy))

    k <- 24
    clrs <- c("blue","red")
    aPsp <- cics_unif_explicit_smooth(xx,yy,k,clrs
                                      ,title=paste0("Smoothing log(AirPassengers) with a ",
                                                    k,"-component uniform CICS"),plotTF = FALSE)
    aPsp
  }
  { ### ggplot of smoothing
    pp <- aPsp$pp;
    pp  # - cics_unif_explicit_smooth
  }
  if(F){
    df <- data.frame(xx=xx, yy=yy)
    pp_ <- pp + geom_line(aes(xx, yy), data = df, color='gray'); pp_
    ##pp <- pp_ + pp; pp
  }

  { ### Add 3 segments with forecast:

    uF <- fd$UU; uF
    yF <- fd$UY; yF
    dF <- fd$UD; dF
    i <- 1; pp <- Plot1CompHerSp(pp, uF[i],uF[i+1],25,yF[i],yF[i+1],dF[i],dF[i+1],'magenta'); pp
    i <- 2; pp <- Plot1CompHerSp(pp, uF[i],uF[i+1],25,yF[i],yF[i+1],dF[i],dF[i+1],'black');   pp
    i <- 3; pp <- Plot1CompHerSp(pp, uF[i],uF[i+1],25,yF[i],yF[i+1],dF[i],dF[i+1],'magenta'); pp

    ### Vertical
    xv <- xx[length(xx)]
    pp <- pp + geom_line(aes(uv,yv), data.frame(uv=c(xv,xv),yv=c(2.5, 2.8)), color='gray')

    tit <- paste0(pp$labels$title,
                  ",\ncorrected last segment and two forecasting segments")
    pp <- pp + ggtitle(tit) +
      theme(plot.title = element_text(size=10)) +
      xlab("XXX") + ylab("YYY")
  }
  print(pp)
} ### forecast_demo
