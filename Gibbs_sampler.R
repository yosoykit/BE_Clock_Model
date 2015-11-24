gibbs <- function(x,L=100,K=100,skip=10,plot=T) {
  if (!is.list(x)) {
	cat("x is not a list! see help file", "\n")
	return()
      }
  names(x)[1] <- "label"
  names(x)[2] <- "est"
  names(x)[3] <- "low"
  names(x)[4] <- "upp"

  npar <- length(x$est)

  if (npar <= 0) stop('no. of parameters must be >= 1')

  # *** initialize graphical output
  if(plot==TRUE) {
      # par(mfrow=c(npar+1,1),bg="grey")
      par(mfrow=c(8,1),bg="grey")
  }

  ## starting point
  y = x$est
  ## data-specific setup
  ## provide: ids, ID, M, a 
  xj = beM[ids,ID]; s = y[M+1]
  
  x.mon <- matrix(0,ncol=npar,nrow=L)

  for (l in 1:L) {
    
    ## full conditionals for b_j's: x[1:M]
    for(m in 1:M) {
      xm = xj[m]
      SD = y[M+2]
      mb = mean_b_i[m] # mean b prior
      sd = sd_b_i[m]   # sd b prior
      r = (xm - alpha_SQ[m] - b_SQ[m]*s)/(a-s)
      denom = (SD^2+((a-s)*sd)^2)
      b.mean = (SD*SD*mb+sd*sd*r*(a-s)^2)/denom
      b.sigm = sd*SD/sqrt(denom)
      y[m] = rnorm(1,mean = b.mean, sd = b.sigm)
    }

    ## full conditional for s: x[M+1]
    test=a+1
    while( test >a || test < 0){
      D = sum((y[1:M]-b_SQ)*(y[1:M]-b_SQ))
      g = alpha_SQ + y[1:M]*a - xj
      E = sum(g*(y[1:M]-b_SQ))
      test = rnorm(1, mean = E/D, sd = y[M+2]/sqrt(D))
    }
    y[M+1] =test 
    s = y[M+1]

    ## full conditional for sigma_BE: x[M+2]
    gam_1=0.6075155
    gam_2=0.1173199
    gam.shape = 0.5*M
    gam.rate = 0.5*sum((xj-alpha_SQ-b_SQ*s-y[1:M]*(a-s))^2)
    dum = rgamma(1, shape=gam_1+gam.shape, rate=gam_2+gam.rate)
    y[M+2] = 1/sqrt(dum)

    x.mon[l,] <- y
    # PLOTTING OF RUNS ##
    if(l%%(100*skip) == 0) {
      if(plot==TRUE) {
        if(l < K) {brncol <- 3} else {brncol <- 2}
        n.skip.1 <- seq(skip,l,skip)
        n.skip.2 <- seq(skip,min(l,K),skip)
        
        par(mar=c(0, 5, 0.6, 4) + 0.1)
        plot(x.mon[n.skip.1,1], type='l', xlab = " ", ylab = x$label[1],col=2, main = ID)
        lines(x.mon[n.skip.2,1],col=brncol) # burn-in cycles

        # for (i in 1:(npar-1)) {
        for (i in 2:6) {
          par(mar=c(0, 5, 0, 4) + 0.1)
          plot(x.mon[n.skip.1,i], type='l', xlab = " ", ylab = x$label[i], col=2)
          lines(x.mon[n.skip.2,i],col=brncol) #pilot cycles
        }
        par(mar=c(0, 5, 0, 4) + 0.1)
        plot(x.mon[n.skip.1,npar-1], type='l', xlab = " ", xaxt='n', ylab = x$label[npar-1], col=2)
        lines(x.mon[n.skip.2,npar-1],col=brncol) #pilot cycles
        
        par(mar=c(0.1, 5, 0, 4) + 0.1)
        plot(x.mon[n.skip.1,npar], type='l', xlab = " ", xaxt='n', ylab = x$label[npar], col=2)
        lines(x.mon[n.skip.2,npar],col=brncol) #pilot cycles
      }
    }
  }
  return(x.mon)
}