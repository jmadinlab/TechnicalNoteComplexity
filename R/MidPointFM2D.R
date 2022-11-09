###########################################################################
## R code for Loke and Chisholm 2022 
## SEE README FILE FOR CITATION INFORMATION
###########################################################################

f3 <- function(delta,x0,x1,x2)
  {
  if ( is.matrix(x0) ) return((x0+x1+x2)/3+t(matrix(rnorm(prod(dim(x0)),0,delta),dim(x0)[2],dim(x0)[1])))
  return((x0+x1+x2)/3+rnorm(length(x0),0,delta))
  }

f4 <- function(delta,x0,x1,x2,x3)
  {
  if ( is.matrix(x0) ) return((x0+x1+x2+x3)/4+t(matrix(rnorm(prod(dim(x0)),0,delta),dim(x0)[2],dim(x0)[1])))
  return((x0+x1+x2+x3)/4+rnorm(length(x0),0,delta))
  }

# from Saupe (1988)
MidPointFM2D <- function(maxlevel,H,sigma)
  {
  #set.seed(200)
  N = 2^maxlevel
  addition = T
  
  delta = sigma
  X = matrix(0,N+1,N+1)
  X[1,1] = rnorm(1,0,delta)
  X[1,N+1] = rnorm(1,0,delta)
  X[N+1,1] = rnorm(1,0,delta)
  X[N+1,N+1] = rnorm(1,0,delta)
  D = N
  d = N/2
  
  for ( stage in 1:maxlevel )
    {
    delta = delta*(0.5^(H/2))
    xs = seq(d,N-d,D)
    ys = seq(d,N-d,D)
    X[xs+1,ys+1] = f4(delta,X[xs+d+1,ys+d+1],X[xs+d+1,ys-d+1],X[xs-d+1,ys+d+1],X[xs-d+1,ys-d+1])

    delta = delta*(0.5^(H/2))
    
    xs = seq(d,N-d,D)
    X0 = X
    X[xs+1,1] = f3(delta,X0[xs+d+1,1],X0[xs-d+1,1],X0[xs+1,d+1])
    X[xs+1,N+1] = f3(delta,X0[xs+d+1,N+1],X0[xs-d+1,N+1],X0[xs+1,N-d+1])
    X[1,xs+1] = f3(delta,X0[1,xs+d+1],X0[1,xs-d+1],X0[d+1,xs+1])
    X[N+1,xs+1] = f3(delta,X0[N+1,xs+d+1],X0[N+1,xs-d+1],X0[N-d+1,xs+1])

    if ( D <= N-d )
      {
      X0 = X

      xs = seq(d,N-d,D)
      ys = seq(D,N-d,D)
      X[xs+1,ys+1] = f4(delta,X0[xs+1,ys+d+1],X0[xs+1,ys-d+1],X0[xs+d+1,ys+1],X0[xs-d+1,ys+1])

      xs = seq(D,N-d,D)
      ys = seq(d,N-d,D)
      X[xs+1,ys+1] = f4(delta,X0[xs+1,ys+d+1],X0[xs+1,ys-d+1],X0[xs+d+1,ys+1],X0[xs-d+1,ys+1])
      }

    D = D/2
    d = d/2
    }
  return(X)
  }

