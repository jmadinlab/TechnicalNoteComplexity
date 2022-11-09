###########################################################################
## R code for Loke and Chisholm 2022 
## SEE README FILE FOR CITATION INFORMATION
###########################################################################

gradient <- function(x)
  {
  n = length(x)
  c(x[2]-x[1],x[-(1:2)]-x[-((n-1):n)],x[n]-x[n-1])
  }

# helper function for computing fractal dimension using box counting
est_D_box_counting_method_helper <- function(X,n_box_factor=1)
  {
  # convert to 2D
  X2 = (X>median(X))

  L = dim(X2)[1]

  # box size
  epsilon = L-1

  epsilons = numeric()
  n_box_cover_all = numeric()

  # do box count
  while ( epsilon >= 1 )
    {
    n_box_cover = 0

    box_x0 = seq(1,L,epsilon) - runif(1,0,epsilon)
    box_x1 = box_x0 + epsilon
    n_box_x = length(box_x0)

    box_x0 = pmax(box_x0,1)
    box_x1 = pmin(box_x1,L)

    box_y0 = seq(1,L,epsilon) - runif(1,0,epsilon)
    box_y1 = box_y0 + epsilon
    n_box_y = length(box_y0)

    box_y0 = pmax(box_y0,1)
    box_y1 = pmin(box_y1,L)

    for ( ix in 1:n_box_x )
      {
      for ( iy in 1:n_box_y )
        {
        temp = X2[box_x0[ix]:box_x1[ix],box_y0[iy]:box_y1[iy]]
        sum_temp = sum(temp)
        d_box = (0<sum_temp & sum_temp<prod(dim(temp)))
        n_box_cover = n_box_cover + d_box
        }
      }

    epsilons = c(epsilons,epsilon)
    n_box_cover_all = c(n_box_cover_all,n_box_cover)
    #cat('epsilon = ', epsilon,'; N(epsilon) = ',n_box_cover,'\n',sep='')
    epsilon = epsilon/2
    }

  return(list(epsilons=epsilons,ns=n_box_cover_all))
  }

# function to compute fractal dimension using box counting
est_D_box_counting_method <- function(X,n_box_factor=1)
  {
  temp = est_D_box_counting_method_helper(X,n_box_factor)
  b = lm(log(temp$ns)~log(temp$epsilons))$coef[2]
  D_est = 1-b
  return(D_est)
  }

# function to compute fractal dimension using box counting at intermediate scales
est_D_box_counting_method_int_scales <- function(X,n_box_factor=1)
  {
  temp = est_D_box_counting_method_helper(X,n_box_factor)
  Ds_est = -gradient(log(temp$ns))/gradient(log(temp$epsilons))
  D_est = Ds_est[ceiling(length(Ds_est)/2)]+1
  return(D_est)
  }

# function to calculate the range of heights at a given resolution
calc_dH <- function(fM,L1)
{
  L = dim(fM)[1]-1
  stopifnot(dim(fM)[2]-1 == L)
  
  k = L/L1
  
  dHs = matrix(NA,k,k)
  
  for ( i in 1:k )
    for ( j in 1:k )
    {
      zs = fM[(i-1)*L1+(1:(L1+1)),(j-1)*L1+(1:(L1+1))]
      dHs[i,j] = diff(range(zs))
    }   
  
  return(mean(dHs))
  }

# helper function for computing fractal dimension using the method of variation
est_D_variation_method_helper <- function(fM,L0,L1)
  {
  L = dim(fM)[1]-1
  stopifnot(dim(fM)[2]-1 == L)
  log2L = round(log2(L))
  stopifnot(2^log2L == L)

  # specify resolutions to use
  Li_s = 2^(log2(L0):log2(L1))
  mean_dHs = numeric(length(Li_s))

  n = length(Li_s)

  for ( i in 1:n )
    {
    dHs = calc_dH(fM,Li_s[i])
    # store mean patch height range at this resolution
    mean_dHs[i] = dHs
    }
  
  return(list(Li_s=Li_s,mean_dHs=mean_dHs))
  }

# function to compute fractal dimension using box counting
est_D_variation_method <- function(X,L0,L1)
  {
  temp = est_D_variation_method_helper(X,L0,L1)
  S = lm(log(temp$mean_dHs)~log(temp$Li_s))$coef[2]
  D_est = 3-S
  return(D_est)
  }

# function to compute fractal dimension using box counting at intermediate scales
est_D_variation_method_int_scales <- function(X,L0,L1)
  {
  temp = est_D_variation_method_helper(X,L0,L1)
  bs = gradient(log(temp$mean_dHs))/gradient(log(temp$Li_s))
  D_est = 3-bs[ceiling(length(bs)/2)]
  return(D_est)
  }

# function to compute fractal dimension using method of Torres-Pulliza et al.;
# only the resolutions L0 to L1 are used
est_D_Torres_Pulliza_method <- function(fM,L0,L1)
  {
  L = dim(fM)[1]-1
  stopifnot(dim(fM)[2]-1 == L)
  log2L = round(log2(L))
  stopifnot(2^log2L == L)

  # specify resolutions to use
  Li_s = 2^c(log2(L0):log2(L1))
  mean_dHs = numeric(length(Li_s))

  n = length(Li_s)

  for ( i in 1:n )
    {
    dHs = calc_dH(fM,Li_s[i])
    # store mean patch height range at this resolution
    mean_dHs[i] = dHs
    }
  
  # calculate fractal dimension from the log-log slope of height range versus resolution
  S = mean(diff(log(mean_dHs))/diff(log(Li_s)))

  D_est = 3-S
  return(list(D_est = D_est, mean_H_min = mean_dHs[1], mean_H_max = mean_dHs[9]))
}

R_func <- function(H0, L0) {
  sqrt((H0^2) / (2 * L0^2) + 1)
}

