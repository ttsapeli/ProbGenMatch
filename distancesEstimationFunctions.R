# input
# N: the number of samples
# vars: the number of deterministic counding variables
# pvars: the number of probabilistic counding variables
# X: deterministic confounding variables of size N x vars
# pX: probabilistic confounding variables of size pvars*N*(N-1)/2 x 10. If Tr is also probabilistic, pX corresponds to the confounding variables devided by Treatment CDF
# Tr: treatment. If it is probabilistic the size is N*(N-1)/2 x 10 otherwise N x 1
#output:
#an array with the distances between each pair of units of size N*(N-1)/2 x (vars+pvars)
getDistances <- function(X, pX, Tr, N, vars, pvars)
{
  print("Get distances with area method!")
  dists = c()
  #if there are some deterministic confounding variables
  if(vars>0)
  {
    if (dim(Tr)[2]==1)
      dists = getDistances.Deterministic(X, Tr, vars)
    else
    {
      xAxis <- buildXAxis(Tr)
      dists = getDistances.Probabilistic.Tr(X, Tr, vars, xAxis)
    }
      
  }
  #if there are some probabilistic confounding variables
  if (pvars>0)
  {
    if (dim(Tr)[2]==1)
    {
      xAxis <- buildXAxis(pX)
      dists = cbind(dists, getDistances.Probabilistic.C(pX, Tr, pvars, xAxis))
    }
    else
    {
      xAxis <- buildXAxis(pX)
      dists = cbind(dists, getDistances.Probabilistic.All(pX, pvars, xAxis))
    }
      
  }
  return(dists)
}


getDistances.Deterministic <- function(X, Tr, vars)
{
  N = dim(Tr)[1]
  dists = cbind()
  for (i in 1:vars)
  {
    d = c()
    for (j in 1:(N-1))
    {
      d = c(d, ((X[j,i] - X[(j+1):N,i])/(Tr[j] - Tr[(j+1):N]))^2)
    }
    dists = cbind(dists, d)
  }
  return(dists)
}

getDistances.Probabilistic.Tr <- function(X, Tr, vars, xAxis)
{
  N = dim(X)[1]
  dists = cbind()
  trDiffAll = c()
  for (j in 1:(N-1))
  {
    # get the the pairs of j with each other of the N-j units
    startPos = (j-1)*N - (j-1)*j/2
    # for each pair (j, (e-startPos+j) ) estimate the distance
    for (e in (startPos+1):(startPos+N-j))
    {
      trDiffAll = c(trDiffAll, cdfDistance(Tr[e,], xAxis))
    }
  }
  trDiffAll = max(trDiffAll) - trDiffAll
  trDiffAll = (trDiffAll - mean(trDiffAll))/sqrt(var(trDiffAll))
  trDiffAll = trDiffAll - min(trDiffAll)
  
  # for each deterministic confounding variable
  for (i in 1:vars)
  {
    d = c()
    for (j in 1:(N-1))
    {
      startPos = (j-1)*N - (j-1)*j/2
      d = c(d, ((X[j,i] - X[(j+1):N,i])*trDiffAll[(startPos+1):(startPos+N-j)])^2)
    }
    dists = cbind(dists, d)
  }
  return(dists)
}

getDistances.Probabilistic.C <- function(pX, Tr, pvars, xAxis)
{
  N = dim(Tr)[1]
  dists = cbind()
  # for each stochastic confounding variable
  for (i in 1:pvars)
  {
    d = c()
    confDistAll = c()
    # get the position of the stochastic variable in the vector
    Xstart = (i-1)*N*(N-1)/2 + 1
    Xend = i*N*(N-1)/2
    # store the CDFs related to the pairs of units of this stochastic variable at variable X
    X = pX[Xstart:Xend,]
    for (j in 1:(N-1))
    {
      startPos = (j-1)*N - (j-1)*j/2
      # for each pair (j, (e-startPos+j) ) estimate the distance
      for (e in (startPos+1):(startPos+N-j))
      {
        confDistAll = c(confDistAll, cdfDistance(pX[e,], xAxis))
      }
    }
    confDistAll = max(confDistAll) - confDistAll
    confDistAll = (confDistAll - mean(confDistAll))/sqrt(var(confDistAll))
    confDistAll = confDistAll - min(confDistAll)
    for (j in 1:(N-1))
    {
      startPos = (j-1)*N - (j-1)*j/2
      d <- c(d, (confDistAll[(startPos+1):(startPos+N-j)]/(Tr[j] - Tr[(j+1):N]))^2)
    }
    dists = cbind(dists, d)
  }
  return(dists)
}

getDistances.Probabilistic.All <- function(pX, pvars, xAxis)
{
  N = dim(Tr)[1]
  dists = cbind()
  # for each stochastic confounding variable
  for (i in 1:pvars)
  {
    d = rbind()
    # get the position of the stochastic variable in the vector
    Xstart = (i-1)*N*(N-1)/2 + 1
    Xend = i*N*(N-1)/2
    # store the CDFs related to the pairs of units of this stochastic variable at variable X
    X = pX[Xstart:Xend,]
    for (j in 1:N)
    {
      startPos = (j-1)*N - (j-1)*j/2
      # for each pair (j, e) estimate the distance
      for (e in (startPos+1):(startPos+N-j))
      {
        d <- c(d, (cdfDistance(pX[e,], xAxis))^2)
      }
    }
    dists = cbind(dists, d)
  }
  return(dists)
}


buildXAxis <- function(pVar)
{
  x <- as.vector(pVar)
  x <- round(x, digits=4)
  x <- seq(min(x), max(x), (max(x) - min(x))/100)
  return(x)
}


getCDFReduced <- function(myCdf)
{
  y <- seq(0.01, 0.99, 0.01)
  reducedY <- c()
  reducedCdf <- unique(myCdf)
  for ( u in reducedCdf)
  {
    reducedY = c(reducedY, y[tail(which(myCdf==u), 1)])
  }
  return(list(myCdf = reducedCdf, y = reducedY))
}
cdfDistance <- function(myCdf, x)
{
  l = getCDFReduced(myCdf)
  myCdf = l$myCdf
  y = l$y
  if(round(myCdf[1], 2)>=round(x[length(x)],1))
  {
    return(1)
  }
  
  h = which.min(abs(myCdf[length(myCdf)] - x))
  if(h<length(x))
  {
    h = h+1
    myCdf = c(myCdf, x[h:length(x)])
    y = c(y, rep(1, length(x[h:length(x)])))
  }
  
  interp <- approx(myCdf, y, x)
  yInterp <- interp$y
  yInterp[is.na(yInterp)] <- 0
  d <- sum(yInterp)
  # d <- (length(x) - d)/length(x)
  return(d)
}

