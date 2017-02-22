
#This method is called when all variables are deterministic
#Tr: treatment values
#rr: matched treatment - control pairs
#X: confounding variables
#returns just the difference on the confounding variables for all the pairs
getDecileStats.Deterministic <- function(Tr, rr, X, N, nvars)
{
	pairsNum = dim(rr)[1]
  	diffDecilesPerVar = rbind()
  	treatmentPrecentilesPerVar = rbind()
  	controlPrecentilesPerVar = rbind()
  	#for each deterministic confounding variable
  	for (j in seq(nvars))
  	{	  
    
    	indexTr = rr[,1]
    	indexCt = rr[,2]
    	treatmentPrecentiles = X[indexTr,j]/(Tr[indexTr]-Tr[indexCt])
    	controlDeciles = X[indexCt,j]/(Tr[indexTr]-Tr[indexCt])
    
    	indx = round(seq(N/10, N, N/10))
    	sortedDiffs = sort(abs(treatmentPrecentiles-controlDeciles))
    	diffDecilesPerVar = rbind(diffDecilesPerVar, sortedDiffs)
  	}
  
  	return(diffDecilesPerVar)
}

# Given a vector of values and their corresponding probabilities, returns the deciles
calculateDeciles <-function(values, probs, N)
{
	sValues = sort.int(values, index.return = TRUE)
  	values = sValues$x
  	iValues = sValues$ix
  	probs = probs[iValues]
  	deciles = c()
  	probsCumSum = cumsum(probs)/N
  	for (p in seq(0.1, 1, 0.1))
	{
    	deciles = c(deciles, values[tail(which(probsCumSum<=p), 1)])
	}
  	return(deciles)
}
#This method is called when all confounding variables are deterministic and treatment is probabilistic
#TrPDF: PDF of 1/(Z_tr - Z_ct)
#rr: matched treatment - control pairs
#X: deterministic confounding variables
getDecileStats.Probabilistic.Tr <- function(TrPDF, rr, X, N, nvars)
{
	pairsNum = dim(rr)[1]
  	treatmentPrecentilesPerVar = rbind()
  	controlPrecentilesPerVar = rbind()
  	#for each confounding variable
  	for (j in seq(nvars))
  	{  
    	diffDeciles.Values =  c()
    	diffDeciles.Probs = c()
    	#for each matched pair
    	for (i in seq(pairsNum))
		{
			indexTr = rr[i,1]
      		indexCt = rr[i,2]
			#get the position (row) of this pair at the TrCDF matrix
			pairPos = (indexTr-1)*N - (indexTr)*(indexTr-1)/2 + indexCt-indexTr;
      		#deciles of X_tr/(Z_tr - Z_ct)
      		diffDeciles.Values = c(diffDeciles.Values, abs(X[indexTr,j]-X[indexCt,j])*TrPDF$values)
      		diffDeciles.Probs = c(diffDeciles.Probs, TrPDF$probs[pairPos,])
			#deciles of X_ct/(Z_tr - Z_ct)
    	}
		treatmentPrecentilesPerVar = rbind(treatmentPrecentilesPerVar, calculateDeciles(diffDeciles.Values, diffDeciles.Probs, pairsNum))
	}
	return(treatmentPrecentilesPerVar)
}

#This method is called when some (or all) confounding variables are probabilistic and treatment is deterministic
#Tr: treatment values
#rr: matched treatment - control pairs
#X: deterministic confounding variables
#PX: PDF of (X_tr - X_c)/(Z_tr - Z_ct) for each probabilistic confounding variable X
getDecileStats.Probabilistic.C <- function(Tr, rr, X, XPDF, N, nvars, npvars)
{
	pairsNum = dim(rr)[1]
  	#for each deterministic confounding variable
  	diffDecilesPerVar <- getDecileStats.Deterministic(Tr, rr, X, N, nvars)
  	biasDeciles = rbind()
  	#for each probabilistic confounding variable
  	for (j in seq(npvars))
	{
    	diffDeciles.Values =  c()
    	diffDeciles.Probs = c()
    	#for each matched pair
    	for (i in seq(pairsNum))
		{
      	  	indexTr = rr[i,1]
      		indexCt = rr[i,2]
      	  	#get the position (row) of this pair at the XPDF matrix
      	  	pairPos = (indexTr-1)*N - (indexTr)*(indexTr-1)/2 + indexCt-indexTr;
      	  	diffDeciles.Values = c(diffDeciles.Values, XPDF$values/(Tr[indexTr,j]-Tr[indexCt,j]))
      	  	diffDeciles.Probs = c(diffDeciles.Probs, XPDF$probs[pairPos,])
      	  	#deciles of X_ct/(Z_tr - Z_ct)
    	}
		diffDecilesPerVar = rbind(diffDecilesPerVar, calculateDeciles(diffDeciles.Values, diffDeciles.Probs, pairsNum))
	}
	return(diffDecilesPerVar)
}

#This method is called when the treatment variable and some (or all) confounding variables are probabilistic
#Tr: CDF of (X_tr - X_c)/(Z_tr - Z_ct)
#rr: matched treatment - control pairs
#X: deterministic confounding variables
#PX: CDF of (X_tr - X_c)/(Z_tr - Z_ct) for each probabilistic confounding variable X
getDecileStats.Probabilistic <- function(Tr, rr, X, PX, N, nvars, npvars)
{
	pairsNum = dim(rr)[1]
  	diffDecilesPerVar = rbind()
 	#for each deterministic confounding variable
  	for (j in seq(nvars))
  	{  
    	#each line of matrixes treatmentPrecentiles and controlDeciles correspond the potential deciles values i.e. 
    	#for each pair of samples (i, j), there is a line in treatmentPrecentiles describing the deciles of (X(i)/(T(i)-T(j))): c1(i), c2(i), ..., c99(i). 
    	#and a line in controlDeciles describing the deciles of (X(j)/(T(i)-T(j))): c1(j), c2(j), ..., c99(j). 
    	treatmentPrecentiles = rbind()
    	controlDeciles = rbind()
    	#for each matched pair
    	for (i in seq(pairsNum))
		{
      	  	indexTr = rr[i,1]
      		indexCt = rr[i,2]
			#get the position (row) of this pair at the Tr matrix
      		pairPos = (indexTr-1)*N - (indexTr)*(indexTr-1)/2 + indexCt-indexTr;
      	  	#deciles of X_tr/(Z_tr - Z_ct)
      		biasDeciles.treatmentPerc = X[indexTr,j]*Tr[pairPos,]
      	  	#deciles of X_ct/(Z_tr - Z_ct)
      	  	biasDeciles.controlPerc = X[indexCt,j]*Tr[pairPos,]
      	  	treatmentPrecentiles = rbind(treatmentPrecentiles, biasDeciles.treatmentPerc)
      	  	controlDeciles = rbind(controlDeciles, biasDeciles.controlPerc)
		}
    	#diffDecilesPerVar contains one line for each confounding variable. Each line will be: mean(c1), mean(c2), .... mean(c99)
    	diffDecilesPerVar = rbind(diffDecilesPerVar, abs(colMeans(treatmentPrecentiles)-colMeans(controlDeciles)))
	}
  
  	#for each probabilistic confounding variable
  	for (j in seq(npvars))
  	{  
    	treatmentPrecentiles = rbind()
    	controlDeciles = rbind()
    	#for each matched pair
    	for (i in seq(pairsNum))
		{
      	  	indexTr = rr[i,1]
      		indexCt = rr[i,2]
      	  	#get the position (row) of this pair at the TrCDF matrix
      	  	pairPos = (indexTr-1)*N - (indexTr)*(indexTr-1)/2 + indexCt-indexTr;
      	  	pairPos = pairPos + (j-1)*(N^2 - N*(N+1)/2)
      	  	#deciles of (X_tr-X_ct)/(Z_tr - Z_ct)
      	  	biasDeciles = rbind(biasDeciles, PX[pairPos,])
		}
    	diffDecilesPerVar = rbind(diffDecilesPerVar, colMeans(biasDeciles))
  	}
  	return(diffDecilesPerVar)
}


probabilisticQQStats <- function(decilesDiffPerVar)
{
  	mediansPerVar = rowMedians(decilesDiffPerVar)
  	meanPerVar = rowMeans(decilesDiffPerVar)
  	maxPerVar = rowMaxs(decilesDiffPerVar)
  	return(list(medianDiff=mediansPerVar, meanDiff=meanPerVar, maxDiff=maxPerVar))
}
