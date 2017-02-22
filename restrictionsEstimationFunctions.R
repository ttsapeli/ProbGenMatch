

#Args: 
# 	CDFPerPair: the CDF for each pair of samples estimated by function getCDFPerPair
# 	probabilityThreshold: a threshold to the probability that the treatment distance of two units is less than distanceThreshold multiplied by N (N is the number of quantiles)
# 	N: the number of quantiles per pair of units
# 	epsilon: the value that is added to the treatment difference in order to avoid NaNs
# Returns:
# 	restrictions: a binary vector of size N*(N-1)/2 x 1. There is a binary value for each units pair. For N units, there N*(N-1)/2 pairs in total. If the i-th element of vector restrictions is '1', the i-th pair of units is not allowed to be matched
getRestrictions.Tr.Prob <- function(CDFPerPair, distanceThreshold, probabilityThreshold, N, epsilon)
{
	probSame = matrix(0, dim(CDFPerPair))
	#Ps[i] the probability that the quantity 1/dTr[i] is less than the threshold, where dTr[i] the treatment difference of pair i,
	# probSame[i] = Ps[i]*N
	for (i in seq(nrow(CDFPerPair)))
    	probSame[i] = which(CDFPerPair[i,]>=((1/(distanceThreshold+epsilon))))[1]
	#set probSame equal to N for the pairs that have 0 probability to have the same treatment
	probSame[is.na(probSame)] <- N
  	#subtract probSame from N i.e. probSame[i] = N - Ps[i] * N
	probSame = N - probSame
	restrictions = matrix(0, nrow(CDFPerPair))
  	#set the restrictions for the pairs with probability to have the same treatment larger than probabilityThreshold equal to 1
  	restrictions[which(probSame>probabilityThreshold)] = 1
    print("Number of restrictions (probabilistic treatment variable):")
    print(length(which(restrictions==1)))
	return(restrictions)
}


#Args: 
# 	X: either a vector with the tretment values or an array with the confounding variables values for all units
# 	distanceThreshold: a threshold defining the maximum/minimum allowed distance on values of the confounding/treatment variables of matched units 
# 	restrictions: a binary vector of size N*(N-1)/2 x 1 with the restrictions that have been already defined from other methods
# 	treatmentFlag: if 1 X corresponds to the treatment variable otherwise to the confounding variables
# Returns:
# 	restrictions: a binary vector of size N*(N-1)/2 x 1. There is a binary value for each units pair. For N units, there N*(N-1)/2 pairs in total. If the i-th element of vector restrictions is '1', the i-th pair of units is not allowed to be matched
getRestrictions.Deterministic <- function(X, distanceThreshold, restrictions, treatmentFlag)
{
	N = length(X)
  	intdxSt = 0
  	indxEnd = 0
  	
  	for (i in seq(N-1))
  	{
    	indxSt = indxEnd+1
    	indxEnd = indxSt+N-i-1
    	xDiff = abs(X[i+1:N] - X[i])
    	xDiff = xDiff[1:(N-i)]
    	if(treatmentFlag)
      	  restrictionIndx = which(xDiff<distanceThreshold)
    	else
      	  restrictionIndx = which(xDiff>distanceThreshold)
    	thisSampleRestrictions = restrictions[indxSt:indxEnd]
    	thisSampleRestrictions[restrictionIndx] = 1
    	restrictions[indxSt:indxEnd] = thisSampleRestrictions
	}
	print("Number of restrictions (deterministic variables):")
	print(length(which(restrictions==1)))
	return(restrictions)
}
#Args: 
#	CDFPerPair: the CDF for each pair of samples estimated by function getCDFPerPair
# 	probabilityThreshold: a threshold to the probability that the treatment distance of two units is less than distanceThreshold multiplied by N (N is the number of quantiles)
# 	N: the number of quantiles per pair of units
# 	epsilon: the value that is added to the treatment difference in order to avoid NaNs
# Returns:
# 	restrictions: a binary vector of size N*(N-1)/2 x 1. There is a binary value for each units pair. For N units, there N*(N-1)/2 pairs in total. If the i-th element of vector restrictions is '1', the i-th pair of units is not allowed to be matched
getRestrictions.C.Prob <- function(CDFPerPair, distanceThreshold, probabilityThreshold, N)
{
	probSmallDistance = matrix(0, dim(CDFPerPair))
	#Ps[i] the probability that the ith pair of units has distance smaller than distanceThreshold on the 
  	#examined confounding variable, probSmallDistance[i] = N - Ps[i]*N
  	for (i in seq(nrow(CDFPerPair)))
    	probSmallDistance[i] = N - which(CDFPerPair[i,]>=distanceThreshold)[1]
  	#set probSmallDistance equal to 0 for the pairs that have 0 probability to have distance more than distanceThreshold
  	probSmallDistance[is.na(probSmallDistance)] <- 0
  	restrictions = matrix(0, nrow(CDFPerPair))
  	#set the restrictions for the pairs with probability to have the same treatment larger than probabilityThreshold equal to 1
  	restrictions[which(probSmallDistance>probabilityThreshold)] = 1
  	print("Number of restrictions (probabilistic confounding variables):")
  	print(length(which(restrictions==1)))
  	return(restrictions)
}



# Estimates pairs of units that are not allowed to be matched based on the thresholds provided by the user
# Args:
# 	Tr: treatment variable (either probabilistic or deterministic)
# 	PX: probabilistic confounding variables
# 	X: deterministic confounding variables
# 	distanceThreshold.Tr: threshold defining the minimum allowed difference on treatment value
# 	distanceProb.Tr: threshold defining the maximum allowed probability that the treatment difference is smaller than distanceThreshold.Tr
# 	distanceThreshold.C: threshold defining the minimum allowed difference on treatment value
# 	distanceProb.C: threshold defining the maximum allowed probability that the treatment difference is smaller than distanceThreshold.Tr
# 	epsilon: 
# Returns:
# 	restrictions: a binary vector of size N*(N-1)/2 x 1. There is a binary value for each units pair. For N units, there N*(N-1)/2 pairs in total. If the i-th element of vector restrictions is '1', the i-th pair of units is not allowed to be matched
getRestrictions <-function(Tr, PX, X, distanceThreshold.Tr, distanceProb.Tr, distanceThreshold.C, distanceProb.C, epsilon)
{
	restrictions = NULL
  	# if both Tr and at least one confounding variable are probabilistic, then Tr is the CDF for both treatment and confounding
  	# then, restrictions are estimated for all the probabilistic confounders and the treatment
  	if ((dim(Tr)[2]>1)&&(!is.null(PX)))
  	{
    	distanceThreshold = NULL
    	if ((is.null(distanceThreshold.Tr))&&(!is.null(distanceThreshold.C)))
    	{
			distanceThreshold = distanceThreshold.C
      	  	distanceProb = distanceProb.C
		}
    	else if((!is.null(distanceThreshold.Tr))&&(is.null(distanceThreshold.C)))
		{
			distanceThreshold = distanceThreshold.Tr
      	  	distanceProb = distanceProb.Tr
		}
    	else if((!is.null(distanceThreshold.Tr))&&(!is.null(distanceThreshold.C)))
		{
			distanceThreshold = distanceThreshold.Tr/distanceThreshold.C
      	  	distanceProb = min(distanceProb.C, distanceProb.Tr)
		}
    	if (!is.null(distanceThreshold))
			restrictions = getRestrictions.Tr.Prob(Tr, distanceThreshold, distanceProb*dim(Tr)[2], dim(Tr)[2], epsilon)
	}
  	# if the treatment is deterministic but there are probabilstic confounders, 
  	# restrictions are estimated for the probabilistic confounding variables and then, for the deterministici treatment
  	else if((!is.null(PX))&&(!is.null(distanceThreshold.C)))
  	{
		restrictions = getRestrictions.C.Prob(PX, distanceThreshold.C, distanceProb.C*dim(PX)[2], dim(PX)[2])
    	if (!is.null(distanceThreshold.Tr))
			restrictions = getRestrictions.Deterministic(Tr, distanceThreshold.Tr, restrictions, 1)
	}
	# if the treatment is probabilistic and all confounders deterministic, restrictions are estimated for the treatment variable
  	else if((dim(Tr)[2]>1)&&(!is.null(distanceThreshold.Tr)))
	{
		restrictions = getRestrictions.Tr.Prob(Tr, distanceThreshold.Tr, distanceProb.Tr*dim(Tr)[2], dim(Tr)[2], epsilon)
	}
	# if everything is determistic, restrictions are estimated for the treatment variable
	else if (!is.null(distanceThreshold.Tr))
	{
		N = length(Tr)*(length(Tr)-1) - length(Tr)*(length(Tr)-1)/2
		restrictions = getRestrictions.Deterministic(Tr, distanceThreshold.Tr, matrix(0, N), 1)
	}
	# estimate restrictions for the deterministic variables (if any)
	if (!is.null(distanceThreshold.C))
	{
    	if (is.null(restrictions))
		{
			N = length(Tr)*(length(Tr)-1) - length(Tr)*(length(Tr)-1)/2
			restrictions = matrix(0, N)
		}
    
		for (i in seq(ncol(X)))
		{
			restrictions = getRestrictions.Deterministic(X[,i], distanceThreshold.C, restrictions, 0)
		}
	}
	return(restrictions)
}
