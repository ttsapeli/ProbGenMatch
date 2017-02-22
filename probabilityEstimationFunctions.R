# Estimates the probability of a value to occur based on the output of a monte carlo simulation
#Args: 
#	values: a list of all the values resulted from a monte carlo simulation
#	uniqueValues: a list of the unique values
#Return:
#	probs: a vector of size equal to uniqueValues size, with the corresponding probabilities for each value
myPDF <- function(values, uniqueValues)
{
	N = length(values)
  	probs = rep(0, length(uniqueValues))
  	for (i in seq(length(uniqueValues)))
  	{
    	probs[i] = length(which(values==uniqueValues[i]))/N  
  	}
  	return(probs)
}


#Args:
#	valuesC: potential values of a probabilistic confounding variable, size: N x P, where P the number of potential values and N the number of samples
#	probsC: probabilisties for each potential value of a probabilistic confounding variable, size: N x P, where P the number of potential values and N the number of samples
#	valuesTr: potential values of treatment variable, size: N x 1,where P the number of potential values and N the number of samples
#	probsTr: probabilisties for each potential value of treatment variable, size: N x 1, where P the number of potential values and N the number of samples
#	flag: 0 if only the treatment variable is probabilistic, 1 if only the confounding variable is probabilistic and 2 if both are. 
#	epsilon: if the treatment is uncertain, we add to the treatment difference a small positive value epsilon in order to avoid NaNs
#Returns: 
#	CDFPerPair: for each pair of samples the CDF of their difference is estimated. 
#		Case 1: flag = 0. For each pair (i, j), with i < j we estimate the CDF of variable 1/(P[i]-P[j]) and we store it at the row (i-1)*N - (i)*(i-1)/2 + j-i
#		Case 2: flag = 1. For each pair (i, j), with i < j we estimate the CDF of variable (X[i]-X[j]) and we store it at the row (i-1)*N - (i)*(i-1)/2 + j-i
#		Case 3: flag = 2. For each pair (i, j), with i < j we estimate the CDF of variable (X[i]-X[j])/(P[i]-P[j]) and we store it at the row (i-1)*N - (i)*(i-1)/2 + j-i
getCDFPerPair <- function(valuesC, probsC, valuesTr, probsTr, flag, epsilon)
{
	SamplesNum = dim(valuesTr)[1]
  	CDFPerPair <- matrix(0, (SamplesNum-1)*SamplesNum/2, 99)
  	counter <- 1
  	sample1ProbDistC<-NULL
  	sample1ProbDistTr<-NULL
  	for (i in seq(SamplesNum-1))
  	{
    	if(flag>0)
      	  sample1ProbDistC <- mcstoc(rempiricalD, values=valuesC[i,], prob=probsC[i,], nsv=10000)
    	if(flag!=1)
      	  sample1ProbDistTr <- mcstoc(rempiricalD, values=valuesTr[i,], prob=probsTr[i,], nsv=10000)
    
    	for (j in seq(i+1, SamplesNum))
    	{
			if (flag==2)
			{
        		sample2ProbDistC <- mcstoc(rempiricalD, values=valuesC[j,], prob=probsC[j,], nsv=10000)
        		sample2ProbDistTr <- mcstoc(rempiricalD, values=valuesTr[j,]+epsilon, prob=probsTr[j,], nsv=10000)
        		valuesDiff <- abs((sample1ProbDistC-sample2ProbDistC)/(sample1ProbDistTr-sample2ProbDistTr))
      	  	}
			else if(flag==1)
			{
				sample2ProbDistC <- mcstoc(rempiricalD, values=valuesC[j,], prob=probsC[j,], nsv=10000)
        		valuesDiff <- abs(sample1ProbDistC-sample2ProbDistC)
			}
      	  	else
      	  	{
				sample2ProbDistTr <- mcstoc(rempiricalD, values=valuesTr[j,]+epsilon, prob=probsTr[j,], nsv=10000)
				valuesDiff <- 1/abs(sample1ProbDistTr-sample2ProbDistTr)
			}
			diffSummary <- summary(valuesDiff, probs=seq(0, 1, 0.01))
			CDFPerPair[counter,] = diffSummary$node[4:102]
			counter = counter + 1
		}
	}
	return(CDFPerPair)
}


#Args:
#	valuesC: potential values of a probabilistic confounding variable, size: N x P, where P the number of potential values and N the number of samples
#	probsC: probabilisties for each potential value of a probabilistic confounding variable, size: N x P, where P the number of potential values and N the number of samples
#	valuesTr: potential values of treatment variable, size: N x 1,where P the number of potential values and N the number of samples
#	probsTr: probabilisties for each potential value of treatment variable, size: N x 1, where P the number of potential values and N the number of samples
#	flag: 0 if only the treatment variable is probabilistic, 1 if only the confounding variable is probabilistic and 2 if both are. 
#	epsilon: if the treatment is uncertain, we add to the treatment difference a small positive value epsilon in order to avoid NaNs
#Returns: 
#	PDFPerPair: for each pair of samples the PDF of their difference is estimated. 
#		Case 1: flag = 0. For each pair (i, j), with i < j we estimate the CDF of variable 1/(P[i]-P[j]) and we store it at the row (i-1)*N - (i)*(i-1)/2 + j-i
#		Case 2: flag = 1. For each pair (i, j), with i < j we estimate the CDF of variable (X[i]-X[j]) and we store it at the row (i-1)*N - (i)*(i-1)/2 + j-i
#		Case 3: flag = 2. For each pair (i, j), with i < j we estimate the CDF of variable (X[i]-X[j])/(P[i]-P[j]) and we store it at the row (i-1)*N - (i)*(i-1)/2 + j-i
getPDFPerPair <- function(valuesC, probsC, valuesTr, probsTr, flag, epsilon)
{
	SamplesNum = dim(valuesTr)[1]
	numOfProbs = 0
 	uniqueValues = c()
  	if(flag!=1)
  	{
    	uniqueValuesT = c()
    	for (i in valuesTr[1,])
    	{
			for (j in (valuesTr[1,]+epsilon))
			{
				uniqueValuesT = c(uniqueValuesT, 1/abs(i-j))
			}
		}
		uniqueValuesT = unique(uniqueValuesT)
		if (flag==0)
			uniqueValues = uniqueValuesT
	}
	
	if(flag!=0)
	{
		uniqueValuesC = c()
    	for (i in valuesC[1,])
    	{
			for (j in valuesC[1,])
			{
				uniqueValuesC = c(uniqueValuesC, abs(i-j))
			}
		}
		uniqueValuesC = unique(uniqueValuesC)
		if (flag==1)
			uniqueValues = uniqueValuesC
	}
	
	if (flag==2)
	{
		for (i in uniqueValuesT)
		{
			for (j in uniqueValuesC)
				uniqueValues = c(uniqueValues,j/i)
		}
		uniqueValues = unique(uniqueValues)
	}
	numOfProbs = length(uniqueValues)
	PDFPerPair <- matrix(0, (SamplesNum-1)*SamplesNum/2, numOfProbs)
	counter <- 1
	sample1ProbDistC<-NULL
	sample1ProbDistTr<-NULL
  	for (i in seq(SamplesNum-1))
  	{
    	if(flag>0)
			sample1ProbDistC <- mcstoc(rempiricalD, values=valuesC[i,], prob=probsC[i,], nsv=1000)
    	if(flag!=1)
      	  	sample1ProbDistTr <- mcstoc(rempiricalD, values=valuesTr[i,], prob=probsTr[i,], nsv=1000)
		for (j in seq(i+1, SamplesNum))
    	{
			if (flag==2)
			{
				sample2ProbDistC <- mcstoc(rempiricalD, values=valuesC[j,], prob=probsC[j,], nsv=1000)
				sample2ProbDistTr <- mcstoc(rempiricalD, values=valuesTr[j,]+epsilon, prob=probsTr[j,], nsv=1000)
				valuesDiff <- abs((sample1ProbDistC-sample2ProbDistC)/(sample1ProbDistTr-sample2ProbDistTr))
			}
			else if(flag==1)
			{
				sample2ProbDistC <- mcstoc(rempiricalD, values=valuesC[j,], prob=probsC[j,], nsv=1000)
				valuesDiff <- abs(sample1ProbDistC-sample2ProbDistC)
			}
			else
			{
				sample2ProbDistTr <- mcstoc(rempiricalD, values=valuesTr[j,]+epsilon, prob=probsTr[j,], nsv=1000)
				valuesDiff <- 1/abs(sample1ProbDistTr-sample2ProbDistTr)
			}
			PDFPerPair[counter,] = myPDF(cbind(valuesDiff), uniqueValues)
			counter = counter + 1
		}
	}
	return(list(values=uniqueValues, probs=PDFPerPair))
}

