#list of all methods:

#MatchGenoudStage1: normalizes treatment and confounding values to [0, 1]

#myPDF: takes a list of the values resulted from a Monte Carlo simulation and returns a probabilisty distribution for each unique value

#getCDFPerPair: estimates an array of all the CDFs for each potential matched pair of units

#getPDFPerPair: estimates an array of all the PDFs for each potential matched pair of units

#getRestrictionsTrProb: restricts the matching of units with probability to have the same treatment higher than a threshold

#getRestrictionsTrDeterm: restricts the matching of units with difference on their treatment less than a threshold

#getPercentileStats.Deterministic: estimate the percentiles of the random variable Z^p_{(u,v)} = Z^{p, U} - Z^{p, V}
#it is called when all variables are deterministic

#getPercentileStats.Probabilistic.Tr: calculated the difference on the percentiles of matched treated and control samples when only
#treatment is a random variable (i.e. the random variable is U_{(u,v)}Z^p_{(u,v)}, where Z^p_{(u,v)} = Z^{p, U} - Z^{p, V})

#getPercentileStats.Probabilistic.C: calculated the difference on the percentiles of matched treated and control samples when only
#confounding variables are probabilistic (i.e. the random variable is U_{(u,v)}Z^p_{(u,v)}, where Z^p_{(u,v)} = Z^{p, U} - Z^{p, V})

#getPercentileStats.Probabilistic: the same as above, but both treatment and confounders are random variables

#calculatePercentiles: helper function used by getPercentileStats.Probabilistic.Tr and getPercentileStats.Probabilistic.C

library(matrixStats)
library(mc2d)

source("probabilityEstimationFunctions.R", local = TRUE)
source("balanceEstimationFunctions.R", local = TRUE)
source("restrictionsEstimationFunctions.R", local = TRUE)
source("distancesEstimationFunctions.R", local = TRUE)

#normalizes all the variables to be from 0 to 1 and returns the result
#Args:
#	Tr: the values of treatment variable (if it is probabilistic, it contains the potential values)
#	X: the values of deterministic confounding variables
#	PX: the potential values of probabilistic confounding variables
#	Tr: the values of treatment variable (if it is probabilistic, it contains the potential values)
#	weights: a weight for each sample
MatchGenoudStage1 <- function(Tr=Tr, X=X, PX = PX, weights=weights)
{
	N  <- nrow(Tr)
  	xvars <- ncol(X)
  	pxvars <- nrow(PX)/N
  	Tr = (Tr - min(Tr))/(max(Tr) - min(Tr))
  
  	if (!is.null(X))
  	{
    	for (k in 1:xvars)
    	{
      	  	X[,k] = weights*X[,k]
      		X[,k] = (X[,k] - min(X[,k]))/(max(X[,k]) - min(X[,k]))
    	}
  	}
	if (!is.null(PX))
  	{
    	for (k in 1:pxvars)
    	{
      	  	thisPX = PX[((k-1)*N+1):(k*N),]
      		thisPX = matrix(rep(weights, each=ncol(thisPX)),nrow=N,ncol=ncol(thisPX),byrow=TRUE)*thisPX
     		PX[((k-1)*N+1):(k*N),] = (thisPX - min(thisPX))/(max(thisPX) - min(thisPX))
		}
	}
	ret <- list(Tr=Tr, X=X, PX=PX)
	return(ret)
} #end of MatchGenoudStage1


# Probabilistic Genetic Matching Algorithm
#Args:
#	Tr:	treatment variable, it can be either a N x 1 matrix for deterministic treatments or a list for probabilistic treatments
#		if Tr is a list it needs to contain a field values (a matrix of potential values) of size N x P, where N the number of samples and P the number of possible values and a field probs (a matrix of probabilities for each value) of size N x P where N the number of samples and P the number of possible values

#	X:	determnistic confounding variables of size N x N_x, where N the number of samples and N_x the number of deterministic confounding variables

#	pX: list of probabilistic confounding variables. It has a field 'values' (a matrix of potential values) of size N_p*N x P, where N the number of samples, P the number of possible values and N_p the number of confounding variables and a field 'probs' (a matrix of probabilities for each value) of size N*N_p x P. Thus, the first N values are for the first confounding variable, the next N for the 2nd etc.

#	M: the maximum number of times each unit can be used

#	weights: N x 1 matrix of weights per each sample

#	treatmentDiffThreshold: the minimum difference that treatment variables need to have in order to allow matching

# 	treatmentDiffProb: it is used only for probabilistic treatment variables. Units cannot be matched if their probability to have treatment difference smaller than treatmentDiffThreshold is larger than treatmentDiffProb

# 	cdd: the maximum difference that confounding variables need to have in order to allow matching

# 	pcdd: it is used only for probabilistic confounding variables. Units cannot be matched if their probability to have difference on a confounding variable smaller than cdd is larger than pcdd

#	epsilon: it is used only when the treatment is probabilistic. It is added to the treatment difference in order to avoid NaNs

#	pop.size, max.generations, wait.generations, hard.generation.limit, fit.func, MemoryMatrix, tolerance, min.weight, max.weight, Domains, print.level, project.path, loss, data.type.integer: have the same use as the original genetic matching

# Returns:
# 	mgr: probabilistic genetic matching result (resulted matched pairs, weights and parameters) 


GenMatchProbabilistic <- function(Tr, X=NULL, pX=NULL, M=1,
                                  weights=NULL,
                                  treatmentDiffThreshold = NULL,
                                  treatmentDiffProb = 0.1,
                                  cdd=NULL,
                                  pcdd = 0.1,
                                  epsilon = 0.0001,
                                  pop.size = 100, max.generations=100,
                                  wait.generations=4, hard.generation.limit=FALSE,
                                  fit.func="qqmean.mean",
                                  MemoryMatrix=TRUE,
                                  tolerance=sqrt(.Machine$double.eps),
                                  min.weight=0,
                                  max.weight=1000,
                                  Domains=NULL,
                                  print.level=2,
                                  project.path=NULL,
                                  loss=1,
                                  data.type.integer=FALSE, ...)
{

	requireNamespace("rgenoud")

  	#if Tr is a list, then Tr is probabilisitic
  	if (is.list(Tr))
  	{
    	#if Tr is a list it needs to contain a field 'values' (a matrix of potential values) of size N x P, where N the number of samples and P the number of possible values
    	#and a field 'probs' (a matrix of probabilities for each value) of size N x P where N the number of samples and P the number of possible values
		if (!'probs' %in% names(Tr))
		{
			stop("Treatment does not contain a key 'probs'!")
		}
		if (!'values' %in% names(Tr))
		{
			stop("Treatment does not contain a key 'values'!")
		}

    	Tr.probs <- as.matrix(Tr$probs)
		Tr <- as.matrix(Tr$values)
		if (sum(dim(Tr)==dim(Tr.probs))!=2)
    	{
      	  	stop("There is not a probability for each treatment value!")
    	}
    	if(length(which(Tr.probs<0))||length(which(Tr.probs>1)))
      	  	stop("The probabilities on the treatment variables are not between 0 and 1")
		if((length(which(rowSums(Tr.probs)<0.99)))||(length(which(rowSums(Tr.probs)>1.01))))
			stop("The sum of the probabilities for each treatment value should be 1.")
	}
  	else
  	{
		Tr <- as.matrix(Tr)
  	}
  	N <- dim(Tr)[1]
  	if ((is.null(X))&&(is.null(pX)))
  	{
    	stop("Confounding Variables have not been provided!")
  	}
  
  	# 0: all covariates are deterministic
  	#1: all covariates are random variables
  	#2: there are both deterministic covariates and random variables
  	covariatesTypeFlag <- -1
  	nvars <-0
  	if (!is.null(X))
  	{
    	X  <- as.matrix(X)
    	nvars <- nvars + dim(X)[2]
    	covariatesTypeFlag <-0
    	if (dim(X)[1]!=N)
    	{
      	  stop("The number of confounding variables values does not match the number of samples!")
    	}
  	}
    
  	npvars <- 0
  	if (!is.null(pX))
  	{
  		if (!'probs' %in% names(pX))
  		{
  			stop("pX does not contain a key 'probs'!")
  		}
  		if (!'values' %in% names(pX))
  		{
  			stop("pX does not contain a key 'values'!")
  		}
    	X.values <- as.matrix(pX$values)
    	X.probs <- as.matrix(pX$probs)
    	if ((dim(X.values)[1]%%N)!=0)
    	{
      	  	print(dim(X.values)[1]%%N)
      		stop("The number of probability distributions for the confounding variables does not match the number of samples!")
    	}
    	npvars <- dim(X.values)[1]/N
    	if (covariatesTypeFlag==0)
			covariatesTypeFlag <- 2
    	else
			covariatesTypeFlag <- 1
    	if(length(which(X.probs<0))||length(which(X.probs>1)))
			stop("The probabilities on the confounding variables are not between 0 and 1")
		if((length(which(rowSums(X.probs)<0.99)))||(length(which(rowSums(X.probs)>1.01))))
			stop("The sum of the probabilities for each confounding value should be 1.")
	}
  

	if (is.null(weights))
  	{
		weights <- rep(1,nrow(Tr))
    	weights.flag <- FALSE
	} else {
		weights.flag <- TRUE
		weights <- as.double(weights)
		if( nrow(Tr) != length(weights))
 	   	{
		   stop("length(Tr) != length(weights)")
	   	}
	}

  	isna  <- sum(is.na(Tr)) + sum(is.na(weights))
  
  	if (covariatesTypeFlag==0)
		isna <- isna + sum(is.na(X))
	if (covariatesTypeFlag==1)
    	isna <- isna + sum(is.na(X.values))
	if (covariatesTypeFlag==2)
    	isna <- isna + sum(is.na(X.values)) + sum(is.na(X))
	if (isna!=0)
  	{
    	stop("GenMatch(): input includes NAs")
    	return(invisible(NULL))
	}

  	#check inputs
  
  	if ((treatmentDiffProb < 0)||(treatmentDiffProb > 1))
  	{
		warning("User set 'treatmentDiffProb' to less than 0 or more than 1.  Resetting to the default which is 0.1.")
    	treatmentDiffProb <- 0.1
	}
  
  	if ((pcdd < 0)||(pcdd > 1))
  	{
		warning("User set 'pcdd' to less than 0 or more than 1.  Resetting to the default which is 0.1.")
    	treatmentDiffProb <- 0.1
	}

  	if (pop.size < 0 | pop.size!=round(pop.size) )
  	{
    	warning("User set 'pop.size' to an illegal value.  Resetting to the default which is 100.")
    	pop.size <- 100
  	}
  	if (max.generations < 0 | max.generations!=round(max.generations) )
  	{
    	warning("User set 'max.generations' to an illegal value.  Resetting to the default which is 100.")
    	max.generations <-100
	}
  	if (wait.generations < 0 | wait.generations!=round(wait.generations) )
  	{
    	warning("User set 'wait.generations' to an illegal value.  Resetting to the default which is 4.")
    	wait.generations <- 4
  	}
  	if (hard.generation.limit != 0 & hard.generation.limit !=1 )
  	{
    	warning("User set 'hard.generation.limit' to an illegal value.  Resetting to the default which is FALSE.")
    	hard.generation.limit <- FALSE
  	}
  	if (data.type.integer != 0 & data.type.integer !=1 )
  	{
    	warning("User set 'data.type.integer' to an illegal value.  Resetting to the default which is TRUE.")
    	data.type.integer <- TRUE
  	}
  	if (MemoryMatrix != 0 & MemoryMatrix !=1 )
  	{
    	warning("User set 'MemoryMatrix' to an illegal value.  Resetting to the default which is TRUE.")
    	MemoryMatrix <- TRUE
  	}
  
  	if (min.weight < 0)
  	{
    	warning("User set 'min.weight' to an illegal value.  Resetting to the default which is 0.")
    	min.weight <- 0
  	}
  	if (max.weight < 0)
  	{
    	warning("User set 'max.weight' to an illegal value.  Resetting to the default which is 1000.")
    	max.weight <- 1000
  	}
  	if (print.level != 0 & print.level !=1 & print.level !=2 & print.level !=3)
  	{
    	warning("User set 'print.level' to an illegal value.  Resetting to the default which is 2.")
    	print.level <- 2
  	}
  
  	##from Match()
  	if (tolerance < 0)
  	{
    	warning("User set 'tolerance' to less than 0.  Resetting to the default which is 0.00001.")
    	tolerance <- 0.00001
  	}
	if (M < 1)
  	{
    	warning("User set 'M' to less than 1.  Resetting to the default which is 1.")
    	M <- 1
  	}
  	if ( M!=round(M) )
  	{
    	warning("User set 'M' to an illegal value.  Resetting to the default which is 1.")
    	M <- 1
  	}
  

  	#print warning if pop.size, max.generations and wait.generations are all set to their original values
  	if(pop.size==100 & max.generations==100 & wait.generations==4)
  	{
    	warning("The key tuning parameters for optimization were are all left at their default values.  The 'pop.size' option in particular should probably be increased for optimal results.  For details please see the help page and http://sekhon.berkeley.edu/papers/MatchingJSS.pdf")
  	}

  	#loss function
  	if (is.double(loss))
  	{
    	if (loss==1)  
		{
			loss.func=sort
      		lexical=nvars+npvars 
		} 
		else if(loss==2)
		{
      	  loss.func=min
		  lexical=0
	  	} 
		else
		{
			stop("unknown loss function")
		}
	} 
	else if (is.function(loss)) {
    loss.func=loss
    lexical=1
	} 
	else 
	{
    	stop("unknown loss function")
	}

	#set lexical for fit.func
  	if (fit.func=="qqmean.max" | fit.func=="qqmedian.max" | fit.func=="qqmax.max")
	{
    	lexical=nvars+npvars
  	} 
	else if (fit.func!="qqmean.mean" & fit.func!="qqmean.max" &
             fit.func!="qqmedian.median" & fit.func!="qqmedian.max"
             & fit.func!="pvals") 
	{
		stop("invalid 'fit.func' argument")
	} 
	else
	{
    	lexical = 0
	}

	isunix  <- .Platform$OS.type=="unix"
	if (is.null(project.path))
	{
    	if (print.level < 3 & isunix)
		{
			project.path="/dev/null"
		} 
		else 
		{
			project.path=paste(tempdir(),"/genoud.pro",sep="")

			#work around for rgenoud bug
			#if (print.level==3)
			#print.level <- 2
		}
	}

  	balancevars <- nvars + npvars
  	starting.values=rep(1,balancevars)
  	if (is.null(Domains))
  	{
    	Domains <- matrix(min.weight, nrow=nvars+npvars, ncol=2)
    	Domains[,2] <- max.weight
  	} 
	else
	{
    	indx <- (starting.values < Domains[,1]) | (starting.values > Domains[,2])
    	starting.values[indx] <- round( (Domains[indx,1]+Domains[indx,2])/2 )
  	}
  

  	#if all covariates are deterministic PX=NULL
  	if (!covariatesTypeFlag)
  	{
    	s1 <- MatchGenoudStage1(Tr=Tr, X=X, PX=NULL, weights=weights);
    	s1.X <- s1$X
    	s1.X.values <- NULL
  	}
  	#if all covariates are probabilistic X=NULL
  	else if(covariatesTypeFlag==1)
  	{
    	s1 <- MatchGenoudStage1(Tr=Tr, X=NULL, PX=X.values, weights=weights);
		s1.X <- NULL
		s1.X.values <- s1$PX
	}
	else
  	{
    	s1 <- MatchGenoudStage1(Tr=Tr, X=X, PX=X.values, weights=weights);
    	s1.X <- s1$X
    	s1.X.values <- s1$PX
  	}
  
  	s1.treatmentDiffThreshold <-treatmentDiffThreshold
  	s1.Tr <- s1$Tr
  	rm(s1)
  
  
  	PX <- NULL
  	PairPDF <-NULL
  	restrictions <- NULL
  	if (covariatesTypeFlag)
  	{
    	#number of probabilistic confounding variables
    	PVarsNum <- nrow(s1.X.values)/N
		#PX: matrix with the CDFs each probabilistic confounding variable for each sample
		PX <- c()
		#for each probabilistic confounding variable:
		for (i in seq(PVarsNum))
		{
			#the PDFs for the ith variable are the positions: (i-1)*N+1 to (i-1)*N+1+N-1
			indx <- (i-1)*N+1
			#if the treatment is also probabilistic put the value 2 on the last flag argument of getCDFPerPair
			#getCDFPerPair estimates the CDF for each pair of samples (see more details on function comments)
			if (dim(s1.Tr)[2]>1)
			{
				PX <- rbind(PX, getCDFPerPair(s1.X.values[indx:(indx+N-1),], X.probs[indx:(indx+N-1),], s1.Tr, Tr.probs, 2, epsilon))
				PairPDF <- getPDFPerPair(s1.X.values[indx:(indx+N-1),], X.probs[indx:(indx+N-1),], s1.Tr, Tr.probs, 2, epsilon)
			}
			else
			{
				PX <- rbind(PX, getCDFPerPair(s1.X.values[indx:(indx+N-1),], X.probs[indx:(indx+N-1),], s1.Tr, NULL, 1, epsilon))
				PairPDF <- getPDFPerPair(s1.X.values[indx:(indx+N-1),], X.probs[indx:(indx+N-1),], s1.Tr, NULL, 1, epsilon)
			}
		}
		PX = PX[,seq(10, 99, 10)]
	  }
	  #if only the treatment is probabilistic estimate the CDF of variable 1/(Tr[i]-Tr[j]). If everything is deterministic, there is no reason to do anything...
	  else if(dim(s1.Tr)[2]>1)
	  {
		  PairPDF <- getPDFPerPair(NULL, NULL, s1.Tr, Tr.probs, 0, epsilon)
		  s1.Tr <- getCDFPerPair(NULL, NULL, s1.Tr, Tr.probs, 0, epsilon)
		  s1.Tr = s1.Tr[,seq(10, 99, 10)]
	  }

	  restrictions = getRestrictions(s1.Tr, PX, s1.X, treatmentDiffThreshold, treatmentDiffProb, cdd, pcdd, epsilon)
	  if(dim(s1.Tr)[2]>1)
	  {
		  if (length(which(restrictions==1))>0.95*dim(s1.Tr)[1])
		  {
			  stop("The restrictions specified do not allow for matching.")
		  }
	  }
	  else if(covariatesTypeFlag)
	  {
		  if (length(which(restrictions==1))>0.9*dim(PX)[1]/PVarsNum)
		  {
			  stop("The restrictions specified do not allow for matching.")
		  }
	  }

  
	  dists <- getDistances(s1.X, PX, s1.Tr, N, nvars, npvars)
  
	  genoudfunc  <- function(x)
	  {
    
    	  wmatrix <- diag(x, nrow=nvars+npvars)
		  if ( min(eigen(wmatrix, symmetric=TRUE, only.values=TRUE)$values) < tolerance )
			  wmatrix <- wmatrix + diag(nvars+npvars)*tolerance

		  ww <- chol(wmatrix)
    
		  rr <- getPairs(dists, diag(ww), restrictions, M, N)
		  for (i in 1:length(rr[,1])){
			  rr[i,1:2] = sort(rr[i,1:2])
		  }
		  rr = unique( rr[ ,1:3] )
		  rmIndx = which(rowSums(rr)!=0)
		  if (length(rmIndx)>0)
			  rr = rr[rmIndx,]
		  
		  GenBalanceQQ.internal <- function(rr, X, XP, Tr, PairPDF, nvars, npvars, N, covariatesTypeFlag, summarystat="mean", summaryfunc="mean")
		  {
			  
			  if (covariatesTypeFlag)
			  {
				  if(dim(Tr)[2]>1)
					  percentilesPerVar <- getPercentileStats.Probabilistic(Tr, rr, X, PairPDF, N, nvars, npvars)
				  else
					  percentilesPerVar <- getPercentileStats.Probabilistic.C(Tr, rr, X, PairPDF, N, nvars, npvars)
			  }
			  else
			  {
				  if(dim(Tr)[2]>1)
					  percentilesPerVar <- getPercentileStats.Probabilistic.Tr(PairPDF, rr, X, N, nvars)
				  else
					  percentilesPerVar <- getPercentileStats.Deterministic(Tr, rr, X, N, nvars)
			  }
      
			  qqSummary <- probabilisticQQStats(percentilesPerVar)
			  # print(qqSummary$meanDiff)
			  qqSummaryStatistic <- c()
      
			  if (summarystat=="median")
			  {
				  qqSummaryStatistic <- qqSummary$mediandiff
			  } else if (summarystat=="max")  {
				  qqSummaryStatistic <- qqSummary$maxDiff
			  } else {
				  qqSummaryStatistic <- qqSummary$meanDiff
			  }
      
			  
			  if (summaryfunc=="median")
			  {
				  return(median(qqSummaryStatistic))
			  } else if (summaryfunc=="max")  {
				  return(max(qqSummaryStatistic))
			  } else {
				  return(mean(qqSummaryStatistic))
			  }
		  } #end of GenBalanceQQ.internal
    

		  if (fit.func=="qqmean.mean") {
			  a <- GenBalanceQQ.internal(rr=rr, X=s1.X, XP=XP, Tr=s1.Tr, PairPDF, nvars, npvars, N = N, covariatesTypeFlag=covariatesTypeFlag, summarystat="mean", summaryfunc="mean")
			  return(a)
		  } else if (fit.func=="qqmean.max") {
			  a <- GenBalanceQQ.internal(rr=rr, X=s1.X, XP=XP, Tr=s1.Tr, PairPDF, nvars, npvars, N = N, covariatesTypeFlag=covariatesTypeFlag, summarystat="mean", summaryfunc="max")
			  return(a)
		  } else if (fit.func=="qqmax.mean") {
			  a <- GenBalanceQQ.internal(rr=rr, X=s1.X, XP=XP, Tr=s1.Tr, PairPDF, nvars, npvars, N = N, covariatesTypeFlag=covariatesTypeFlag, summarystat="max", summaryfunc="mean")
			  return(a)
		  } else if (fit.func=="qqmax.max") {
			  a <- GenBalanceQQ.internal(rr=rr, X=s1.X, XP=XP, Tr=s1.Tr, PairPDF, nvars, npvars, N = N, covariatesTypeFlag=covariatesTypeFlag, summarystat="max", summaryfunc="max")
			  return(a)
		  } else if (fit.func=="qqmedian.median") {
			  a <- GenBalanceQQ.internal(rr=rr, X=s1.X, XP=XP, Tr=s1.Tr, PairPDF, nvars, npvars, N = N, covariatesTypeFlag=covariatesTypeFlag, summarystat="median", summaryfunc="median")
			  return(a)
		  } else if (fit.func=="qqmedian.max") {
			  a <- GenBalanceQQ.internal(rr=rr, X=s1.X, XP=XP, Tr=s1.Tr, PairPDF, nvars, npvars, N = N, covariatesTypeFlag=covariatesTypeFlag, summarystat="median", summaryfunc="max")
			  return(a)
		  }
	  } #end genoudfunc

  
	  cl.genoud <- FALSE
 
	  do.max <- FALSE
	  rr <- rgenoud::genoud(genoudfunc, nvars=nvars+npvars, starting.values=starting.values,
                        pop.size=pop.size, max.generations=max.generations,
                        wait.generations=wait.generations, hard.generation.limit=hard.generation.limit,
                        Domains=Domains,
                        MemoryMatrix=MemoryMatrix,
                        max=do.max, gradient.check=FALSE, data.type.int=data.type.integer,
                        hessian=FALSE,
                        BFGS=FALSE, project.path=project.path, print.level=print.level,
                        lexical=lexical,
                        cluster=cl.genoud,
                        balance=TRUE,
                        ...)
						
	  wmatrix <- diag(rr$par, nrow=nvars+npvars)
	  if ( min(eigen(wmatrix, symmetric=TRUE, only.values=TRUE)$values) < tolerance )
		  wmatrix <- wmatrix + diag(nvars+npvars)*tolerance
	  ww <- chol(wmatrix)

	  mout <- getPairs(dists, diag(ww), restrictions, M, N)
	  for (i in 1:length(mout[,1])){
		  mout[i,1:2] = sort(mout[i,1:2])
	  }
	  mout = unique( mout[ ,1:3] )
	  rmIndx = which(rowSums(mout)!=0)
	  mout = mout[rmIndx,]
	  gmr <- list(value=rr$value, par=rr$par, Weight.matrix=wmatrix, matches=mout, ecaliper=NULL)
	  class(gmr) <- "GenMatch"
	  return(gmr)
} #end of GenMatch


# find the best match for each unit based on a weighted mahalanobis distance
# Args:
# 	dists: a matrix of size N*(N-1)/2 x (vars+pvars) with the distances between each pair of units (with N the number of units, vars the number of confounding variables and pvars the number of probabilistic confounding variables)
# 	ww: an array of weights of size (vars+pvars) x 1 (i.e. one weight for each confounding variable)
# 	restrictions: a binary vector of size N*(N-1)/2 x 1. '1' indicates pairs of units that are not allowed to be matched
# 	M: the maximum number of times each unit can be used
# 	N: the number of samples
# Returns:
# 	matches: an array with the matched pairs. The first two columns correspond to the pair of matched units and the 3rd to the weights of each pair
getPairs <- function(dists, ww, restrictions, M, N)
{
	MAXDIST = 100000000
  	vNum = length(ww)
  	matches = cbind()
  	# multiply with dists with ww
  	for (i in 1:vNum)
  	{
    	dists[,i] = ww[i]^2 * dists[,i]
  	}
  	dists = rowSums(dists)
	dists = sqrt(dists)
	dists[which(restrictions==1)] = MAXDIST
  	distsMatrix = matrix(0, N, N)
  	indx = 1
  	for (i in 1:(N-1))
	{
    	distsMatrix[i,(i+1):N] = dists[indx:(indx-1+N-i)]
    	distsMatrix[(i+1):N,i] = dists[indx:(indx-1+N-i)]
    	indx = indx + N-i
  	}
  	diag(distsMatrix) = rep(MAXDIST, N)
  	distsMatrix[which(is.nan(distsMatrix))] = MAXDIST
  	# keep how many times each unit has been used at vector elUseTimes. If a unit has been used M times cannot be used again. 
  	elUseTimes = rep(0, N)
  	for (i in 1:N)
  	{
    	if (elUseTimes[i]<M)
    	{
			minD = min(distsMatrix[i,])
			if (minD < MAXDIST)
		   	{
        		minDpos = which.min(distsMatrix[i,])
        		matches = rbind(matches, c(i, minDpos, 0))
        		elUseTimes[i] = elUseTimes[i] + 1
        		elUseTimes[minDpos] = elUseTimes[minDpos] + 1
        		distsMatrix[i,minDpos] = MAXDIST
        		distsMatrix[minDpos,i] = MAXDIST
        		if (elUseTimes[i]==M)
				{
					distsMatrix[i,] = MAXDIST
          		  	distsMatrix[,i] = MAXDIST
        		}
				if (elUseTimes[minDpos]==M)
				{
					distsMatrix[minDpos,] = MAXDIST
					distsMatrix[,minDpos] = MAXDIST
				}
			}
		}
	}
	matches = getPairFrequency(matches)
  	return(matches)
}

# estimates a weight for each pair of units
# Args:
# 	matches: array of matched pairs. This is a two-column array (each row corresponds to a pair of units)
# Returns:
# 	matches: a 3rd column with a weight for each pair has been added to the initial matches array
getPairFrequency <- function(matches)
{
	matchesAsVector = c(matches[,1], matches[,2])
  	maxEl = max(matchesAsVector)
  	matchesFreq = rep(0, maxEl)
  	for (i in seq(maxEl))
	{
    	matchesFreq[i] = length(which(matchesAsVector==i))
  	}
  	matches[,3] = 0.5*(1/matchesFreq[matches[,1]] + 1/matchesFreq[matches[,2]])
  	return(matches)
}