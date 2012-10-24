# (C) Copyright David Hastie and Silvia Liverani, 2012.

# DiPBaC++ is free software; you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License as published by the Free Software
# Foundation; either version 3 of the License, or (at your option) any later
# version.

# DiPBaC++ is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with DiPBaC++ in the documentation directory. If not, see
# <http://www.gnu.org/licenses/>.

# The external linear algebra library Eigen, parts of which are included  in the
# lib directory is released under the LGPL3+ licence. See comments in file headers
# for details.

# The Boost C++ header library, parts of which are included in the  lib directory
# is released under the Boost Software Licence, Version 1.0, a copy  of which is
# included in the documentation directory.

profRegr<-function(covNames, fixedEffectsNames=-999, outcome="outcome", outcomeT=-999, data, output="output", hyper, predict, nSweeps=1000, nBurn=1000, nProgress=500, nFilter=1, nClusInit, seed, yModel="Bernoulli", xModel="Discrete", sampler="SliceDependent", alpha=-1, excludeY, extraYVar, varSelect, entropy){

	nCovariates<-length(covNames)
	
	if (!is.data.frame(data)) stop("Input data must be a data.frame with outcome, covariates and fixed effect names as column names.")

	# open file to write output
	fileName<-paste(output,"_input.txt",sep="")
	# make big data matrix with outcome, covariates and fixed effects	
	# outcome
	dataMatrix<-data$outcome

	# covariates
	covIndeces<-vector()
	for (i in 1:nCovariates){
		tmpIndex<-which(colnames(data)==covNames[i])
		if (length(tmpIndex)==0) stop("ERROR: covariate names in data.frame provided do not correspond to list of covariates for profile regression")
		covIndeces<-append(covIndeces,tmpIndex)
	}
	dataMatrix<-cbind(dataMatrix,data[,covIndeces])

	# fixed effects
	if (fixedEffectsNames!=-999) {
		nFixedEffects<-length(fixedEffectsNames)
		FEIndeces<-vector()
		for (i in 1:nFixedEffects){
			tmpIndex<-which(colnames(data)==fixedEffectsNames[i])
			if (length(tmpIndex)==0) stop("ERROR: fixed effects names in data.frame provided do not correspond to list of fixed effects for profile regression")
			FEIndeces<-append(FEIndeces,tmpIndex)
		}
		dataMatrix<-cbind(dataMatrix,data[,FEIndeces])
	} else {
		nFixedEffects<-0
	}

	#  extra outcome data
	if (yModel=="Poisson"||yModel=="Binomial") {
		if(outcomeT==-999){
			stop ("It is required to set outcomeT for Poisson (offset) or Binomial (number of trials) outcome.")
		} else {
			dataMatrix<-cbind(dataMatrix,outcomeT)	
		}
	} else {
		if(outcomeT!=-999) stop ("It is only required to set outcomeT for Poisson and Binomial outcome.")
	}		

	# print number of subjects
	write(as.character(dim(dataMatrix)[1]), fileName,ncolumns=1)
	# print number of covariates and their names
	write(as.character(nCovariates),fileName,append=T,ncolumns=1)
	write(t(covNames), fileName,append=T,ncolumns=1)
	# print number of fixed effects and their names
	write(nFixedEffects, fileName,append=T,ncolumns=1)
	if (nFixedEffects>0){
		write(t(fixedEffectsNames), fileName,append=T,ncolumns=1)
	}

	is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

	if (yModel=="Categorical"||yModel=="Bernoulli"){
		outcomeFactor<-as.factor(dataMatrix[,1])
		yLevels<-length(levels(outcomeFactor))	
		if (yModel=="Categorical") write(yLevels,fileName,append=T,ncolumns=length(yLevels))
		if (!(min(outcome)==0&&max(outcome)==(yLevels-1)&&sum(!is.wholenumber(outcome))==0)) {
			print("Recoding of the outcome as follows")
			tmpLevels<-levels(outcomeFactor)
			print(paste("Replacing level ",levels(outcomeFactor)," with ",c(0:(yLevels-1)),sep=""))
			levels(outcomeFactor)<-c(0:(yLevels-1))
			dataMatrix[,1]<-outcomeFactor
		}
	}

	if (xModel=="Discrete"){
		xLevels<-vector()
		for (k in 1:nCovariates){
			tmpCov<-dataMatrix[,(1+k)]
			xLevels[k]<-length(levels(as.factor(tmpCov)))	
			if (!(min(tmpCov)==0&&max(tmpCov)==(xLevels[k]-1)&&sum(!is.wholenumber(tmpCov))==0)) {
				print(paste("Recoding of covariate number ",colnames(dataMatrix)[k+1]," as follows",sep=""))
				tmpCovFactor<-as.factor(tmpCov)
				tmpLevels<-levels(tmpCovFactor)
				print(paste("Replacing level ",levels(tmpCovFactor)," with ",c(0:(xLevels[k]-1)),sep=""))
				levels(tmpCovFactor)<-c(0:(xLevels[k]-1))	
				dataMatrix[,(1+k)]<-tmpCovFactor
			}
		}
		write(xLevels,fileName,append=T,ncolumns=length(xLevels))
	}

	write(t(dataMatrix), fileName,append=T,ncolumns=dim(dataMatrix)[2])

	# other checks to ensure that there are no errors when calling the program
	if (xModel!="Discrete"&xModel!="Continuous") stop("Error in xModel")
	if (yModel!="Poisson"&yModel!="Binomial"&yModel!="Bernoulli"&yModel!="Normal"&yModel!="Categorical") stop("Error in xModel")

	inputString<-paste("--xModel=",xModel," --yModel=",yModel," --input=",fileName," --output=",output,sep="")

	if (!missing(alpha)) inputString<-paste(inputString," --alpha=",alpha,sep="")
	if (!missing(sampler)) inputString<-paste(inputString," --sampler=",sampler,sep="")
	if (!missing(hyper)) inputString<-paste(inputString," --hyper=",hyper,sep="")
	if (!missing(predict)) inputString<-paste(inputString," --predict=",predict,sep="")
	if (!missing(nSweeps)) inputString<-paste(inputString," --nSweeps=",nSweeps,sep="")
	if (!missing(nBurn)) inputString<-paste(inputString," --nBurn=",nBurn,sep="")
	if (!missing(nProgress)) inputString<-paste(inputString," --nProgress=",nProgress,sep="")
	if (!missing(nFilter)) inputString<-paste(inputString," --nFilter=",nFilter,sep="")
	if (!missing(nClusInit)) inputString<-paste(inputString," --nClusInit=",nClusInit,sep="")
	if (!missing(seed)) inputString<-paste(inputString," --seed=",seed,sep="")
	if (!missing(excludeY)) inputString<-paste(inputString," --excludeY",sep="")
	if (!missing(extraYVar)) inputString<-paste(inputString," --extraYVar",sep="")
	if (!missing(varSelect)) inputString<-paste(inputString," --varSelect=",varSelect,sep="")
	if (!missing(entropy)) inputString<-paste(inputString," --entropy",sep="")

	.Call('profRegr', inputString, PACKAGE = 'DiPBaC')
}


profRegrOld<-function(nCovariates, covNames, nFixedEffects, fixedEffectsNames, xLevels, yLevel, data="./input.txt", output="./output", hyper, predict, nSweeps=1000, nBurn=1000, nProgress=500, nFilter=1, nClusInit, seed, yModel="Bernoulli", xModel="Discrete", sampler="SliceDependent", alpha=-1, excludeY, extraYVar, varSelect, entropy){
	
	# prepare input file
	nColsData<-nCovariates+nFixedEffects+1
	if (yModel=="Poisson"||yModel=="Binomial") nColsData<-nColsData+1
	fileName<-paste(output,"_input.txt",sep="")
	dataMatrix<-as.matrix(read.table(data))
	# print number of subjects
	write(as.character(dim(dataMatrix)[1]), fileName,ncolumns=1)
	# print number of covariates and their names
	# check that the dimension of the data matrix against other inputs given	
	if (nColsData==dim(dataMatrix)[2]) {
		write(as.character(nCovariates), fileName,append=T,ncolumns=1)
	} else {
		stop("The number of covariates, fixed effects and other columns in data does not match the other inputs.")
	}
	write(t(covNames), fileName,append=T,ncolumns=1)
	# print number of fixed effects and their names
	if (nFixedEffects>0){
		write(nFixedEffects, fileName,append=T,ncolumns=1)
		write(t(fixedEffectsNames), fileName,append=T,ncolumns=1)
	} else {
		write(as.character(0), fileName,append=T,ncolumns=1)
	}
	if (xModel=="Categorical") write(yLevels,fileName,append=T,ncolumns=length(yLevels))
	if (xModel=="Discrete") write(xLevels,fileName,append=T,ncolumns=length(xLevels))
	write(t(dataMatrix), fileName,append=T,ncolumns=nColsData)

	# other checks to ensure that there are no errors when calling the program
	if (xModel!="Discrete"&xModel!="Continuous") stop("Error in xModel")
	if (yModel!="Poisson"&yModel!="Binomial"&yModel!="Bernoulli"&yModel!="Normal"&yModel!="Categorical") stop("Error in xModel")

	inputString<-paste("--xModel=",xModel," --yModel=",yModel," --input=",fileName," --output=",output,sep="")

	if (!missing(alpha)) inputString<-paste(inputString," --alpha=",alpha,sep="")
	if (!missing(sampler)) inputString<-paste(inputString," --sampler=",sampler,sep="")
	if (!missing(hyper)) inputString<-paste(inputString," --hyper=",hyper,sep="")
	if (!missing(predict)) inputString<-paste(inputString," --predict=",predict,sep="")
	if (!missing(nSweeps)) inputString<-paste(inputString," --nSweeps=",nSweeps,sep="")
	if (!missing(nBurn)) inputString<-paste(inputString," --nBurn=",nBurn,sep="")
	if (!missing(nProgress)) inputString<-paste(inputString," --nProgress=",nProgress,sep="")
	if (!missing(nFilter)) inputString<-paste(inputString," --nFilter=",nFilter,sep="")
	if (!missing(nClusInit)) inputString<-paste(inputString," --nClusInit=",nClusInit,sep="")
	if (!missing(seed)) inputString<-paste(inputString," --seed=",seed,sep="")
	if (!missing(excludeY)) inputString<-paste(inputString," --excludeY",sep="")
	if (!missing(extraYVar)) inputString<-paste(inputString," --extraYVar",sep="")
	if (!missing(varSelect)) inputString<-paste(inputString," --varSelect=",varSelect,sep="")
	if (!missing(entropy)) inputString<-paste(inputString," --entropy",sep="")

	.Call('profRegr', inputString, PACKAGE = 'DiPBaC')
}


readRunInfo<-function(directoryPath,fileStem='output'){

	# Read the number of sweeps, the length of burn in and filter from the log file
	runData<-readLines(file.path(directoryPath,paste(fileStem,'_log.txt',sep='')))

	# Number of sweeps
	nSweeps<-runData[grep('Number of sweeps',runData)]
	nSweeps<-substr(nSweeps,regexpr(':',nSweeps)+1,nchar(nSweeps))
	nSweeps<-gsub(' ','',nSweeps)
	nSweeps<-gsub('\t','',nSweeps)
	nSweeps<-as.integer(nSweeps)
 
	# Length of burn in
	nBurn<-runData[grep('Burn in sweeps',runData)]
	nBurn<-substr(nBurn,regexpr(':',nBurn)+1,nchar(nBurn))
	nBurn<-gsub(' ','',nBurn)
	nBurn<-gsub('\t','',nBurn)
	nBurn<-as.integer(nBurn)
	 
	# Output filter
	nFilter<-runData[grep('Output filter',runData)]
	nFilter<-substr(nFilter,regexpr(':',nFilter)+1,nchar(nFilter))
	nFilter<-gsub(' ','',nFilter)
	nFilter<-gsub('\t','',nFilter)
	nFilter<-as.integer(nFilter)
	 
	# Number of subjects
	nSubjects<-runData[grep('Number of subjects',runData)]
	nSubjects<-substr(nSubjects,regexpr(':',nSubjects)+1,nchar(nSubjects))
	nSubjects<-gsub(' ','',nSubjects)
	nSubjects<-gsub('\t','',nSubjects)	
	nSubjects<-as.integer(nSubjects)
 
	# Number of prediction subjects
	nPredictSubjects<-runData[grep('Number of prediction subjects',runData)]
	if(length(nPredictSubjects)==0){
		# This is for legacy runs
		nPredictSubjects<-0
	}else{
		nPredictSubjects<-substr(nPredictSubjects,regexpr(':',nPredictSubjects)+1,
		nchar(nPredictSubjects))
		nPredictSubjects<-gsub(' ','',nPredictSubjects)
		nPredictSubjects<-gsub('\t','',nPredictSubjects)
		nPredictSubjects<-as.integer(nPredictSubjects)
	}
	
	# Model for X
	xModel<-as.character(runData[grep('Model for X',runData)])
	xModel<-substr(xModel,regexpr(':',xModel)+1,nchar(xModel))
	xModel<-gsub(' ','',xModel)
	xModel<-gsub('\t','',xModel)

	if(length(grep('Include response: False',runData))>0){
		includeResponse<-F
		yModel<-NULL
		readNCategoriesY<-as.character(runData[grep('Model for Y',runData)])
	}else{
		includeResponse<-T
		# Model for Y
		yModel<-as.character(runData[grep('Model for Y',runData)])
		yModel<-substr(yModel,regexpr(':',yModel)+1,nchar(yModel))
		yModel<-gsub(' ','',yModel)
		yModel<-gsub('\t','',yModel)
		readNCategoriesY<-yModel
	}
	
	# Variable selection type
	if(length(grep("Variable selection",runData))==0){
		varSelect<-F
		varSelectType<-NULL
	}else{
		if(length(grep('Variable selection: None',runData))>0){
			varSelect<-F
			varSelectType<-NULL
		}else{
			varSelect<-T
			varSelectType<-as.character(runData[grep('Variable selection',runData)])
			varSelectType<-substr(varSelectType,regexpr(':',varSelectType)+1,nchar(varSelectType))
			varSelectType<-gsub(' ','',varSelectType)
			varSelectType<-gsub('\t','',varSelectType)
		}
	}
	# Input data file name
	inputFileName<-runData[grep('Data file path',runData)]
	inputFileName<-substr(inputFileName,regexpr(':',inputFileName)+1,nchar(inputFileName))
	inputFileName<-gsub(' ','',inputFileName)
	inputFileName<-gsub('\t','',inputFileName)

	# Get covariate and fixed effect information
	inputData<-readLines(inputFileName)
	inputData<-gsub("\t"," ",inputData)
	nCovariates<-as.integer(inputData[2])
	nFixedEffects<-as.integer(inputData[3+nCovariates])
	nExtraRows<-0   
	nCategoriesY<-1
	if ((includeResponse && yModel=='Categorical')||readNCategoriesY=="Categorical"){
		nCategoriesY<-as.integer(inputData[(4+nCovariates+nFixedEffects)])
		nExtraRows<-nExtraRows+1
	}
  	nCategories<-NULL
	if(xModel=='Discrete'){
		nCategories<-inputData[(4+nCovariates+nFixedEffects+nExtraRows)]
		nCategories<-as.integer(unlist(strsplit(nCategories," ")))
		nExtraRows<-nExtraRows+1
	}
	# Covariate names
	covNames<-inputData[3:(2+nCovariates)]

	xMat<-inputData[(4+nCovariates+nFixedEffects+nExtraRows):length(inputData)]
	xMat<-matrix(as.numeric(unlist(strsplit(xMat," "))),nrow=nSubjects,byrow=T)
	yMat<-NULL
	wMat<-NULL
	if(includeResponse){
		yMat<-matrix(xMat[,1],ncol=1)
		if(yModel=='Poisson'){
			offset<-xMat[,ncol(xMat)]
			yMat<-cbind(yMat,offset)
		}else if(yModel=='Binomial'){
			nTrials<-xMat[,ncol(xMat)]
			yMat<-cbind(yMat,nTrials)
		}
	
		if(nFixedEffects>0){
			wMat<-as.matrix(xMat[,(2+nCovariates):(1+nCovariates+nFixedEffects)])
		}
	}
	xMat<-matrix(xMat[,2:(nCovariates+1)],nrow=nSubjects)
	
	return(list("directoryPath"=directoryPath,
		"fileStem"=fileStem,
		"inputFileName"=inputFileName,
		"nSweeps"=nSweeps,
		"nBurn"=nBurn,
		"nFilter"=nFilter,
		"nSubjects"=nSubjects,
		"nPredictSubjects"=nPredictSubjects,
		"covNames"=covNames,
		"xModel"=xModel,
		"includeResponse"=includeResponse,
		"yModel"=yModel,
		"varSelect"=varSelect,
		"varSelectType"=varSelectType,
		"nCovariates"=nCovariates,
		"nFixedEffects"=nFixedEffects,
		"nCategoriesY"=nCategoriesY,
		"nCategories"=nCategories,
		"xMat"=xMat,"yMat"=yMat,"wMat"=wMat))

}


# Function to take the output from the C++ run and return an average dissimilarity
# matrix
calcDissimilarityMatrix<-function(runInfoObj){

   attach(runInfoObj)
   fileName <- file.path(directoryPath,paste(fileStem,'_z.txt',sep=''))
      
   # Call the C++ to compute the dissimilarity matrix
   disSimList<-.Call('calcDisSimMat',fileName,nSweeps,nBurn,nFilter,nSubjects,
                       nPredictSubjects, PACKAGE = 'DiPBaC')

   disSimMat<-disSimList$disSimMat
   lsOptSweep<-disSimList$lsOptSweep
   disSimMatPred<-NULL              
   if(nPredictSubjects>0){
      disSimMatPred<-disSimMat[(1+(nSubjects*(nSubjects-1)/2)):length(disSimMat)]
      disSimMat<-disSimMat[1:(nSubjects*(nSubjects-1)/2)]
   }   
   detach(runInfoObj)              
   disSimObj<-list('disSimRunInfoObj'=runInfoObj,'disSimMat'=disSimMat,
                     'disSimMatPred'=disSimMatPred,'lsOptSweep'=lsOptSweep)              
   return(disSimObj)
}

# Given a dissimilarity matrix (or list of dissimilarity matrices)
# run partitioning around medoids clustering
calcOptimalClustering<-function(disSimObj,maxNClusters=NULL,useLS=F){
   
   attach(disSimObj)
   attach(disSimRunInfoObj)

   if(useLS){
      zFileName <- file.path(directoryPath,paste(fileStem,'_z.txt',sep=''))
      zFile<-file(zFileName,open="r")
            
      optZ<-scan(zFile,what=integer(),skip=lsOptSweep-1,n=nSubjects+nPredictSubjects,quiet=T)
      optZFit<-optZ[1:nSubjects]
      if(nPredictSubjects>0){
         optZPredict<-optZ[(nSubjects+1):(nSubjects+nPredictSubjects)]
      }
      uZFit<-unique(optZFit)
      chosenNClusters<-length(unique(uZFit))
      clustVec<-match(optZFit,uZFit)
      clustSizes<-rep(0,chosenNClusters)
      for(c in 1:chosenNClusters){
         clustSizes[c]<-length(which(clustVec==c))
      }

      clusteringPred<-NULL
      if(nPredictSubjects>0){
         clusteringPred<-match(optZPredict,uZFit)
      }
      avgSilhouetteWidth<-NULL
      close(zFile)
      
	}else{
   
		if(is.null(maxNClusters)){
			# Determine the maximum number of clusters
			nMembersFileName<-file.path(directoryPath,paste(fileStem,'_nMembers.txt',sep=''))
			nMembersFile<-file(nMembersFileName,open="r")
			nClustersFileName<- file.path(directoryPath,paste(fileStem,'_nClusters.txt',sep=''))
			nClustersFile<-file(nClustersFileName,open="r")
   	
			# Restrict to sweeps after burn in
			firstLine<-2+nBurn/nFilter
			lastLine<-1+(nSweeps+nBurn)/nFilter
			maxNClusters<-0
   	
			for(sweep in firstLine:lastLine){
				if(sweep==firstLine){
					skipVal<-firstLine-1
				}else{
					skipVal<-0
				}
	   
				# Get the current number of members for each cluster
				nClusters<-scan(nClustersFile,what=integer(),skip=skipVal,n=1,quiet=T)	
				currNMembers<-scan(nMembersFile,what=integer(),skip=skipVal,n=nClusters+1,quiet=T)
				currNMembers<-currNMembers[1:nClusters]
				# Find the number of non-empty clusters
				nNotEmpty<-sum(currNMembers>0)
				if(nNotEmpty>maxNClusters){
					maxNClusters<-nNotEmpty
				}
			}   
			# Add on another 5 just to make sure bound is safe (but don't let it exceed no. of subjects)
			maxNClusters<-min(maxNClusters+5,nSubjects)

			close(nMembersFile)
			close(nClustersFile)
		}   
   
		# If the input was a list of dissimilarity matrices then take the average
		if(is.list(disSimMat)){
			for(i in 1:length(disSimMat)){
				if(i==1){
					tmpMat<-disSimMat[[i]]      
				}else{
					tmpMat<-tmpMat+disSimMat[[i]]
				}
			}
			tmpMat<-tmpMat/length(disSimMat)
			disSimMat<-tmpMat
		}
   
		# Loop over the possible number of clusters
		avgSilhouetteWidth<--1.0;
		cat(paste("Max no of possible clusters:",maxNClusters,"\n"))
		for(c in 2:maxNClusters){
			cat(paste("Trying",c,"clusters\n"))
			tmpObj<-pam(disSimMat,k=c,diss=T)
			# Check whether the silhouette width from this clustering improves previous best
			if(avgSilhouetteWidth<tmpObj$silinfo$avg.width){
				avgSilhouetteWidth<-tmpObj$silinfo$avg.width
				chosenNClusters<-c
				clustVec<-tmpObj$clustering
				clustSizes<-tmpObj$clusinfo[,1]
				# The id of the objects chosen as the medoids
				clustMedoids<-tmpObj$id.med
			}
		}
   
		# Work out the clustering of the prediction objects
		clusteringPred<-NULL
		if(nPredictSubjects>0){
			disSimMatPred<-matrix(disSimMatPred,nrow=nPredictSubjects,byrow=T)
			clusteringPred<-rep(0,nPredictSubjects)
			for(i in 1:nPredictSubjects){
				tmpVec<-disSimMatPred[i,clustMedoids]
				clusteringPred[i]<-which(tmpVec==min(tmpVec))
			}
		}
	}
	
	detach(disSimRunInfoObj)
	detach(disSimObj)
	return(list("clusObjRunInfoObj"=disSimObj$disSimRunInfoObj,
		"clusObjDisSimMat"=disSimObj$disSimMat,
		"nClusters"=chosenNClusters,
		"clusterSizes"=clustSizes,
		"clustering"=clustVec,
		"avgSilhouetteWidth"=avgSilhouetteWidth,
		"clusteringPred"=clusteringPred))
}

# Function to take the optimal clustering and computing the risk and probability
# profile
calcAvgRiskAndProfile<-function(clusObj,includeFixedEffects=F){

	attach(clusObj)
	attach(clusObjRunInfoObj)
	 
	# Construct the number of clusters file name
	nClustersFileName<-file.path(directoryPath,paste(fileStem,'_nClusters.txt',sep=''))
	nClustersFile<-file(nClustersFileName,open="r")
	
	# Construct the allocation file name
	zFileName <- file.path(directoryPath,paste(fileStem,'_z.txt',sep=''))
	zFile<-file(zFileName,open="r")
	
	if(xModel=="Discrete"){
		# Construct the allocation file name
		phiFileName <- file.path(directoryPath,paste(fileStem,'_phi.txt',sep=''))
		phiFile<-file(phiFileName,open="r")
		# Get the maximum number of categories
		maxNCategories<-max(nCategories)
		if(varSelect){
			nullPhiFileName <- file.path(directoryPath,paste(fileStem,'_nullPhi.txt',sep=''))
			if(varSelectType=="Continuous"){
				gammaFileName <- file.path(directoryPath,paste(fileStem,'_rho.txt',sep=''))
				gammaFile<-file(gammaFileName,open="r")
			}else{
				gammaFileName <- file.path(directoryPath,paste(fileStem,'_gamma.txt',sep=''))
				gammaFile<-file(gammaFileName,open="r")
			}
		}
	
	}else if(xModel=="Normal"){
		muFileName <- file.path(directoryPath,paste(fileStem,'_mu.txt',sep=''))
		muFile<-file(muFileName,open="r")
		SigmaFileName <- file.path(directoryPath,paste(fileStem,'_Sigma.txt',sep=''))
		SigmaFile<-file(SigmaFileName,open="r")
		if(varSelect){
			nullMuFileName <- file.path(directoryPath,paste(fileStem,'_nullMu.txt',sep=''))
			if(varSelectType=="Continuous"){
				gammaFileName <- file.path(directoryPath,paste(fileStem,'_rho.txt',sep=''))
				gammaFile<-file(gammaFileName,open="r")
			}else{
				gammaFileName <- file.path(directoryPath,paste(fileStem,'_gamma.txt',sep=''))
				gammaFile<-file(gammaFileName,open="r")
			}
		}
	}
	
	if(includeResponse){
		# Construct the theta file name
		thetaFileName <- file.path(directoryPath,paste(fileStem,'_theta.txt',sep=''))
		thetaFile<-file(thetaFileName,open="r")
		
		if(nFixedEffects>0){
			# Construct the fixed effect coefficient file name
			betaFileName <-file.path(directoryPath,paste(fileStem,'_beta.txt',sep=''))
			betaFile<-file(betaFileName,open="r")
		}
	}
	
	# Restrict to sweeps after burn in
	firstLine<-2+nBurn/nFilter
	lastLine<-1+(nSweeps+nBurn)/nFilter
	nSamples<-lastLine-firstLine+1
	
	# Make a list of the subjects in each of the optimal clusters
	optAlloc<-vector("list",nClusters)
	for(c in 1:nClusters){
		optAlloc[[c]]<-which(clustering==c)
	}
	
	if(includeResponse){
		# Initialise the object for storing the risks
		riskArray<-array(0,dim=c(nSamples,nClusters,nCategoriesY))
		thetaArray<-array(0,dim=c(nSamples,nClusters,nCategoriesY))
		if(nFixedEffects>0){
			betaArray<-array(0,dim=c(nSamples,nFixedEffects,nCategoriesY))
		}
	}else{
		riskArray<-NULL
	}
	
	# Initialise the object for storing the profiles
	if(xModel=='Discrete'){
		phiArray<-array(dim=c(nSamples,nClusters,nCovariates,maxNCategories))
		if(varSelect){
			phiStarArray<-array(dim=c(nSamples,nClusters,nCovariates,maxNCategories))
			tmpCurrNullPhi<-scan(nullPhiFileName,what=double(),quiet=T)
			tmpCurrNullPhi<-array(tmpCurrNullPhi,dim=c(maxNCategories,nCovariates))
			currNullPhi<-array(dim=c(1,maxNCategories,nCovariates))
			currNullPhi[1,,]<-tmpCurrNullPhi
		}else{
			phiStarArray<-NULL
		}
	}else if(xModel=='Normal'){
		muArray<-array(dim=c(nSamples,nClusters,nCovariates))
		if(varSelect){
			muStarArray<-array(dim=c(nSamples,nClusters,nCovariates))
			currNullMu<-scan(nullMuFileName,what=double(),quiet=T)
			currNullMu<-array(currNullMu,dim=c(nCovariates,1))
			currNullMu<-t(currNullMu)
		}else{
			muStarArray<-NULL
		}
		sigmaArray<-array(dim=c(nSamples,nClusters,nCovariates,nCovariates))
	}
	
	
	for(sweep in firstLine:lastLine){
		if(sweep==firstLine){
			skipVal<-firstLine-1
		}else{
			skipVal<-0
		}
	
		if(sweep-firstLine==0||(sweep-firstLine+1)%%1000==0){
			cat(paste("Processing sweep",sweep-firstLine+1,"of ",lastLine-firstLine+1,"\n"))
		}
	
		currMaxNClusters<-scan(nClustersFile,what=integer(),skip=skipVal,n=1,quiet=T)
	
		# Get the current allocation data for this sweep
		currZ<-scan(zFile,what=integer(),skip=skipVal,n=nSubjects+nPredictSubjects,quiet=T)
		currZ<-1+currZ
	
		if(includeResponse){
			# Get the risk data corresponding to this sweep
			if (yModel=="Categorical") {
				currThetaVector<-scan(thetaFile,what=double(),skip=skipVal,n=currMaxNClusters*(nCategoriesY-1),quiet=T)
				currTheta<-matrix(currThetaVector,ncol=(nCategoriesY-1),byrow=T)
				currTheta<-cbind(rep(0,dim(currTheta)[1]),currTheta)
			} else {
				currThetaVector<-scan(thetaFile,what=double(),skip=skipVal,n=currMaxNClusters*nCategoriesY,quiet=T)
				currTheta<-matrix(currThetaVector,ncol=nCategoriesY,byrow=T)
			}
			if(nFixedEffects>0){
				if (yModel=="Categorical") {
					currBetaVector<-scan(betaFile,what=double(),skip=skipVal,n=nFixedEffects*(nCategoriesY-1),quiet=T)
					currBeta<-matrix(currBetaVector,ncol=(nCategoriesY-1),byrow=T)
					currBeta<-cbind(rep(0,dim(currBeta)[1]),currBeta)
				} else {
					currBetaVector<-scan(betaFile,what=double(),skip=skipVal,n=nFixedEffects*nCategoriesY,quiet=T)
					currBeta<-matrix(currBetaVector,ncol=nCategoriesY,byrow=T)
				}
				betaArray[sweep-firstLine+1,,]<-currBeta
			}
			# Calculate the average risk (over subjects) for each cluster
			for(c in 1:nClusters){
				currLambdaVector<-currTheta[currZ[optAlloc[[c]]],]
				currLambda<-matrix(currLambdaVector,ncol=nCategoriesY)
				if(includeFixedEffects&&nFixedEffects>0){
					if (yModel=="Categorical"){
						for (i in 1:length(optAlloc[[c]])){
							for (k in 1:nCategoriesY) currLambda[i,k]<-currLambda[i,k]+
								t(as.matrix(wMat[optAlloc[[c]][i],]))%*%
								currBeta[,yMat[optAlloc[[c]][i]]+1]
						}	
					} else {
						currLambda<-currLambda+as.matrix(wMat[optAlloc[[c]],])%*%currBeta
					}
				}
				if(yModel=="Poisson"){
					currRisk<-exp(currLambda)
				}else if(yModel=="Bernoulli"||yModel=="Binomial"){
					currRisk<-1.0/(1.0+exp(-currLambda))
				}else if(yModel=="Normal"){
					currRisk<-currLambda
				}else if(yModel=="Categorical"){
					currRisk<-matrix(0,ncol=length(optAlloc[[c]]),nrow=nCategoriesY)
					currRisk<-exp(currLambda)/rowSums(exp(currLambda))
				}
				riskArray[sweep-firstLine+1,c,]<-apply(currRisk,2,mean)
				thetaArray[sweep-firstLine+1,c,]<-apply(as.matrix(currTheta[currZ[optAlloc[[c]]],],ncol=nCategoriesY),2,mean)
			}
		}
	
		# Calculate the average profile (over subjects) for each cluster
		if(xModel=='Discrete'){
			currPhi<-scan(phiFile,what=double(),
				skip=skipVal,n=currMaxNClusters*maxNCategories*nCovariates,quiet=T)
				# This is slightly convoluted, because of the way that R reads in by column
				# I switched the order of categories and covariates in column below, and then
				# take the transpose to correct in the loop
				currPhi<-array(currPhi,dim=c(currMaxNClusters,maxNCategories,nCovariates))
				if(varSelect){
					# We increase dimensions of currGamma using duplicates, to
					# enable easier calculation of phiStar
					if(varSelectType=='BinaryCluster'){
						tmpCurrGamma<-scan(gammaFile,what=double(),
							skip=skipVal,n=nCovariates*currMaxNClusters,quiet=T)
						tmpCurrGamma<-array(tmpCurrGamma,dim=c(currMaxNClusters,nCovariates))
					}else{
						tmpCurrGamma<-scan(gammaFile,what=double(),skip=skipVal,n=nCovariates,quiet=T)
						tmpCurrGamma<-array(tmpCurrGamma,dim=c(nCovariates,currMaxNClusters))
						tmpCurrGamma<-t(tmpCurrGamma)
					}
					currGamma<-array(dim=c(currMaxNClusters,maxNCategories,nCovariates))
				for(p in 1:maxNCategories){
					currGamma[,p,]<-tmpCurrGamma
				}
			}	
			for(c in 1:nClusters){
				phiArray[sweep-firstLine+1,c,,]<-t(apply(array(currPhi[currZ[optAlloc[[c]]],,],
					dim=c(length(optAlloc[[c]]),
					dim(currPhi)[2],dim(currPhi)[3])),2:3,mean))
				if(varSelect){
					phiStarArray[sweep-firstLine+1,c,,]<-t(apply(array(currGamma[currZ[optAlloc[[c]]],,],
						dim=c(length(optAlloc[[c]]),dim(currGamma)[2],
						dim(currGamma)[3]))*array(currPhi[currZ[optAlloc[[c]]],,],
						dim=c(length(optAlloc[[c]]),dim(currPhi)[2],dim(currPhi)[3]))+
						(1-array(currGamma[currZ[optAlloc[[c]]],,],
						dim=c(length(optAlloc[[c]]),dim(currGamma)[2],dim(currGamma)[3])))*
						array(currNullPhi[rep(1,length(optAlloc[[c]])),,],
						dim=c(length(optAlloc[[c]]),dim(currNullPhi)[2],
						dim(currNullPhi)[3])),2:3,mean))
				}
			}
		}else if(xModel=='Normal'){
			# mu stored like phi
			currMu<-scan(muFile,what=double(),skip=skipVal,n=currMaxNClusters*nCovariates,quiet=T)
			currMu<-array(currMu,dim=c(currMaxNClusters,nCovariates))
			if(varSelect){
				# We increase dimensions of nullPhi and currGamma using duplicates, to
				# enable easier calculation of phiStar
				# We increase dimensions of currGamma using duplicates, to
				# enable easier calculation of phiStar
				if(varSelectType=='BinaryCluster'){
					tmpCurrGamma<-scan(gammaFile,what=double(),skip=skipVal,n=currMaxNClusters*nCovariates,quiet=T)
					currGamma<-array(tmpCurrGamma,dim=c(currMaxNClusters,nCovariates))
				}else{
					tmpCurrGamma<-scan(gammaFile,what=double(),skip=skipVal,n=nCovariates,quiet=T)
					tmpCurrGamma<-array(tmpCurrGamma,dim=c(nCovariates,currMaxNClusters))
					currGamma<-t(tmpCurrGamma)
				}
	
			}
			for(c in 1:nClusters){
				muArray[sweep-firstLine+1,c,]<-apply(matrix(currMu[currZ[optAlloc[[c]]],],ncol=nCovariates),2,mean)
				if(varSelect){
					muStarArray[sweep-firstLine+1,c,]<-apply(matrix(currGamma[currZ[optAlloc[[c]]],],
						ncol=nCovariates)*matrix(currMu[currZ[optAlloc[[c]]],],ncol=nCovariates)+
						matrix(1-currGamma[currZ[optAlloc[[c]]],],ncol=nCovariates)*
						matrix(currNullMu[rep(1,length(optAlloc[[c]])),],ncol=nCovariates),2,mean)
				}
			}
	
			currSigma<-scan(SigmaFile,what=double(),skip=skipVal,n=currMaxNClusters*nCovariates*nCovariates,quiet=T)
			currSigma<-array(currSigma,dim=c(currMaxNClusters,nCovariates,nCovariates))
			for(c in 1:nClusters){
				sigmaArray[sweep-firstLine+1,c,,]<-apply(array(currSigma[currZ[optAlloc[[c]]],,],
					dim=c(length(optAlloc[[c]]),dim(currSigma)[2],dim(currSigma)[3])),2:3,mean)
			}
		}
	}
	
	
	# Calculate the empiricals
	empiricals<-rep(0,nClusters)
	if(!is.null(yModel)){
		for(c in 1:nClusters){
			if(yModel=='Bernoulli'||yModel=='Normal'){
				empiricals[c]<-mean(yMat[optAlloc[[c]],1])
			}else if(yModel=='Binomial'){
				empiricals[c]<-mean(yMat[optAlloc[[c]],1]/yMat[optAlloc[[c]],2])
			}else if(yModel=='Poisson'){
				empiricals[c]<-mean(yMat[optAlloc[[c]],1]/yMat[optAlloc[[c]],2])
			}
		}
	}

	if(xModel=='Discrete'){
		out<-list('riskProfClusObj'=clusObj,'risk'=riskArray,'profile'=phiArray,'profileStar'=phiStarArray,'empiricals'=empiricals)
	}else if(xModel=='Normal'){
		out<-list('riskProfClusObj'=clusObj,'risk'=riskArray,
			'profile'=muArray,'profileStar'=muStarArray,
			'profileStdDev'=sigmaArray,'empiricals'=empiricals)
	}
	
	close(zFile)
	close(nClustersFile)
	if(xModel=="Discrete"){
		close(phiFile)
		if(varSelect){
			close(gammaFile)
		}
	}else if(xModel=="Normal"){
		close(muFile)
		close(SigmaFile)
		if(varSelect){
			close(gammaFile)
		}
	}
	
	if(includeResponse){
		close(thetaFile)
		if(nFixedEffects>0){
			close(betaFile)
		}
	}
	
	detach(clusObjRunInfoObj)
	detach(clusObj)
		
	return(out)
}
	
		

# Plot output values
plotRiskProfile<-function(riskProfObj,outFile,showRelativeRisk=F,orderBy=NULL,whichClusters=NULL,whichCovariates=NULL,useProfileStar=F){

	attach(riskProfObj)
	attach(riskProfClusObj)
	attach(clusObjRunInfoObj)

	if(includeResponse){
		if(yModel=="Normal"){
			showRelativeRisk<-F
		}
	}

	if(useProfileStar){
		profile<-profileStar
	}
	if(!is.null(whichCovariates)){
		if(xModel=='Discrete'){
			profile<-profile[,,whichCovariates,]
			nCategories<-nCategories[whichCovariates]
		}else if(xModel=='Normal'){
			profile<-profile[,,whichCovariates]
			profileStdDev<-profileStdDev[,,whichCovariates,whichCovariates]
		}
		covNames<-covNames[whichCovariates]
		nCovariates<-length(whichCovariates)
	}

	png(outFile,width=1200,height=800)
	
	if(!is.null(orderBy)){
		if(!includeResponse){
			if(orderBy!='Empirical'&&orderBy!='ClusterSize'&&!orderBy%in%covNames){
				orderBy<-NULL
			}
		}else{
			if(orderBy!='Risk'&&orderBy!='Empirical'&&orderBy!='ClusterSize'&&!orderBy%in%covNames){
				orderBy<-NULL
			}
		}
	}
	
	# Set up the layout for the plot
	plotLayout<-grid.layout(ncol = nCovariates+2, nrow = 6)
	grid.newpage()
	pushViewport(viewport(layout = plotLayout))
	
	orderProvided<-F
	if(is.numeric(orderBy)){
		if(length(orderBy)==nClusters){
			orderProvided<-T
			meanSortIndex<-orderBy
		}else{
			cat("Order vector provided not of same length as number of clusters. Reverting to default ordering.\n")
			orderBy<-NULL
		}
	}

	if(!orderProvided){
		if(!is.null(risk)){
			if(is.null(orderBy)){
				# Default is to order by posterior theta risk
				# Compute the means
				orderStat<-apply(risk,2,median)
			}else{
				if(orderBy=='Risk'){
					orderStat<-apply(risk,2,median)
				}else if(orderBy=='Empirical'){
					orderStat<-empiricals
				}else if(orderBy=='ClusterSize'){
					orderStat<-clusterSizes
				}else{
					whichCov<-match(orderBy,covNames)
					if(xModel=='Normal'){
						orderStat<-apply(profile[,,whichCov],2,median)
					}else{
						# This assumes that there is some order to the categories
						# and then uses an expected value
						tmpMat<-profile[,,whichCov,1]
						if(nCategories[whichCov]>1){
							for(k in 2:nCategories[whichCov]){
								tmpMat<-tmpMat+k*profile[,,whichCov,k]
							}
						}
						orderStat<-apply(tmpMat,2,median)
					}
				}
			}
		}else{
			if(is.null(orderBy)){
				# Default is to order by empirical risk
				orderStat<-empiricals
			}else{
				if(orderBy=='Empirical'){
					orderStat<-empiricals
				}else if(orderBy=='ClusterSize'){
					orderStat<-clusterSizes
				}else{
					whichCov<-match(orderBy,covNames)
					if(xModel=='Normal'){
						orderStat<-apply(profile[,,whichCov],2,median)
					}else{
						# This assumes that there is some order to the categories
						# and then uses an expected value
						tmpMat<-profile[,,whichCov,1]
						if(nCategories[whichCov]>1){
							for(k in 2:nCategories[whicCov]){
								tmpMat<-tmpMat+k*profile[,,whichCov,k]
							}
						}
						orderStat<-apply(tmpMat,2,median)
					}
				}
			}	
		}
		# Sort into ascending mean size
		meanSortIndex<-order(orderStat,decreasing=F)
	}
	if(includeResponse){
		# Reorder the risk matrix
		riskDim<-dim(risk)
		risk<-array(risk[,meanSortIndex,],dim=riskDim)
		if(showRelativeRisk){
			for(c in nClusters:1){
				risk[,c,]<-risk[,c,]/risk[,1,]
			}
		}
	}

	# Reorder the cluster sizes
	clusterSizes<-clusterSizes[meanSortIndex]
	# Reorder the empiricals
	empiricals<-empiricals[meanSortIndex]
	meanEmpirical<-sum(empiricals*clusterSizes)/sum(clusterSizes)

	if(includeResponse){
		# Recompute the means and now also credible intervals
		riskMeans<-apply(risk,2,mean,trim=0.005)
		riskMean<-sum(riskMeans*clusterSizes)/sum(clusterSizes)
		riskLower<-apply(risk,2,quantile,0.05)
		riskUpper<-apply(risk,2,quantile,0.95)
		# The next line is to avoid outliers spoiling plot scales
		plotMax<-2*max(riskUpper)-riskMean
		
		# Get the plot colors
		riskColor<-ifelse(riskLower>rep(riskMean,nClusters),"high",
		ifelse(riskUpper<rep(riskMean,nClusters),"low","avg"))
		if (yModel=="Categorical"){		
			riskDF<-data.frame("risk"=c(),"category"=c(),"cluster"=c(),"meanRisk"=c(),
				"lowerRisk"=c(),"upperRisk"=c(),"fillColor"=c())
		} else {
			riskDF<-data.frame("risk"=c(),"cluster"=c(),"meanRisk"=c(),
				"lowerRisk"=c(),"upperRisk"=c(),"fillColor"=c())
		}
	}else{
		riskColor<-ifelse(empiricals>rep(meanEmpirical,nClusters),"high",
		ifelse(empiricals<rep(meanEmpirical,nClusters),"low","avg"))
	}

	if(is.null(whichClusters)){
		whichClusters<-1:nClusters
	}
	nClusters<-length(whichClusters)
	
	empiricalDF<-data.frame("empiricals"=c(),"meanEmpirical"=c(),"cluster"=c(),"fillColor"=c())
	sizeDF<-data.frame("clusterSize"=c(),"cluster"=c(),"fillColor"=c())
	
	# Restructure the data for plotting
	for(c in whichClusters){
		if(includeResponse){
			if (yModel=="Categorical"){
				plotRisk<-risk[,c,]
				nPoints<-dim(plotRisk)[1]
				for (k in 1:nCategoriesY){				
					riskDF<-rbind(riskDF,data.frame("risk"=plotRisk[,k],
						"category"=rep(k,nPoints),
						"cluster"=rep(c,nPoints),
						"meanRisk"=rep(riskMean,nPoints),
						"lowerRisk"=rep(riskLower[c],nPoints),
						"upperRisk"=rep(riskUpper[c],nPoints),
						"fillColor"=rep(riskColor[c],nPoints)))
				}
			} else {
				plotRisk<-risk[,c,]
				plotRisk<-plotRisk[plotRisk<plotMax]
				nPoints<-length(plotRisk)
				riskDF<-rbind(riskDF,data.frame("risk"=plotRisk,"cluster"=rep(c,nPoints),
					"meanRisk"=rep(riskMean,nPoints),
					"lowerRisk"=rep(riskLower[c],nPoints),
					"upperRisk"=rep(riskUpper[c],nPoints),
					"fillColor"=rep(riskColor[c],nPoints)))
			}
		}
		empiricalDF<-rbind(empiricalDF,
			data.frame("empiricals"=empiricals[c],
			"meanEmpirical"=meanEmpirical,"cluster"=c,"fillColor"=riskColor[c]))
		sizeDF<-rbind(sizeDF,
			data.frame("clusterSize"=clusterSizes[c],"cluster"=c,"fillColor"=riskColor[c]))
	}
	
	if(includeResponse){
		if(yModel=='Categorical'){
			riskDF<-data.frame("prob"=c(),"cluster"=c(),"category"=c(),"meanProb"=c(),
				"lowerProb"=c(),"upperProb"=c(),"fillColor"=c())
			for(k in 1:nCategoriesY){
		
				probMat<-risk[,meanSortIndex,k]
				nPoints<-nrow(probMat)
				probMeans<-apply(probMat,2,mean)
				probMean<-sum(probMeans*clusterSizes)/sum(clusterSizes)
				probLower<-apply(probMat,2,quantile,0.05)
				probUpper<-apply(probMat,2,quantile,0.95)
		
				# Get the plot colors
				probColor<-ifelse(probLower>rep(probMean,nClusters),"high",
				ifelse(probUpper<rep(probMean,nClusters),"low","avg"))

				for(c in whichClusters){
					riskDF<-rbind(riskDF,data.frame("prob"=probMat[,c],"cluster"=rep(c,nPoints),
						"category"=rep(k-1,nPoints),
						"meanProb"=rep(probMean,nPoints),
						"lowerProb"=rep(probLower[c],nPoints),
						"upperProb"=rep(probUpper[c],nPoints),
						"fillColor"=rep(probColor[c],nPoints)))
						 rownames(riskDF)<-seq(1,nrow(riskDF),1)
				
				}
			}

			plotObj<-ggplot(riskDF)
			plotObj<-plotObj+facet_wrap(~category,ncol=1,as.table=F,scales="free_y")
			plotObj<-plotObj+geom_hline(aes(x=as.factor(cluster),y=prob,yintercept=meanProb))
			plotObj<-plotObj+geom_boxplot(aes(x=as.factor(cluster),y=prob,fill=as.factor(fillColor)),outlier.size=0.5)
			plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=lowerProb,colour=as.factor(fillColor)),size=1.5)
			plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=upperProb,colour=as.factor(fillColor)),size=1.5)
			plotObj<-plotObj+
				scale_fill_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
				scale_colour_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
				theme(legend.position="none")+labs(x="Cluster")+theme(axis.title.x=element_text(size=10))
			plotObj<-plotObj+labs(y="Probability")+theme(axis.title.y=element_text(size=10,angle=90))
			plotObj<-plotObj+labs(title=ifelse(showRelativeRisk,'Relative Risk','Risk'),plot.title=element_text(size=10))
			# Margin order is (top,right,bottom,left)
			plotObj<-plotObj+theme(plot.margin=unit(c(0.5,0.15,0.5,0.15),'lines'))+
				theme(plot.margin=unit(c(0,0,0,0),'lines'))
			print(plotObj,vp=viewport(layout.pos.row=1:6,layout.pos.col=2))
	   }else{
			rownames(riskDF)<-seq(1,nrow(riskDF),1)


			# Create the risk plot
			plotObj<-ggplot(riskDF)
			plotObj<-plotObj+geom_hline(aes(x=as.factor(cluster),y=risk,yintercept=meanRisk))
			plotObj<-plotObj+geom_boxplot(aes(x=as.factor(cluster),y=risk,fill=as.factor(fillColor)),outlier.size=0.5)
			plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=lowerRisk,colour=as.factor(fillColor)),size=1.5)
			plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=upperRisk,colour=as.factor(fillColor)),size=1.5)
			plotObj<-plotObj+scale_fill_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
				scale_colour_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
				theme(legend.position="none")+
				labs(x="Cluster",y=ifelse(showRelativeRisk,'RR',
				ifelse(yModel=="Categorical"||yModel=="Bernoulli"||yModel=="Binomial","Probability","E[Y]")))
			plotObj<-plotObj+theme(axis.title.y=element_text(size=10,angle=90),axis.title.x=element_text(size=10))
			plotObj<-plotObj+labs(title=ifelse(showRelativeRisk,'Relative Risk','Risk'),plot.title=element_text(size=10))
			# Margin order is (top,right,bottom,left)
			plotObj<-plotObj+theme(plot.margin=unit(c(0,0,0,0),'lines'))+theme(plot.margin=unit(c(0.5,0.15,0.5,0.15),'lines'))
			print(plotObj,vp=viewport(layout.pos.row=1:6,layout.pos.col=2))
		}	
	}
	
	# Create a bar chart of cluster empiricals
	if(!is.null(yModel)){
		plotObj<-ggplot(empiricalDF)
		plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=empiricals,colour=as.factor(fillColor)),size=3)
		plotObj<-plotObj+geom_hline(aes(x=as.factor(cluster),y=empiricals,yintercept=meanEmpirical))
		plotObj<-plotObj+scale_colour_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
			theme(legend.position="none")
		plotObj<-plotObj+labs(title='Empirical Data',plot.title=element_text(size=10))
		plotObj<-plotObj+theme(axis.title.x=element_text(size=10),axis.title.y=element_text(size=10,angle=90))
		plotObj<-plotObj+
			labs(y=ifelse(yModel=="Bernoulli","Proportion of cases",
			ifelse(yModel=="Binomial","Avg Proportion of occurrence",
			ifelse(yModel=="Poisson","Avg Count",
			ifelse(yModel=="Categorical","Avg Proportion of occurrence","Avg Y")))),x="Cluster")
		plotObj<-plotObj+theme(plot.margin=unit(c(0,0,0,0),'lines'))+theme(plot.margin=unit(c(0.15,0.5,0.5,1),'lines'))
		print(plotObj,vp=viewport(layout.pos.row=1:3,layout.pos.col=1))
	}
	# Create a bar chart of cluster sizes
	plotObj<-ggplot(sizeDF)
	plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=clusterSize,colour=as.factor(fillColor)),size=3)
	plotObj<-plotObj+scale_colour_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+theme(legend.position="none")
	plotObj<-plotObj+labs(title="Size",plot.title=element_text(size=10))
	plotObj<-plotObj+theme(axis.title.x=element_text(size=10),axis.title.y=element_text(size=10,angle=90))
	plotObj<-plotObj+labs(y="No. of Subjects",x="Cluster")
	plotObj<-plotObj+theme(plot.margin=unit(c(0,0,0,0),'lines'))+theme(plot.margin=unit(c(0.15,0.5,0.5,1),'lines'))
	print(plotObj,vp=viewport(layout.pos.row=4:6,layout.pos.col=1))

	# Loop over the covariates
	for(j in 1:nCovariates){
		if(xModel=='Discrete'){
			profileDF<-data.frame("prob"=c(),"cluster"=c(),"category"=c(),"meanProb"=c(),
				"lowerProb"=c(),"upperProb"=c(),"fillColor"=c())
			for(k in 1:nCategories[j]){
				probMat<-profile[,meanSortIndex,j,k]
				nPoints<-nrow(probMat)
				probMeans<-apply(probMat,2,mean)
				probMean<-sum(probMeans*clusterSizes)/sum(clusterSizes)
				probLower<-apply(probMat,2,quantile,0.05)
				probUpper<-apply(probMat,2,quantile,0.95)
		
				# Get the plot colors
				probColor<-ifelse(probLower>rep(probMean,nClusters),"high",
				ifelse(probUpper<rep(probMean,nClusters),"low","avg"))
			
	
				for(c in whichClusters){
					profileDF<-rbind(profileDF,data.frame("prob"=probMat[,c],"cluster"=rep(c,nPoints),
						"category"=rep(k-1,nPoints),
						"meanProb"=rep(probMean,nPoints),
						"lowerProb"=rep(probLower[c],nPoints),
						"upperProb"=rep(probUpper[c],nPoints),
						"fillColor"=rep(probColor[c],nPoints)))
						 rownames(profileDF)<-seq(1,nrow(profileDF),1)
				
				}
			}
		
			plotObj<-ggplot(profileDF)
			plotObj<-plotObj+facet_wrap(~category,ncol=1,as.table=F,scales="free_y")
			plotObj<-plotObj+geom_hline(aes(x=as.factor(cluster),y=prob,yintercept=meanProb))
			plotObj<-plotObj+geom_boxplot(aes(x=as.factor(cluster),y=prob,fill=as.factor(fillColor)),outlier.size=0.5)
			plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=lowerProb,colour=as.factor(fillColor)),size=1.5)
			plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=upperProb,colour=as.factor(fillColor)),size=1.5)
			plotObj<-plotObj+
				scale_fill_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
				scale_colour_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
				theme(legend.position="none")+labs(x="Cluster")+theme(axis.title.x=element_text(size=10))
			if(j==1){
				plotObj<-plotObj+labs(y="Probability")+theme(axis.title.y=element_text(size=10,angle=90))
			}else{
				plotObj<-plotObj+theme(axis.title.y=element_blank())
			}
			plotObj<-plotObj+labs(title=covNames[j],plot.title=element_text(size=10))
			plotObj<-plotObj+theme(plot.margin=unit(c(0.5,ifelse(j==nCovariates,1,0),0.5,ifelse(j==1,0.5,0)),'lines'))+
				theme(plot.margin=unit(c(0,0,0,0),'lines'))
		
			print(plotObj,vp=viewport(layout.pos.row=1:6,layout.pos.col=j+2))
		}else if(xModel=='Normal'){
			# Plot the means
			profileDF<-data.frame("mu"=c(),"cluster"=c(),"muMean"=c(),
				"lowerMu"=c(),"upperMu"=c(),"fillColor"=c())
			muMat<-profile[,meanSortIndex,j]
			muMeans<-apply(muMat,2,mean)
			muMean<-sum(muMeans*clusterSizes)/sum(clusterSizes)
			muLower<-apply(muMat,2,quantile,0.05)
			muUpper<-apply(muMat,2,quantile,0.95)
			# The next line is to avoid outliers spoiling plot scales
			plotMax<-2*max(muUpper)-muMean
			plotMin<-2*min(muLower)-muMean
			
			# Get the plot colors
			muColor<-ifelse(muLower>rep(muMean,nClusters),"high",
			ifelse(muUpper<rep(muMean,nClusters),"low","avg"))
			for(c in whichClusters){
				plotMu<-muMat[,c]
				plotMu<-plotMu[plotMu<plotMax&plotMu>plotMin]
				nPoints<-length(plotMu)
				profileDF<-rbind(profileDF,data.frame("mu"=plotMu,"cluster"=rep(c,nPoints),
					"meanMu"=rep(muMean,nPoints),
					"lowerMu"=rep(muLower[c],nPoints),
					"upperMu"=rep(muUpper[c],nPoints),
					"fillColor"=rep(muColor[c],nPoints)))
			}
			rownames(profileDF)<-seq(1,nrow(profileDF),1)
			
			plotObj<-ggplot(profileDF)
			plotObj<-plotObj+geom_hline(aes(x=as.factor(cluster),y=mu,yintercept=meanMu))
			plotObj<-plotObj+geom_boxplot(aes(x=as.factor(cluster),y=mu,fill=as.factor(fillColor)),outlier.size=0.5)
			plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=lowerMu,colour=as.factor(fillColor)),size=1.5)
			plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=upperMu,colour=as.factor(fillColor)),size=1.5)
			plotObj<-plotObj+
				scale_fill_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
				scale_colour_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
				theme(legend.position="none")+labs(x="Cluster")+theme(axis.title.x=element_text(size=10))
			if(j==1){
				plotObj<-plotObj+labs(y="Mean")+theme(axis.title.y=element_text(size=10,angle=90))
			}else{
				plotObj<-plotObj+theme(axis.title.y=element_blank())
			}
				plotObj<-plotObj+labs(title=covNames[j],plot.title=element_text(size=10))
				plotObj<-plotObj+
					theme(plot.margin=unit(c(0.5,ifelse(j==nCovariates,1,0),0.5,ifelse(j==1,0.5,0)),'lines'))+
					theme(plot.margin=unit(c(0,0,0,0),'lines'))
		
			print(plotObj,vp=viewport(layout.pos.row=1:3,layout.pos.col=j+2))
		
			# Plot the variances
			profileDF<-data.frame("sigma"=c(),"cluster"=c(),"sigmaMean"=c(),
				"lowerSigma"=c(),"upperSigma"=c(),"fillColor"=c())
			sigmaMat<-profileStdDev[,meanSortIndex,j,j]
			sigmaMeans<-apply(sigmaMat,2,mean)
			sigmaMean<-sum(sigmaMeans*clusterSizes)/sum(clusterSizes)
			sigmaLower<-apply(sigmaMat,2,quantile,0.05)
			sigmaUpper<-apply(sigmaMat,2,quantile,0.95)
			# The next line is to avoid outliers spoiling plot scales
			plotMax<-2*max(sigmaUpper)-sigmaMean
	
			# Get the plot colors
			sigmaColor<-ifelse(sigmaLower>rep(sigmaMean,nClusters),"high",
			ifelse(sigmaUpper<rep(sigmaMean,nClusters),"low","avg"))
			for(c in whichClusters){
				plotSigma<-sigmaMat[,c]
				plotSigma<-plotSigma[plotSigma<plotMax]
				nPoints<-length(plotSigma)
				profileDF<-rbind(profileDF,data.frame("sigma"=plotSigma,"cluster"=rep(c,nPoints),
					"meanSigma"=rep(sigmaMean,nPoints),
					"lowerSigma"=rep(sigmaLower[c],nPoints),
					"upperSigma"=rep(sigmaUpper[c],nPoints),
					"fillColor"=rep(sigmaColor[c],nPoints)))
			}
			rownames(profileDF)<-seq(1,nrow(profileDF),1)
			plotObj<-ggplot(profileDF)
			plotObj<-plotObj+geom_hline(aes(x=as.factor(cluster),y=sigma,yintercept=meanSigma))
			plotObj<-plotObj+geom_boxplot(aes(x=as.factor(cluster),y=sigma,fill=as.factor(fillColor)),outlier.size=0.5)
			plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=lowerSigma,colour=as.factor(fillColor)),size=1.5)
			plotObj<-plotObj+geom_point(aes(x=as.factor(cluster),y=upperSigma,colour=as.factor(fillColor)),size=1.5)
			plotObj<-plotObj+
				scale_fill_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
				scale_colour_manual(values = c(high ="#CC0033",low ="#0066CC", avg ="#33CC66"))+
				theme(legend.position="none")+labs(x="Cluster")+theme(axis.title.x=element_text(size=10))
			if(j==1){
				plotObj<-plotObj+labs(y="Std Dev")+theme(axis.title.y=element_text(size=10,angle=90))
			}else{
				plotObj<-plotObj+theme(axis.title.y=element_blank())
			}
			plotObj<-plotObj+
				theme(plot.margin=unit(c(0.5,ifelse(j==nCovariates,1,0),0.5,ifelse(j==1,0.5,0)),'lines'))+
				theme(plot.margin=unit(c(0,0,0,0),'lines'))
		
			print(plotObj,vp=viewport(layout.pos.row=4:6,layout.pos.col=j+2))
	
	
		}
	}
	dev.off()
	detach(clusObjRunInfoObj)
	detach(riskProfClusObj)
	detach(riskProfObj)
	
	return(meanSortIndex)
	
}
	
plotClustering<-function(clusObj,outFile,clusterPlotOrder=NULL,whichCovariates=NULL){
	
	attach(clusObj)
	attach(clusObjRunInfoObj)
	
	png(outFile,width=1200,height=800)
	
	if(is.null(clusterPlotOrder)){
		clusterPlotOrder<-1:nClusters
	}
	# Read in the raw X data
	rawXFileName<-gsub('.txt','_RawX.txt',inputFileName)
	rawXData<-readLines(rawXFileName)
	rawXData<-rawXData[(nCovariates+3):length(rawXData)]
	rawXData<-matrix(as.numeric(unlist(strsplit(rawXData," "))),ncol=nCovariates,byrow=T)
	if('CessationTime'%in%covNames){
		relInd<-match('CessationTime',covNames)
		meanCat1<-mean(rawXData[xMat[,relInd]==1,relInd])
		rawXData[xMat[,relInd]==0,relInd]<-meanCat1+10
	}
	rawXData[rawXData==-999]<-NA
	includeVec<-rep(T,nSubjects)
	for(i in 1:nSubjects){
		if(any(is.na(rawXData[i,]))){
			includeVec[i]<-F
		}
	}
	rawXData<-rawXData[includeVec,]
	if(!is.null(whichCovariates)){
		rawXData<-rawXData[,whichCovariates]
	}
	clustering<-clustering[includeVec]
	d<-dist(rawXData)
	
	principalComp<-cmdscale(d,k=2)
	plot(principalComp[,1],principalComp[,2],pch=match(clustering,clusterPlotOrder),
		col=match(clustering,clusterPlotOrder),xlab="Principal Component 1",ylab="Principal Component 2")
	title(main="2D visualization of how subjects cluster")
	legend("topleft",legend=1:nClusters,pch=1:nClusters,col=1:nClusters,bg="white")
	box()
	dev.off()
	
	detach(clusObjRunInfoObj)
	detach(clusObj)
	
	return(list("rawXDist"=d,"rawXprincipalComponents"=principalComp,"rawXIncludeVec"=includeVec))
}
	
	
# Calculate predictions, and if possible assess predictive performance
calcPredictions<-function(riskProfObj,predictResponseFileName=NULL,doRaoBlackwell=F,fullSweepPredictions=F,fullSweepLogOR=F){
	
	attach(riskProfObj)
	attach(riskProfClusObj)
	attach(clusObjRunInfoObj)
	
	if(yModel=="Poisson"||yModel=="Normal"){
		fullSweepLogOR=F
		cat("Log odds ratio does not make sense for Poisson or Normal response\n")
	}
	
	firstLine<-2+nBurn/nFilter
	lastLine<-1+(nSweeps+nBurn)/nFilter
	nSamples<-lastLine-firstLine+1
	
	# First of all we see if there a covariate file has been supplied
	if(nPredictSubjects==0){
		detach(clusObjRunInfoObj)
		detach(riskProfClusObj)
		detach(riskProfObj)
		stop("No prediction subjects processed by C++\n")
	}
	
	# Get the response and fixed effects data if available
	responseProvided<-F
	extraInfoProvided<-F
	fixedEffectsProvided<-F
	if(!is.null(predictResponseFileName)){
		predictResponseData<-scan(predictResponseFileName,quiet=T)
		predictResponseMat<-matrix(predictResponseData[2:length(predictResponseData)],
		nrow=nPredictSubjects,byrow=T)
		predictYMat<-matrix(predictResponseMat[,1],ncol=1)
		if(yModel=="Poisson"||yModel=="Binomial"){
			predictYMat<-cbind(predictYMat,predictResponseMat[,ncol(predictResponseMat)])
			if(all(predictYMat[,2]>-999)){
				extraInfoProvided<-T
			}
		}
		if(all(predictYMat[,1]>-999)){
			responseProvided<-T
		}
		if(nFixedEffects>0){
			fixedEffectsProvided<-T
			predictWMat<-matrix(predictResponseMat[,2:(nFixedEffects+1)],nrow=nPredictSubjects)
			# Set missing values to their average or reference value
			predictWMat[predictWMat==-999]<-0
		}
	}

	if(fixedEffectsProvided){
		betaFileName <-file.path(directoryPath,paste(fileStem,'_beta.txt',sep=''))
		betaFile<-file(betaFileName,open="r")
		betaArray<-array(0,dim=c(nSamples,nFixedEffects,nCategoriesY))
		for(sweep in firstLine:lastLine){
			if(sweep==firstLine){
				skipVal<-firstLine-1
			}else{
				skipVal<-0
			}
			if (yModel=="Categorical") {
				currBetaVector<-scan(betaFile,what=double(),skip=skipVal,n=nFixedEffects*(nCategoriesY-1),quiet=T)
				currBeta<-matrix(currBetaVector,ncol=(nCategoriesY-1),byrow=T)
				currBeta<-cbind(rep(0,dim(currBeta)[1]),currBeta)
			} else {
				currBetaVector<-scan(betaFile,what=double(),skip=skipVal,n=nFixedEffects,quiet=T)
				currBeta<-matrix(currBetaVector,ncol=nCategoriesY,byrow=T)
			}
			betaArray[sweep-firstLine+1,,]<-currBeta
		}
	}

	# Already done the allocation in the C++
	if(doRaoBlackwell){
		# Construct the RB theta file name
		thetaFileName<-file.path(directoryPath,paste(fileStem,'_predictThetaRaoBlackwell.txt',sep=''))
		thetaArray<-array(0,dim=c(nSamples,nPredictSubjects,nCategoriesY))
		# Read the RB theta data
		if (yModel=="Categorical"){
			thetaMat<-matrix(scan(thetaFileName,what=double(),quiet=T),byrow=T,ncol=nPredictSubjects*(nCategoriesY-1))
			for (sweep in firstLine:lastLine){
				thetaArrayRow<-cbind(rep(0,nPredictSubjects),matrix(thetaMat[sweep,],ncol=nCategoriesY-1,byrow=T))

				thetaArray[sweep-firstLine+1,,]<-thetaArrayRow
			}
		} else {
			thetaMat<-matrix(scan(thetaFileName,what=double(),quiet=T),byrow=T,ncol=nPredictSubjects)
			thetaArray[,,1]<-thetaMat[firstLine:nrow(thetaMat),]
		}
	}else{
		# Construct the file names
		zFileName <- file.path(directoryPath,paste(fileStem,'_z.txt',sep=''))
		zFile<-file(zFileName,open="r")
		thetaFileName<-file.path(directoryPath,paste(fileStem,'_theta.txt',sep=''))
		thetaFile<-file(thetaFileName,open="r")
		nClustersFileName<-file.path(directoryPath,paste(fileStem,'_nClusters.txt',sep=''))
		nClustersFile<-file(nClustersFileName,open="r")
		# initialise theta and beta arrays
		thetaArray<-array(0,dim=c(nSamples,nPredictSubjects,nCategoriesY))
		for(sweep in firstLine:lastLine){
			if(sweep==firstLine){
				skipVal<-firstLine-1
			}else{
				skipVal<-0
			}
			currZ<-1+scan(zFile,what=integer(),skip=skipVal,n=nSubjects+nPredictSubjects,quiet=T)
			tmpCurrZ<-currZ
			currZ<-currZ[(nSubjects+1):length(currZ)]
			currNClusters<-scan(nClustersFile,skip=skipVal,n=1,what=integer(),quiet=T)
			if (yModel=="Categorical") {
				currThetaVector<-scan(thetaFile,what=double(),skip=skipVal,n=currNClusters*(nCategoriesY-1),quiet=T)
				currTheta<-matrix(currThetaVector,ncol=(nCategoriesY-1),byrow=T)
				currTheta<-cbind(rep(0,dim(currTheta)[1]),currTheta)
			} else {
				currThetaVector<-scan(thetaFile,what=double(),skip=skipVal,n=currNClusters,quiet=T)
				currTheta<-matrix(currThetaVector,ncol=nCategoriesY,byrow=T)
			}
			thetaArray[sweep-firstLine+1,,]<-as.matrix(currTheta[currZ,],ncol=nCategoriesY)
		}
		close(zFile)
		close(thetaFile)
		close(nClustersFile)
	}
	
	# Use the values of theta to derive the predicted values
	predictedY<-array(0,dim=c(nSamples,nPredictSubjects,nCategoriesY))
	for(sweep in 1:nSamples){
		if(sweep==1||sweep%%1000==0){
			cat(paste("Processing sweep",sweep,"of ",nSamples,"\n"))
		}
	
		lambda<-thetaArray[sweep,,]
		if(fixedEffectsProvided){
			currBeta<-betaArray[sweep,,]
			lambda<-lambda+predictWMat%*%currBeta
		}
	
		if(yModel=='Poisson'){
			if(extraInfoProvided){
				# Add in the offset
				lambda<-lambda+log(predictYMat[,2])
			}
		}
	
		if(yModel=='Bernoulli'){
			predictedY[sweep,,]<-exp(lambda)/(1+exp(lambda))
		}else if(yModel=='Binomial'){
			if(extraInfoProvided){
				predictedY[sweep,,]<-predictYMat[,2]*exp(lambda)/(1+exp(lambda))
			}else{
				predictedY[sweep,,]<-exp(lambda)/(1+exp(lambda))
			}
		}else if(yModel=='Poisson'){
			predictedY[sweep,,]<-exp(lambda)
		}else if(yModel=='Normal'){
			predictedY[sweep,,]<-lambda
		}else if(yModel=="Categorical"){
			predictedY[sweep,,]<-exp(lambda)/rowSums(exp(lambda))
		}
	}

	if(responseProvided){
		bias<-apply(predictedY,2,mean)-predictYMat[,1]
		rmse<-sqrt(mean(bias^2))
		bias<-mean(bias)
		output<-list("bias"=bias,"rmse"=rmse,
			"observedY"=predictYMat[,1],
			"predictedY"=apply(predictedY,c(2,3),mean),
			"doRaoBlackwell"=doRaoBlackwell)
	}else{
		output<-list("bias"=NA,"rmse"=NA,"observedY"=NA,"predictedY"=apply(predictedY,c(2,3),mean,2),"doRaoBlackwell"=doRaoBlackwell)
	}
	if(fullSweepPredictions){
		output$predictedYPerSweep<-predictedY
	}
	if(fullSweepLogOR){
		logORPerSweep<-matrix(0,ncol=ncol(predictedY),nrow=nrow(predictedY))
		for(i in 1:dim(thetaArray)[1]){
			# Always relative to the first prediction subject
			logORPerSweep[i,]<-thetaArray[i,,]-thetaArray[i,1,]
		}
		output$logORPerSweep<-logORPerSweep
	}
	detach(clusObjRunInfoObj)
	detach(riskProfClusObj)
	detach(riskProfObj)
	
	return(output)
}
	
# Compute Ratio of variances (for extra variation case)
computeRatioOfVariance<-function(runInfoObj){
	
	attach(runInfoObj)
	# Construct the number of clusters file name
	nClustersFileName <- file.path(directoryPath,paste(fileStem,'_nClusters.txt',sep=''))
	# Construct the allocation file name
	zFileName <- file.path(directoryPath,paste(fileStem,'_z.txt',sep=''))
	# Construct the allocation file name
	thetaFileName <- file.path(directoryPath,paste(fileStem,'_theta.txt',sep=''))
	# Construct the allocation file name
	epsilonFileName <- file.path(directoryPath,paste(fileStem,'_epsilon.txt',sep=''))
	
	# Restrict to sweeps after burn in
	firstLine<-2+nBurn/nFilter
	lastLine<-1+(nBurn+nSweeps)/nFilter
	
	ratioOfVariance<-rep(0,length(lastLine-firstLine+1))
	for(sweep in firstLine:lastLine){
		currMaxNClusters<-scan(nClustersFileName,what=integer(),skip=sweep-1,n=1,quiet=T)
		zCurr<-1+scan(zFileName,what=integer(),skip=sweep-1,n=nSubjects+nPredictSubjects,quiet=T)
		zCurr<-zCurr[1:nSubjects]
		thetaCurr<-scan(thetaFileName,what=double(),skip=sweep-1,n=currMaxNClusters,quiet=T)
		thetaCurr<-thetaCurr[zCurr]
		vTheta<-var(thetaCurr)
		epsilonCurr<-scan(epsilonFileName,what=double(),skip=sweep-1,n=nSubjects,quiet=T)
		vEpsilon<-var(epsilonCurr)
		ratioOfVariance[sweep-firstLine+1]<-vTheta/(vTheta+vEpsilon)
		
	}
	detach(runInfoObj)
	return(ratioOfVariance)
	
}
	
# Show the continuous hyperparameter for variable selection
summariseVarSelectRho<-function(runInfoObj){
	
	attach(runInfoObj)
	# Rho file name
	rhoFileName <- file.path(directoryPath,paste(fileStem,'_rho.txt',sep=''))
	
	rhoMat<-matrix(scan(rhoFileName,what=double(),quiet=T),ncol=nCovariates,byrow=T)
	
	# Restrict to after burn in
	firstLine<-2+nBurn/nFilter
	lastLine<-1+(nBurn+nSweeps)/nFilter
	
	rhoMat<-rhoMat[firstLine:lastLine,]
	
	rhoMean<-apply(rhoMat,2,mean)
	rhoMedian<-apply(rhoMat,2,median)
	rhoLowerCI<-apply(rhoMat,2,quantile,0.05)
	rhoUpperCI<-apply(rhoMat,2,quantile,0.95)
	
	output<-list("rho"=rhoMat,"rhoMean"=rhoMean,"rhoMedian"=rhoMedian,"rhoLowerCI"=rhoLowerCI,"rhoUpperCI"=rhoUpperCI)
	detach(runInfoObj)
	return(output)
}
	
	
# For a particular clustering, show the cluster values of a separate variable
computeAssociatedVariable<-function(subjectValues,clusObj,clusterPlotOrder,latexFile=NULL){
	
	attach(clusObj)
	outData<-data.frame('cluster'=rep(0,nClusters),'Mean'=rep(NA,nClusters),
		'SD'=rep(NA,nClusters),'Q1'=rep(NA,nClusters),
		'Med'=rep(NA,nClusters),'Q3'=rep(NA,nClusters))
	for(c in 1:nClusters){
		relData<-na.omit(subjectValues[clustering==clusterPlotOrder[c]])
		tmpDF<-data.frame('cluster'=i,'Mean'=mean(relData),
			'SD'=sd(relData),'Q1'=quantile(relData,0.25),
			'Med'=quantile(relData,0.5),'Q3'=quantile(relData,0.75))
		outData[i,]<-tmpDF
	  
	}
	rownames(outData)<-seq(1,nClusters,1)
	  
	if(!is.null(latexFile)){
		outStr<-c()
		for(c in 1:nClusters){
			tmpStr<-paste(format(round(outData[c,],4),nsmall=4),collapse=" & ")
			if(c<nClusters){
				tmpStr<-paste(tmpStr,"\\\\",sep="")
			}
			outStr<-c(outStr,tmpStr)
		}
		write(outStr,file=latexFile,append=F)
	}
	detach(clusObj)
	return(outData)
}
	
compareClustering<-function(riskProfileObjA,riskProfileObjB,clusterOrder=NULL){
	
	clusterObjA<-riskProfileObjA$riskProfClusObj
	clusterObjB<-riskProfileObjB$riskProfClusObj
	disSimMatA<-clusterObjA$clusObjDisSimMat
	disSimMatB<-clusterObjB$clusObjDisSimMat
	runInfoObjA<-clusterObjA$clusObjRunInfoObj
	runInfoObjB<-clusterObjB$clusObjRunInfoObj
	riskMeansA<-apply(riskProfileObjA$risk,2,mean)
	riskMeansB<-apply(riskProfileObjB$risk,2,mean)
	nSubjects<-runInfoObjA$nSubjects
	if(runInfoObjB$nSubjects!=nSubjects){
		exit("Number of subjects not the same in comparison")
	}
	
	optSimObjA<-matrix(0,nSubjects,nSubjects)
	optSimObjB<-matrix(0,nSubjects,nSubjects)
	rawSimObjA<-matrix(0,nSubjects,nSubjects)
	rawSimObjB<-matrix(0,nSubjects,nSubjects)
	r<-1
	dOpt<-0
	dRaw<-0
	dRisk<-0
	dOptClusA<-rep(0,clusterObjA$nClusters)
	dOptClusB<-rep(0,clusterObjB$nClusters)
	for(i in 1:nSubjects){
		if(i<(nSubjects-1)){
			for(j in (i+1):nSubjects){
				optSimObjA[i,j]<-ifelse(clusterObjA$clustering[i]==clusterObjA$clustering[j],1,0)
				optSimObjB[i,j]<-ifelse(clusterObjB$clustering[i]==clusterObjB$clustering[j],1,0)
				rawSimObjA[i,j]<-1-disSimMatA[r]
				rawSimObjB[i,j]<-1-disSimMatB[r]
				r<-r+1
				if(optSimObjA[i,j]!=optSimObjB[i,j]){
					dOpt<-dOpt+1
					dOptClusA[clusterObjA$clustering[i]]<-dOptClusA[clusterObjA$clustering[i]]+1
					dOptClusB[clusterObjB$clustering[i]]<-dOptClusB[clusterObjB$clustering[i]]+1
					dOptClusA[clusterObjA$clustering[j]]<-dOptClusA[clusterObjA$clustering[j]]+1
					dOptClusB[clusterObjB$clustering[j]]<-dOptClusB[clusterObjB$clustering[j]]+1
				}
				dRaw<-dRaw+abs(rawSimObjA[i,j]-rawSimObjB[i,j])
			}
		}
		dRisk<-dRisk+abs(riskMeansA[clusterObjA$clustering[i]]-riskMeansB[clusterObjB$clustering[i]])
	}
	
	dOpt<-2*dOpt/(nSubjects*(nSubjects-1))
	dOptClusA<-dOptClusA/(nSubjects*clusterObjA$clusterSizes)
	dOptClusB<-dOptClusB/(nSubjects*clusterObjB$clusterSizes)
	if(!is.null(clusterOrder)){
		dOptClusA<-dOptClusA[clusterOrder]
		dOptClusB<-dOptClusB[clusterOrder]
	}
	dRaw<-2*dRaw/(nSubjects*(nSubjects-1))
	dRisk<-dRisk/nSubjects
	
	return(list('dOpt'=dOpt,'dOptClusA'=dOptClusA,'dOptClusB'=dOptClusB,'dRaw'=dRaw,'dRisk'=dRisk))
	
}
	


