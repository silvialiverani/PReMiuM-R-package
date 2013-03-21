plotClustering<-function(clusObj,outFile,clusterPlotOrder=NULL,whichCovariates=NULL){
	
	for (i in 1:length(clusObj)) assign(names(clusObj)[i],clusObj[[i]])
	for (i in 1:length(clusObjRunInfoObj)) assign(names(clusObjRunInfoObj)[i],clusObjRunInfoObj[[i]])
	
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

	return(list("rawXDist"=d,"rawXprincipalComponents"=principalComp,"rawXIncludeVec"=includeVec))
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
		stop("Number of subjects not the same in comparison")
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
	
	
# For a particular clustering, show the cluster values of a separate variable
computeAssociatedVariable<-function(subjectValues,clusObj,clusterPlotOrder,latexFile=NULL){
	
	for (i in 1:length(clusObj)) assign(names(clusObj)[i],clusObj[[i]])
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
	return(outData)
}

	
# Compute Ratio of variances (for extra variation case)
computeRatioOfVariance<-function(runInfoObj){
	
	for (i in 1:length(runInfoObj)) assign(names(runInfoObj)[i],runInfoObj[[i]])
	# Construct the number of clusters file name
	nClustersFileName <- file.path(directoryPath,paste(fileStem,'_nClusters.txt',sep=''))
	# Construct the allocation file name
	zFileName <- file.path(directoryPath,paste(fileStem,'_z.txt',sep=''))
	# Construct the allocation file name
	thetaFileName <- file.path(directoryPath,paste(fileStem,'_theta.txt',sep=''))
	# Construct the allocation file name
	epsilonFileName <- file.path(directoryPath,paste(fileStem,'_epsilon.txt',sep=''))
	
	# Restrict to sweeps after burn in
	firstLine<-ifelse(reportBurnIn,nBurn/nFilter+2,1)
	lastLine<-(nSweeps+ifelse(reportBurnIn,nBurn+1,0))/nFilter
	
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
	return(ratioOfVariance)
	
}
	
