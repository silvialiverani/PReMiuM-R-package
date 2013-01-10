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
	
	
