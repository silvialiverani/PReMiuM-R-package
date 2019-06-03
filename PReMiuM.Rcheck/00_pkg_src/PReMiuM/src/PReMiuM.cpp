/// \file PReMiuM.cpp
/// \author David Hastie
/// \brief Main file for running PReMiuM

/// \note (C) Copyright David Hastie and Silvia Liverani, 2012.

/// PReMiuM++ is free software; you can redistribute it and/or modify it under the
/// terms of the GNU Lesser General Public License as published by the Free Software
/// Foundation; either version 3 of the License, or (at your option) any later
/// version.

/// PReMiuM++ is distributed in the hope that it will be useful, but WITHOUT ANY
/// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
/// PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

/// You should have received a copy of the GNU Lesser General Public License
/// along with PReMiuM++ in the documentation directory. If not, see
/// <http://www.gnu.org/licenses/>.

/// The external linear algebra library Eigen, parts of which are included  in the
/// lib directory is released under the LGPL3+ licence. See comments in file headers
/// for details.

/// The Boost C++ header library, parts of which are included in the  lib directory
/// is released under the Boost Software Licence, Version 1.0, a copy of which is
/// included in the documentation directory.


// Standard includes
#include<cmath>
#include<cstdio>
#include<iostream>
#include<numeric>
#include<vector>
#include<string>
#include<cstdlib>
#include<ctime>
#include<sstream>

// Custom includes
#include<MCMC/sampler.h>
#include<MCMC/model.h>
#include<MCMC/proposal.h>
#include<PReMiuMOptions.h>
#include<PReMiuMModel.h>
#include<PReMiuMData.h>
#include<PReMiuMProposals.h>
#include<PReMiuMIO.h>

#include "include/postProcess.h"

#include<Rcpp.h>

using std::vector;
using std::ostringstream;
using std::time;
using std::string;

RcppExport SEXP profRegr(SEXP inputString) {
	string inputStr = Rcpp::as<string>(inputString);

	/* ---------- Start the timer ------------------*/
	time_t beginTime,currTime;
	beginTime = time(NULL);
	/* -----------Process the command line ---------*/
	pReMiuMOptions options = processCommandLine(inputStr);

	/* ---------- Set up the sampler object--------*/
	// Initialise the sampler object
	mcmcSampler<pReMiuMParams,pReMiuMOptions,
				pReMiuMPropParams,pReMiuMData> pReMiuMSampler;

	// Set the options
	pReMiuMSampler.options(options);

	// Set the model
	pReMiuMSampler.model(&importPReMiuMData,&initialisePReMiuM,
							&pReMiuMLogPost,true);

	// Set the missing data function
	pReMiuMSampler.updateMissingDataFn(&updateMissingPReMiuMData);

	// Add the function for writing output
	pReMiuMSampler.userOutputFn(&writePReMiuMOutput);

	// Seed the random number generator
	pReMiuMSampler.seedGenerator(options.seed());

	// Set the sampler specific variables
	pReMiuMSampler.nSweeps(options.nSweeps());
	pReMiuMSampler.nBurn(options.nBurn());
	pReMiuMSampler.nFilter(options.nFilter());
	pReMiuMSampler.nProgress(options.nProgress());
	pReMiuMSampler.reportBurnIn(options.reportBurnIn());

	/* ---------- Read in the data -------- */
	pReMiuMSampler.model().dataset().outcomeType(options.outcomeType());
	pReMiuMSampler.model().dataset().covariateType(options.covariateType());
	pReMiuMSampler.model().dataset().includeCAR(options.includeCAR());
	pReMiuMSampler.importData(options.inFileName(),options.predictFileName(),options.neighbourFileName());
	pReMiuMData dataset = pReMiuMSampler.model().dataset();

	/* ---------- Add the proposals -------- */

	// Set the proposal parameters
	pReMiuMPropParams proposalParams(options.nSweeps(),dataset.nCovariates(),
										dataset.nFixedEffects(),dataset.nCategoriesY());
	pReMiuMSampler.proposalParams(proposalParams);

	// The gibbs update for the active V
	pReMiuMSampler.addProposal("gibbsForVActive",1.0,1,1,&gibbsForVActive);

	if(options.covariateType().compare("Discrete")==0){
		// For discrete X data we do a mixture of Categorical and ordinal updates
		//  Update for the active phi parameters
		pReMiuMSampler.addProposal("updateForPhiActive",1.0,1,1,&updateForPhiActive);

	}else if(options.covariateType().compare("Normal")==0){
		// Need to add the proposals for the normal case
		// Update for the active mu parameters
		if (options.useNormInvWishPrior()){
		// If Inverse Normal Inverse Wishart Prior is used
		pReMiuMSampler.addProposal("gibbsForMuActiveNIWP",1.0,1,1,&gibbsForMuActiveNIWP);
		}else if (options.useIndependentNormal()) {
		// If Independent Normal conditional likelihood is used
		pReMiuMSampler.addProposal("gibbsForMuActiveIndep", 1.0, 1, 1, &gibbsForMuActiveIndep);
		}else{
		// If independant Normal and Inverse Wishart priors are used
		pReMiuMSampler.addProposal("gibbsForMuActive",1.0,1,1,&gibbsForMuActive);
		}

		// Update for the active Sigma parameters
		if (options.useIndependentNormal()) {
		//If the independent normal conditional likelihood is used.
		pReMiuMSampler.addProposal("gibbsForTauActiveIndep", 1.0, 1, 1, &gibbsForTauActiveIndep);
		}else if (options.useSeparationPrior()) {
			pReMiuMSampler.addProposal("gibbsForTauRActive", 1.0, 1, 1, &gibbsForTauRActive);
			pReMiuMSampler.addProposal("metropolisHastingsForTauS", 1.0, 1, 1, &metropolisHastingsForTauS);
		}else {
		//If the multivariate normal conditional likelihood is used.
		pReMiuMSampler.addProposal("gibbsForTauActive", 1.0, 1, 1, &gibbsForTauActive);
		}
		
		// Update for R1
		if(options.useIndependentNormal()){
			pReMiuMSampler.addProposal("gibbsForR1Indep", 1.0, 1, 1, &gibbsForR1Indep);
		}else if (options.useHyperpriorR1()){
			pReMiuMSampler.addProposal("metropolisHastingsForKappa1", 1.0, 1, 1, &metropolisHastingsForKappa1);
			pReMiuMSampler.addProposal("gibbsForR1",1.0,1,1,&gibbsForR1);
			pReMiuMSampler.addProposal("gibbsForMu00", 1.0, 1, 1, &gibbsForMu00);
			pReMiuMSampler.addProposal("gibbsForTau00", 1.0, 1, 1, &gibbsForTau00);
		}
		else if (options.useSeparationPrior()) {
			pReMiuMSampler.addProposal("gibbsForBetaTauS", 1.0, 1, 1, &gibbsForBetaTauS);
			pReMiuMSampler.addProposal("metropolisHastingsForKappa1SP", 1.0, 1, 1, &metropolisHastingsForKappa1SP);
			pReMiuMSampler.addProposal("gibbsForMu00", 1.0, 1, 1, &gibbsForMu00);
		}

	}else if(options.covariateType().compare("Mixed")==0){
		// For discrete X data we do a mixture of Categorical and ordinal updates
		//  Update for the active phi parameters
		pReMiuMSampler.addProposal("updateForPhiActive",1.0,1,1,&updateForPhiActive);

		// Need to add the proposals for the normal case
		// Update for the active mu parameters
		if (options.useNormInvWishPrior()){
		// If Inverse Normal Inverse Wishart Prior is used
		pReMiuMSampler.addProposal("gibbsForMuActiveNIWP",1.0,1,1,&gibbsForMuActiveNIWP);
		}else if (options.useIndependentNormal()) {
			// If Independent Normal conditional likelihood is used
			pReMiuMSampler.addProposal("gibbsForMuActiveIndep", 1.0, 1, 1, &gibbsForMuActiveIndep);
		}else{
		// If independant Normal and Inverse Wishart priors are used
		pReMiuMSampler.addProposal("gibbsForMuActive",1.0,1,1,&gibbsForMuActive);
		}

		// Update for the active Sigma parameters
		if (options.useIndependentNormal()) {
			//If the independent normal conditional likelihood is used.
			pReMiuMSampler.addProposal("gibbsForTauActiveIndep", 1.0, 1, 1, &gibbsForTauActiveIndep);
		}
		else if (options.useSeparationPrior()) {
			pReMiuMSampler.addProposal("gibbsForTauRActive", 1.0, 1, 1, &gibbsForTauRActive);
			pReMiuMSampler.addProposal("metropolisHastingsForTauS", 1.0, 1, 1, &metropolisHastingsForTauS);
		}
		else {
			//If the multivariate normal conditional likelihood is used.
			pReMiuMSampler.addProposal("gibbsForTauActive", 1.0, 1, 1, &gibbsForTauActive);
		}
		
		// Update for R1
		if (options.useIndependentNormal()) {
			pReMiuMSampler.addProposal("gibbsForR1Indep", 1.0, 1, 1, &gibbsForR1Indep);
		}else if (options.useHyperpriorR1()){
			pReMiuMSampler.addProposal("metropolisHastingsForKappa1", 1.0, 1, 1, &metropolisHastingsForKappa1);
			pReMiuMSampler.addProposal("gibbsForR1",1.0,1,1,&gibbsForR1);
			pReMiuMSampler.addProposal("gibbsForMu00", 1.0, 1, 1, &gibbsForMu00);
			pReMiuMSampler.addProposal("gibbsForTau00", 1.0, 1, 1, &gibbsForTau00);
		}else if (options.useSeparationPrior()) {
			pReMiuMSampler.addProposal("gibbsForBetaTauS", 1.0, 1, 1, &gibbsForBetaTauS);
			pReMiuMSampler.addProposal("metropolisHastingsForKappa1SP", 1.0, 1, 1, &metropolisHastingsForKappa1SP);
			pReMiuMSampler.addProposal("gibbsForMu00", 1.0, 1, 1, &gibbsForMu00);
		}
	}

	if(options.varSelectType().compare("None")!=0){
		// Add the variable selection moves
		unsigned int firstSweep;
		firstSweep=1+(unsigned int)(options.nBurn()/10);
		if(options.varSelectType().compare("Continuous")!=0){
			// Gibbs update for gamma
			pReMiuMSampler.addProposal("gibbsForGammaActive",1.0,1,firstSweep,&gibbsForGammaActive);
		}

	}

	if(options.includeResponse()){
		// The Metropolis Hastings update for the active theta
		pReMiuMSampler.addProposal("metropolisHastingsForThetaActive",1.0,1,1,&metropolisHastingsForThetaActive);

		// Adaptive MH for beta
		if(dataset.nFixedEffects()>0){
			pReMiuMSampler.addProposal("metropolisHastingsForBeta",1.0,1,1,&metropolisHastingsForBeta);
		}

		if(options.responseExtraVar()){
			// Adaptive MH for lambda
			pReMiuMSampler.addProposal("metropolisHastingsForLambda",1.0,1,1,&metropolisHastingsForLambda);

			// Gibbs for tauEpsilon
			pReMiuMSampler.addProposal("gibbsForTauEpsilon",1.0,1,1,&gibbsForTauEpsilon);
		}

		//if spatial random term
		if (options.includeCAR()){
		  if(options.outcomeType().compare("Normal")==0){
		    pReMiuMSampler.addProposal("gibbsForUCARNormal", 1.0,1,1,&gibbsForUCARNormal);
		  }
		  if(options.outcomeType().compare("Poisson")==0){
		    //Adaptive rejection sampling for uCAR
		    if (options.PoissonCARadaptive()){
  		    pReMiuMSampler.addProposal("adaptiveRejectionSamplerForUCARPoisson", 1.0,1,1,&adaptiveRejectionSamplerForUCARPoisson);
		    } else {
		    pReMiuMSampler.addProposal("metropolisForUCARPoisson", 1.0,1,1,&metropolisForUCARPoisson);
		    }
		  }
		      
			//Gibbs for TauCAR
			pReMiuMSampler.addProposal("gibbsForTauCAR", 1.0,1,1,&gibbsForTauCAR);
		}

		if(options.outcomeType().compare("Survival")==0){
			if (!options.weibullFixedShape()) {
				// Gibbs for shape parameter (cluster specific) of Weibull for survival response model
				//pReMiuMSampler.addProposal("gibbsForNu",1.0,1,1,&gibbsForNu);
			}
		}
	}


	// Gibbs for U
	if(options.samplerType().compare("Truncated")!=0){
		pReMiuMSampler.addProposal("gibbsForU",1.0,1,1,&gibbsForU);
	}

	// The Metropolis Hastings update for alpha
	if(options.fixedAlpha()<=-1){
		pReMiuMSampler.addProposal("metropolisHastingsForAlpha",1.0,1,1,&metropolisHastingsForAlpha);
	}

	// The gibbs update for the inactive V
	pReMiuMSampler.addProposal("gibbsForVInActive",1.0,1,1,&gibbsForVInActive);

	if(options.covariateType().compare("Discrete")==0){
		// For discrete X data we do a mixture of Categorical and ordinal updates
		//  Update for the inactive phi parameters
		pReMiuMSampler.addProposal("gibbsForPhiInActive",1.0,1,1,&gibbsForPhiInActive);

	}else if(options.covariateType().compare("Normal")==0){
		// Need to add the proposals for the normal case
		// Update for the active mu parameters
		if (options.useNormInvWishPrior()){
		// If Inverse Normal Inverse Wishart Prior is used
		pReMiuMSampler.addProposal("gibbsForMuInActiveNIWP",1.0,1,1,&gibbsForMuInActiveNIWP);
		}else if (options.useIndependentNormal()) {
		// If the independent normal conditional likelihood is used
		pReMiuMSampler.addProposal("gibbsForMuInActiveIndep", 1.0, 1, 1, &gibbsForMuInActiveIndep);
		}else{
		// If independant Normal and Inverse Wishart priors are used
		pReMiuMSampler.addProposal("gibbsForMuInActive",1.0,1,1,&gibbsForMuInActive);
		}


		// Update for the active Sigma parameters
		if (options.useIndependentNormal()) {
			// If the independent normal conditional likelihood is used
			pReMiuMSampler.addProposal("gibbsForTauInActiveIndep", 1.0, 1, 1, &gibbsForTauInActiveIndep);
		}
		else if (options.useSeparationPrior()) {
			pReMiuMSampler.addProposal("gibbsForTauRInActive", 1.0, 1, 1, &gibbsForTauRInActive);
			pReMiuMSampler.addProposal("gibbsForTauSInActive", 1.0, 1, 1, &gibbsForTauSInActive);

		}else {
			// If the multivariate normal conditional likelihood is used
			pReMiuMSampler.addProposal("gibbsForTauInActive", 1.0, 1, 1, &gibbsForTauInActive);
		}
		

	}else if(options.covariateType().compare("Mixed")==0){

		// For discrete X data we do a mixture of Categorical and ordinal updates
		//  Update for the inactive phi parameters
		pReMiuMSampler.addProposal("gibbsForPhiInActive",1.0,1,1,&gibbsForPhiInActive);

		// Need to add the proposals for the normal case
		// Update for the active mu parameters
		if (options.useNormInvWishPrior()){
		// If Inverse Normal Inverse Wishart Prior is used
		pReMiuMSampler.addProposal("gibbsForMuInActiveNIWP",1.0,1,1,&gibbsForMuInActiveNIWP);
		}else if (options.useIndependentNormal()) {
		// If the independent normal conditional likelihood is used
		pReMiuMSampler.addProposal("gibbsForMuInActiveIndep", 1.0, 1, 1, &gibbsForMuInActiveIndep);
		}else{
		// If independant Normal and Inverse Wishart priors are used
		pReMiuMSampler.addProposal("gibbsForMuInActive",1.0,1,1,&gibbsForMuInActive);
		}


		// Update for the active Sigma parameters
		if (options.useIndependentNormal()) {
			// If the independent normal conditional likelihood is used
			pReMiuMSampler.addProposal("gibbsForTauInActiveIndep", 1.0, 1, 1, &gibbsForTauInActiveIndep);
		}
		else if (options.useSeparationPrior()) {
			pReMiuMSampler.addProposal("gibbsForTauRInActive", 1.0, 1, 1, &gibbsForTauRInActive);
			pReMiuMSampler.addProposal("gibbsForTauSInActive", 1.0, 1, 1, &gibbsForTauSInActive);

		}
		else {
			// If the multivariate normal conditional likelihood is used
			pReMiuMSampler.addProposal("gibbsForTauInActive", 1.0, 1, 1, &gibbsForTauInActive);
		}
		
	}

	if(options.varSelectType().compare("None")!=0){
		// Add the variable selection moves
		unsigned int firstSweep;
		firstSweep=1+(unsigned int)(options.nBurn()/10);
		if(options.varSelectType().compare("Continuous")!=0){
			// Gibbs update for gamma
			pReMiuMSampler.addProposal("gibbsForGammaInActive",1.0,1,firstSweep,&gibbsForGammaInActive);
		}

	}

	// The Metropolis Hastings update for labels
	if(options.whichLabelSwitch().compare("123")==0){
		pReMiuMSampler.addProposal("metropolisHastingsForLabels123",1.0,1,1,&metropolisHastingsForLabels123);
	} else if (options.whichLabelSwitch().compare("12")==0){
		pReMiuMSampler.addProposal("metropolisHastingsForLabels12",1.0,1,1,&metropolisHastingsForLabels12);
	} else if(options.whichLabelSwitch().compare("3")==0){
		pReMiuMSampler.addProposal("metropolisHastingsForLabels3",1.0,1,1,&metropolisHastingsForLabels3);
	} 


	if(options.includeResponse()){
		// The Metropolis Hastings update for the inactive theta
		pReMiuMSampler.addProposal("gibbsForThetaInActive",1.0,1,1,&gibbsForThetaInActive);
		if(options.outcomeType().compare("Survival")==0&&!options.weibullFixedShape()) {
			pReMiuMSampler.addProposal("gibbsForNu",1.0,1,1,&gibbsForNu);
			pReMiuMSampler.addProposal("gibbsForNuInActive",1.0,1,1,&gibbsForNuInActive);

		}
	}

	if(options.varSelectType().compare("None")!=0){
		// Add the variable selection moves
		// Metropolis Hastings for joint update of rho and omega
		unsigned int firstSweep;
		firstSweep=1+(unsigned int)(options.nBurn()/10);

		pReMiuMSampler.addProposal("metropolisHastingsForRhoOmega",1.0,1,firstSweep,&metropolisHastingsForRhoOmega);
	}


	if(options.includeResponse()){
		if(options.outcomeType().compare("Normal")==0){
			// Gibbs for sigmaSqY for Normal response model
			pReMiuMSampler.addProposal("gibbsForSigmaSqY",1.0,1,1,&gibbsForSigmaSqY);
		}
                if(options.outcomeType().compare("Quantile")==0){
			// Gibbs for sigmaSqY for Quantile response model
			pReMiuMSampler.addProposal("gibbsForSigmaSqYQuantile",1.0,1,1,&gibbsForSigmaSqYQuantile);
		}
		if(options.outcomeType().compare("Survival")==0){
      			if (options.weibullFixedShape()) {
				// Gibbs for shape parameter of Weibull for survival response model
				pReMiuMSampler.addProposal("gibbsForNu",1.0,1,1,&gibbsForNu);
			}
		}
	}


	// Gibbs update for the allocation parameters
	pReMiuMSampler.addProposal("gibbsForZ",1.0,1,1,&gibbsForZ);


	/* ---------- Initialise the output files -----*/
	pReMiuMSampler.initialiseOutputFiles(options.outFileStem());

	/* ---------- Write the log file ------------- */
	// The standard log file
	pReMiuMSampler.writeLogFile();

	/* ---------- Initialise the chain ---- */
	pReMiuMSampler.initialiseChain();

	pReMiuMHyperParams hyperParams = pReMiuMSampler.chain().currentState().parameters().hyperParams();
	unsigned int nClusInit = pReMiuMSampler.chain().currentState().parameters().workNClusInit();

	// The following is only used if the sampler type is truncated
	unsigned int maxNClusters = pReMiuMSampler.chain().currentState().parameters().maxNClusters();
	/* ---------- Run the sampler --------- */
	// Note: in this function the output gets written
	pReMiuMSampler.run();

	/* -- End the clock time and write the full run details to log file --*/
	currTime = time(NULL);
    	double timeInSecs=(double)currTime-(double)beginTime;
	string tmpStr = storeLogFileData(options,dataset,hyperParams,nClusInit,maxNClusters,timeInSecs);
	pReMiuMSampler.appendToLogFile(tmpStr);

	/* ---------- Clean Up ---------------- */
	pReMiuMSampler.closeOutputFiles();

	//int err = 0;
	return Rcpp::wrap(0);
	// alternative output
	// return Rcpp::List::create(Rcpp::Named("yModel")=options.outcomeType());


}

