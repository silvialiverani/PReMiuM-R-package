/// \file PReMiuMIO.h
/// \author David Hastie
/// \brief Header file for handling input and output for PReMiuM

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
/// is released under the Boost Software Licence, Version 1.0, a copy  of which is
/// included in the documentation directory.


#ifndef DIPBACIO_H_
#define DIPBACIO_H_

// Standard includes
#include<cmath>
#include<vector>
#include<iostream>
#include<fstream>
#include<string>
#include<numeric>
#include<limits>
#include<cstdlib>
#include<sstream>
#include<algorithm>
#include<iterator>
#include<cstdint>

#include<Eigen/Core>
#include<Eigen/Cholesky>
#include<Eigen/LU>

#include<Rcpp.h>

// Custom includes
#include<MCMC/chain.h>
#include<MCMC/model.h>
#include<MCMC/sampler.h>
#include<PReMiuMOptions.h>
#include<PReMiuMData.h>
#include<PReMiuMModel.h>

using namespace Eigen;

using std::vector;
using std::ostringstream;
using std::istringstream;
using std::string;
using std::endl;
using std::stringstream;

// Process the command line run time options
pReMiuMOptions processCommandLine(string inputStr){

	/* ---------- Handle run time options ----------*/
	pReMiuMOptions options;

	vector<string> inputStrings;
	istringstream iss(inputStr);
	copy(std::istream_iterator<string>(iss),
		std::istream_iterator<string>(),
		std::back_inserter<vector<string> >(inputStrings));

	int argc = inputStrings.size();

	int currArg=1;
	bool wasError=false;
	if(argc>1){
		string inString = inputStrings[currArg];
		if(inString.compare("--help")==0){
			// Output help if requested
			Rprintf("\n\n ### PReMiuMpp Help Page. ### \n\n");
			Rprintf("Possible arguments (defaults in parentheses):\n\n");
			Rprintf("--help\n\tShow this help page\n");
			Rprintf("--input=<string>\n\tThe full file path to the input data (./input.txt)\n");
			Rprintf("--output=<string>\n\tThe file stem (including full file path) where the data\n\tshould be written (./output)\n");
			Rprintf("--hyper=<string>\n\tThe full file path to the file containing hyper\n\t parameters. (Hyper parameter file not used)\n");
			Rprintf("--predict=<string>\n\tThe full file path to the file containing prediction\n\t covariates. (Prediction file not used)\n");
			Rprintf("--nSweeps=<unsigned int>\n\tThe number of sweeps (after burn in) to run the\n\t sampler for (10000)\n");
			Rprintf("--nBurn=<unsigned int>\n\tThe number of sweeps in the burn in period (1000)\n");
			Rprintf("--reportBurnIn=<bool>\n\tIt enables reporting in the output files of the burn-in period (true).\n");
			Rprintf("--nProgress=<unsigned int>\n\tThe number of sweeps at which to print a\nprogress update (500)\n");
			Rprintf("--nFilter=<unsigned int>\n\tThe frequency (in sweeps) with which to write\n\tthe output to file (1)\n");
			Rprintf("--nClusInit=<unsigned int>\n\tThe number of clusters individuals should be\n\tinitially randomly assigned to (Unif[50,60])\n");
			Rprintf("--seed=<unsigned int>\n\tThe value for the seed for the random number\n\tgenerator (current time)\n");
			Rprintf("--yModel=<string>\n\tThe model type for the outcome variable. Options are\n\tcurrently 'Bernoulli','Poisson','Binomial', 'Categorical', 'Survival', 'Normal' and 'Quantile' (Bernoulli)\n");
			Rprintf("--xModel=<string>\n\tThe model type for the covariates. Options are\n\tcurrently 'Discrete', 'Normal' and 'Mixed' (Discrete)\n");
			Rprintf("--sampler=<string>\n\tThe sampler type to be used. Options are\n\tcurrently 'SliceDependent', 'SliceIndependent' and 'Truncated' (SliceDependent)\n");
			Rprintf("--alpha=<double>\n\tThe value to be used if alpha is to remain fixed.\n\tIf a negative value is used then alpha is updated (-2)\n");
			Rprintf("--dPitmanYor=<double>\n\tThe value to be used for the discount parameter of the Pitman-Yor process prior.\n\tThe default corresponds to the Dirichlet process prior (0)\n");
			Rprintf("--pQuantile=<double>\n\tThe quantile for the Quantile y-model.\n");
			Rprintf("--excludeY\n\tIf included only the covariate data X is modelled (not included)\n");
			Rprintf("--extraYVar\n\tIf included extra Gaussian variance is included in the\n\tresponse model (not included).\n");
			Rprintf("--varSelect=<string>\n\tThe type of variable selection to be used 'None',\n\t'BinaryCluster' or 'Continuous' (None)\n");
			Rprintf("--entropy\n\tIf included then we compute allocation entropy (not included)\n");
			Rprintf("--predictType=<string>\n\tThe type of predictions to be used 'RaoBlackwell' or 'random' (RaoBlackwell)\n");
			Rprintf("--weibullFixedShape=<bool>\n\tWhether the shape parameter of the Weibull distribution is fixed.\n");
			Rprintf("--PoissonCARadaptive=<bool>\n\tWhether the adaptive rejection sampler or the random walk Metropolis are used for the Poisson outcome with spatial random effects.\n");
		}else{
			while(currArg < argc){
				inString.assign(inputStrings[currArg]);
				if(inString.find("--input")!=string::npos){
					size_t pos = inString.find("=")+1;
					string inFileName = inString.substr(pos,inString.size()-pos);
					options.inFileName(inFileName);
				}else if(inString.find("--output")!=string::npos){
					size_t pos = inString.find("=")+1;
					string outFileStem = inString.substr(pos,inString.size()-pos);
					options.outFileStem(outFileStem);
				}else if(inString.find("--hyper")!=string::npos){
					size_t pos = inString.find("=")+1;
					string hyperParamFileName = inString.substr(pos,inString.size()-pos);
					options.hyperParamFileName(hyperParamFileName);
				}else if(inString.find("--predict")!=string::npos){
					size_t pos = inString.find("=")+1;
					string predictFileName = inString.substr(pos,inString.size()-pos);
					options.predictFileName(predictFileName);
				}else if(inString.find("--nSweeps")!=string::npos){
					size_t pos = inString.find("=")+1;
					string tmpStr = inString.substr(pos,inString.size()-pos);
					unsigned int nSweeps = (unsigned int)atoi(tmpStr.c_str());
					options.nSweeps(nSweeps);
				}else if(inString.find("--nBurn")!=string::npos){
					size_t pos = inString.find("=")+1;
					string tmpStr = inString.substr(pos,inString.size()-pos);
					unsigned int nBurn=(unsigned int)atoi(tmpStr.c_str());
					options.nBurn(nBurn);
				}else if(inString.find("--reportBurnIn")!=string::npos){
					options.reportBurnIn(true);
				}else if(inString.find("--nProgress")!=string::npos){
					size_t pos = inString.find("=")+1;
					string tmpStr = inString.substr(pos,inString.size()-pos);
					unsigned int nProgress=(unsigned int)atoi(tmpStr.c_str());
					options.nProgress(nProgress);
				}else if(inString.find("--nFilter")!=string::npos){
					size_t pos = inString.find("=")+1;
					string tmpStr = inString.substr(pos,inString.size()-pos);
					unsigned int nFilter=(unsigned int)atoi(tmpStr.c_str());
					options.nFilter(nFilter);
				}else if(inString.find("--nClusInit")!=string::npos){
					size_t pos = inString.find("=")+1;
					string tmpStr = inString.substr(pos,inString.size()-pos);
					unsigned int nClusInit = (unsigned int)atoi(tmpStr.c_str());
					options.nClusInit(nClusInit);
				}else if(inString.find("--seed")!=string::npos){
					size_t pos = inString.find("=")+1;
					string tmpStr = inString.substr(pos,inString.size()-pos);
					uint_fast32_t rndSeed=(uint_fast32_t)atoi(tmpStr.c_str());
					options.seed(rndSeed);
				}else if(inString.find("--yModel")!=string::npos){
					size_t pos = inString.find("=")+1;
					string outcomeType = inString.substr(pos,inString.size()-pos);
					if(outcomeType.compare("Poisson")!=0&&outcomeType.compare("Bernoulli")!=0&&
							outcomeType.compare("Categorical")!=0&&outcomeType.compare("Survival")!=0&&
							outcomeType.compare("Binomial")!=0&&
							outcomeType.compare("Normal")!=0&&outcomeType.compare("Quantile")!=0){
						// Illegal outcome model entered
						wasError=true;
						break;
					}
					options.outcomeType(outcomeType);
					if(outcomeType.compare("Normal")==0&&options.responseExtraVar()){
						Rprintf("Response extra variation not permitted with Normal response\n");			
						options.responseExtraVar(false);
					}
					if(outcomeType.compare("Survival")==0&&options.responseExtraVar()){
						Rprintf("Response extra variation not permitted with Survival response\n");			
						options.responseExtraVar(false);
					}
					if(outcomeType.compare("Quantile")==0&&options.responseExtraVar()){
						Rprintf("Response extra variation not permitted with Quantile response\n");			
						options.responseExtraVar(false);
					}
				}else if(inString.find("--xModel")!=string::npos){
					size_t pos = inString.find("=")+1;
					string covariateType = inString.substr(pos,inString.size()-pos);
					if(covariateType.compare("Discrete")!=0&&covariateType.compare("Normal")!=0&&covariateType.compare("Mixed")!=0){
						// Illegal covariate type entered
						wasError=true;
						break;
					}
					options.covariateType(covariateType);
				}else if(inString.find("--whichLabelSwitch")!=string::npos){
					size_t pos = inString.find("=")+1;
					string whichLabelSwitch = inString.substr(pos,inString.size()-pos);
					if(whichLabelSwitch.compare("123")!=0&&whichLabelSwitch.compare("12")!=0&&whichLabelSwitch.compare("3")!=0){
						// Illegal covariate type entered
						wasError=true;
						break;
					}
					options.whichLabelSwitch(whichLabelSwitch);
				}else if(inString.find("--sampler")!=string::npos){
					size_t pos = inString.find("=")+1;
					string samplerType = inString.substr(pos,inString.size()-pos);
					if(samplerType.compare("SliceDependent")!=0&&samplerType.compare("SliceIndependent")!=0
							&&samplerType.compare("Truncated")!=0){
						// Illegal sampler type entered
						wasError=true;
						break;
					}
					options.samplerType(samplerType);
				}else if(inString.find("--alpha")!=string::npos){
					size_t pos = inString.find("=")+1;
					string tmpStr = inString.substr(pos,inString.size()-pos);
					double alpha=(double)atof(tmpStr.c_str());
					options.fixedAlpha(alpha);
				}else if(inString.find("--dPitmanYor")!=string::npos){
					size_t pos = inString.find("=")+1;
					string tmpStr = inString.substr(pos,inString.size()-pos);
					double dPitmanYor=(double)atof(tmpStr.c_str());
					options.dPitmanYor(dPitmanYor);
				}else if(inString.find("--excludeY")!=string::npos){
					options.includeResponse(false);
                		}else if(inString.find("--includeCAR")!=string::npos){
					options.includeCAR(true);
        		        }else if(inString.find("--neighbours")!=string::npos){
					size_t pos = inString.find("=")+1;
					string neighboursFile = inString.substr(pos,inString.size()-pos);
					options.neighbourFileName(neighboursFile);
        		        }else if(inString.find("--uCARinit")!=string::npos){
					options.includeuCARinit(true);
					size_t pos = inString.find("=")+1;
					string uCARinitFile = inString.substr(pos,inString.size()-pos);
					options.uCARinitFileName(uCARinitFile);
				}else if(inString.find("--extraYVar")!=string::npos){
					if(options.outcomeType().compare("Normal")!=0){
						options.responseExtraVar(true);
					}else{
						Rprintf("Response extra variation not permitted with Normal response\n");
					}
					if(options.outcomeType().compare("Survival")==0) Rprintf("Response extra variation not permitted with Survival response\n");
					if(options.outcomeType().compare("Quantile")==0) Rprintf("Response extra variation not permitted with Quantile response\n");

				}else if(inString.find("--varSelect")!=string::npos){
					size_t pos = inString.find("=")+1;
					string varSelectType = inString.substr(pos,inString.size()-pos);
					if(varSelectType.compare("None")!=0&&
						varSelectType.compare("BinaryCluster")!=0&&varSelectType.compare("Continuous")!=0){
						// Illegal type for variable selection entered
						wasError=true;
						break;
					}
					options.varSelectType(varSelectType);
				}else if(inString.find("--entropy")!=string::npos){
					options.computeEntropy(true);
				}else if(inString.find("--predType")!=string::npos){
					size_t pos = inString.find("=")+1;
					string predictType = inString.substr(pos,inString.size()-pos);
					if(predictType.compare("RaoBlackwell")!=0&&predictType.compare("random")!=0){
						// Illegal predictType type entered
						wasError=true;
						break;
					}
					options.predictType(predictType);
				}else if(inString.find("--PoissonCARadaptive")!=string::npos){
					options.PoissonCARadaptive(true);
				}else if(inString.find("--weibullFixedShape")!=string::npos){
					options.weibullFixedShape(true);
		                }else if(inString.find("--useNormInvWishPrior")!=string::npos){
					options.useNormInvWishPrior(true);
		                }else if(inString.find("--useHyperpriorR1")!=string::npos){
					options.useHyperpriorR1(true);
						}else if (inString.find("--useSeparationPrior") != string::npos) {
					options.useSeparationPrior(true);
						}else if (inString.find("--useIndependentNormal") != string::npos) {
					options.useIndependentNormal(true);
				}else{
					Rprintf("Unknown command line option.\n");
					wasError=true;
					break;
				}
				currArg++;
			}
		}
	}

	// Return if there was an error
	if(wasError){
		Rprintf("There is a mistake in the arguments provided in profRegr.\n");
		Rprintf("The code will be run with default values.\n");
	//	Rprintf("Please use:\n");
	//	Rprintf("\t profileRegression --help\n");
	//	Rprintf("to get help on correct usage.\n");
	//	exit(-1);
	}

	return options;

}

// Read the PReMiuM data set
void importPReMiuMData(const string& fitFilename,const string& predictFilename, const string& neighboursFilename, pReMiuMData& dataset){

	ifstream inputFile,predictFile;
	inputFile.open(fitFilename.c_str());
	if(!inputFile.is_open()){
		Rprintf("Input file not found\n");
	//	exit(-1);
	}
	if(predictFilename.compare("")!=0){
		predictFile.open(predictFilename.c_str());
		if(!predictFile.is_open()){
			Rprintf("Prediction covariate file not found\n");
	//		exit(-1);
		}
	}
	unsigned int& nSubjects=dataset.nSubjects();
	unsigned int& nCovariates=dataset.nCovariates();
	unsigned int& nDiscreteCovs=dataset.nDiscreteCovs();
	unsigned int& nContinuousCovs=dataset.nContinuousCovs();
	unsigned int& nFixedEffects=dataset.nFixedEffects();
	unsigned int& nCategoriesY=dataset.nCategoriesY();
	unsigned int& nPredictSubjects=dataset.nPredictSubjects();
	vector<unsigned int>& nCategories=dataset.nCategories();
	vector<unsigned int>& discreteY=dataset.discreteY();
	vector<double>& continuousY=dataset.continuousY();
	vector<vector<int> >& discreteX=dataset.discreteX();
	vector<vector<double> >& continuousX=dataset.continuousX();
	vector<string>& covNames=dataset.covariateNames();
	vector<vector<bool> >& missingX=dataset.missingX();
	vector<unsigned int>& nContinuousCovariatesNotMissing=dataset.nContinuousCovariatesNotMissing();
	vector<vector<double> >& W=dataset.W();
	vector<string>& confNames=dataset.fixedEffectNames();
	string outcomeType = dataset.outcomeType();
	string covariateType = dataset.covariateType();
	vector<double>& logOffset=dataset.logOffset();
	vector<unsigned int>& nTrials=dataset.nTrials();
	vector<unsigned int>& censoring=dataset.censoring();
	vector<vector<unsigned int> >& neighbours=dataset.neighbours();
	vector<unsigned int>& nNeighbours=dataset.nNeighbours();
	bool& includeCAR=dataset.includeCAR();

	bool wasError=false;

	// Get the number of subjects
	inputFile >> nSubjects;
	// Get the number of covariates
	inputFile >> nCovariates;
	covNames.resize(nCovariates);
	if(covariateType.compare("Mixed")==0){
		inputFile >> nDiscreteCovs;
		inputFile >> nContinuousCovs;
		if(nDiscreteCovs+nContinuousCovs!=nCovariates){
				Rprintf("Illegal number of covariates, discrete covariates or continuous covariates\n");
				// Illegal number of covariates, discrete covariates or continuous covariates
				wasError=true;
		}
		if(nDiscreteCovs==0 || nContinuousCovs==0){
				Rprintf("If xModel=Mixed a positive number of discrete and continuous covariates must be provided\n");
				// Illegal number of discrete covariates or continuous covariates
				wasError=true;
		}
	} else {
		nDiscreteCovs = 0;
		nContinuousCovs = 0;
	}
	for(unsigned int i=0;i<nCovariates;i++){
		inputFile >> covNames[i];
	}
	// Get the number of fixed effects
	inputFile >> nFixedEffects;
	confNames.resize(nFixedEffects);
	for(unsigned int i=0;i<nFixedEffects;i++){
		inputFile >> confNames[i];
	}
	// Get the number of categories of outcome Y
	if(outcomeType.compare("Categorical")==0){
		inputFile >> nCategoriesY;
		nCategoriesY--;
	} else {
		nCategoriesY=1;
	}

	nCategories.resize(nCovariates);
	if(covariateType.compare("Discrete")==0){
		// Get the number of categories for each covariate
		for(unsigned int j=0;j<nCovariates;j++){
			inputFile >> nCategories[j];
		}
	}else if(covariateType.compare("Normal")==0){
		for(unsigned int j=0;j<nCovariates;j++){
			nCategories[j]=0;
		}
	}else if(covariateType.compare("Mixed")==0){
		for(unsigned int j=0;j<nDiscreteCovs;j++){
			inputFile >> nCategories[j];
		}
		for(unsigned int j=nDiscreteCovs;j<nCovariates;j++){
			nCategories[j]=0;
		}
	}

	if(predictFile.is_open()){
		predictFile >> nPredictSubjects;
	}

	// Get the data
	discreteY.resize(nSubjects);
	continuousY.resize(nSubjects);
	discreteX.resize(nSubjects+nPredictSubjects);
	continuousX.resize(nSubjects+nPredictSubjects);
	W.resize(nSubjects);
	if(outcomeType.compare("Poisson")==0){
		logOffset.resize(nSubjects);
	}
	if(outcomeType.compare("Binomial")==0){
		nTrials.resize(nSubjects);
	}
	if(outcomeType.compare("Survival")==0){
		censoring.resize(nSubjects);
	}
	missingX.resize(nSubjects+nPredictSubjects);
	nContinuousCovariatesNotMissing.resize(nSubjects+nPredictSubjects);
	vector<double> meanX(nCovariates,0);
	vector<unsigned int> nXNotMissing(nCovariates,0);
	for(unsigned int i=0;i<nSubjects;i++){
		if(outcomeType.compare("Normal")==0||outcomeType.compare("Survival")==0||outcomeType.compare("Quantile")==0){
			inputFile >> continuousY[i];
		}else{
			inputFile >> discreteY[i];
		}
		if(covariateType.compare("Discrete")==0 || covariateType.compare("Normal")==0){
			discreteX[i].resize(nCovariates);
			continuousX[i].resize(nCovariates);
		} else if(covariateType.compare("Mixed")==0){
			discreteX[i].resize(nDiscreteCovs);
			continuousX[i].resize(nContinuousCovs);
		}
		missingX[i].resize(nCovariates);
		for(unsigned int j=0;j<nCovariates;j++){
			missingX[i][j]=true;
			if(covariateType.compare("Discrete")==0){
				inputFile >> discreteX[i][j];
				// -999 is missing data indicator
				if(discreteX[i][j]!=-999){
					meanX[j]+=(double)discreteX[i][j];
					nXNotMissing[j]+=1;
					missingX[i][j]=false;
				}
			}else if(covariateType.compare("Normal")==0){
				inputFile >> continuousX[i][j];
				// -999 is missing data indicator
				if(fabs(continuousX[i][j]+999)>0.00000000001){
					nContinuousCovariatesNotMissing[i]++;
					meanX[j]+=continuousX[i][j];
					nXNotMissing[j]+=1;
					missingX[i][j]=false;
				}
			}else if(covariateType.compare("Mixed")==0){
				if (j < nDiscreteCovs) {
					inputFile >> discreteX[i][j];
					if(discreteX[i][j]!=-999){
						meanX[j]+=(double)discreteX[i][j];
						nXNotMissing[j]+=1;
						missingX[i][j]=false;
					}
				} else {
					inputFile >> continuousX[i][j-nDiscreteCovs];
					if(fabs(continuousX[i][j-nDiscreteCovs]+999)>0.00000000001){
						nContinuousCovariatesNotMissing[i]++;
						meanX[j]+=continuousX[i][j-nDiscreteCovs];
						nXNotMissing[j]+=1;
						missingX[i][j]=false;
					}
				}
			}
		}
		W[i].resize(nFixedEffects);
		for(unsigned int j=0;j<nFixedEffects;j++){
			inputFile >> W[i][j];
		}
		if(outcomeType.compare("Poisson")==0){
			double tmp;
			inputFile >> tmp;
			logOffset[i]=log(tmp);
		}
		if(outcomeType.compare("Binomial")==0){
			inputFile >> nTrials[i];
		}
		if(outcomeType.compare("Survival")==0){
			inputFile >> censoring[i];
		}
	}

	for(unsigned int i=nSubjects;i<nSubjects+nPredictSubjects;i++){
		if(covariateType.compare("Discrete")==0 || covariateType.compare("Normal")==0){
			discreteX[i].resize(nCovariates);
			continuousX[i].resize(nCovariates);
		} else if(covariateType.compare("Mixed")==0){
			discreteX[i].resize(nDiscreteCovs);
			continuousX[i].resize(nContinuousCovs);
		}
		missingX[i].resize(nCovariates);
		for(unsigned int j=0;j<nCovariates;j++){
			missingX[i][j]=true;
			if(covariateType.compare("Discrete")==0){
				predictFile >> discreteX[i][j];
				// -999 is missing data indicator
				if(discreteX[i][j]!=-999){
					missingX[i][j]=false;
				}
			}else if(covariateType.compare("Normal")==0){
				predictFile >> continuousX[i][j];
				// -999 is missing data indicator
				if(fabs(continuousX[i][j]+999)>0.00000000001){
					nContinuousCovariatesNotMissing[i]++;
					missingX[i][j]=false;
				}
			}else if(covariateType.compare("Mixed")==0){
				if (j < nDiscreteCovs) {
					predictFile >> discreteX[i][j];
					if(discreteX[i][j]!=-999){
						missingX[i][j]=false;
					}
				} else {
					predictFile >> continuousX[i][j-nDiscreteCovs];
					if(fabs(continuousX[i][j-nDiscreteCovs]+999)>0.00000000001){
						nContinuousCovariatesNotMissing[i]++;
						missingX[i][j]=false;
					}
				}
			}
		}
	}


	/// Initially we just replace missing values by their means
	for(unsigned int j=0;j<nCovariates;j++){
		meanX[j]=meanX[j]/(double)nXNotMissing[j];
	}

	for(unsigned int i=0;i<nSubjects+nPredictSubjects;i++){
		for(unsigned int j=0;j<nCovariates;j++){
			if(missingX[i][j]){
				if(covariateType.compare("Discrete")==0){
					discreteX[i][j]=(int)meanX[j];
				}else if(covariateType.compare("Normal")==0){
					continuousX[i][j]=meanX[j];
				}else if(covariateType.compare("Mixed")==0){
					if (j < nDiscreteCovs) {
						discreteX[i][j]=(int)meanX[j];
					} else {
						continuousX[i][j-nDiscreteCovs]=(double)meanX[j-nDiscreteCovs];
					}
				}
			}
		}
	}

	inputFile.close();
	if(predictFile.is_open()){
		predictFile.close();
	}

	//Fill nNeighbours and Neighbours
	if (includeCAR){
        ifstream neighFile;
        neighFile.open(neighboursFilename.c_str());

        if (!neighFile.is_open()){
            Rprintf("Neighbourhood structure file not found\n");
		wasError = true;
        }
        if (neighFile.good()){
            string line;
            getline(neighFile, line);
            stringstream streamline(line);
            unsigned int nsub;
            streamline >> nsub;
            nNeighbours.resize(nSubjects);
            neighbours.resize(nSubjects);
        }
        int i=0;
        while (neighFile.good()){
            string line;
            getline(neighFile, line);
            stringstream streamline(line);
            int j;
            streamline>>j;
            streamline>>nNeighbours[j-1];
            neighbours[j-1].resize(nNeighbours[j-1]);
            int k=0;
            while (streamline.good()){
                streamline>>neighbours[j-1][k];
                k++;
            }
          i++;
        }
        neighFile.close();
	}

	// Return if there was an error
	if(wasError){
		Rprintf("Please use:\n");
		Rprintf("\t profileRegression --help\n");
		Rprintf("to get help on correct usage.\n");
	//	exit(-1);
	}

}

// Function to read the hyper parameters from file
void readHyperParamsFromFile(const string& filename,pReMiuMHyperParams& hyperParams){

	ifstream inputFile;
	inputFile.open(filename.c_str());
	if(!inputFile.is_open()){
		Rprintf("Parameter file not found\n");
	//	exit(-1);
	}

	string inString;

	bool wasError=false;

	while(!inputFile.eof()){
		getline(inputFile,inString);
		if(inString.find("shapeAlpha")==0){
			size_t pos = inString.find("=")+1;
			string tmpStr = inString.substr(pos,inString.size()-pos);
			double shapeAlpha = (double)atof(tmpStr.c_str());
			hyperParams.shapeAlpha(shapeAlpha);
		}else if(inString.find("rateAlpha")==0){
			size_t pos = inString.find("=")+1;
			string tmpStr = inString.substr(pos,inString.size()-pos);
			double rateAlpha = (double)atof(tmpStr.c_str());
			hyperParams.rateAlpha(rateAlpha);

//		}else if(inString.find("useReciprocalNCatsPhi")==0){
//			size_t pos = inString.find("=")+1;
//			string tmpStr = inString.substr(pos,inString.size()-pos);
//			bool useRecip = false;
//			if(tmpStr.compare("true")==0){
//				useRecip = true;
//			}
//			hyperParams.useReciprocalNCatsPhi(useRecip);
		}
		else if (inString.find("shapeKappa1") == 0) {
			size_t pos = inString.find("=") + 1;
			string tmpStr = inString.substr(pos, inString.size() - pos);
			double shapeKappa1 = (double)atof(tmpStr.c_str());
			hyperParams.shapeKappa1(shapeKappa1);
		}else if (inString.find("scaleKappa1") == 0) {
			size_t pos = inString.find("=") + 1;
			string tmpStr = inString.substr(pos, inString.size() - pos);
			double scaleKappa1 = (double)atof(tmpStr.c_str());
			hyperParams.scaleKappa1(scaleKappa1);
		}else if(inString.find("aPhi")==0){
			size_t pos = inString.find("=")+1;
			string tmpStr = inString.substr(pos,inString.size()-pos);
			vector<double> aVec;
			while(tmpStr.find(" ")!=string::npos){
				pos = tmpStr.find(" ");
				if(pos==(tmpStr.size()-1)){
					string elem = tmpStr.substr(0,pos);
					aVec.push_back((double)atof(elem.c_str()));
					tmpStr = tmpStr.substr(pos+1,tmpStr.size());
					break;
				}
				string elem = tmpStr.substr(0,pos);
				aVec.push_back((double)atof(elem.c_str()));
				tmpStr = tmpStr.substr(pos+1,tmpStr.size());
			}
			hyperParams.aPhi(aVec);
		}else if(inString.find("mu0")==0){
			size_t pos = inString.find("=")+1;
			string tmpStr = inString.substr(pos,inString.size()-pos);
			vector<double> muVec;
			while(tmpStr.find(" ")!=string::npos){
				pos = tmpStr.find(" ");
				if(pos==(tmpStr.size()-1)){
					string elem = tmpStr.substr(0,pos);
					muVec.push_back((double)atof(elem.c_str()));
					tmpStr = tmpStr.substr(pos+1,tmpStr.size());
					break;
				}
				string elem = tmpStr.substr(0,pos);
				muVec.push_back((double)atof(elem.c_str()));
				tmpStr = tmpStr.substr(pos+1,tmpStr.size());
			}
			VectorXd mu0=VectorXd::Zero(muVec.size());
			for(unsigned int j=0;j<muVec.size();j++){
				mu0(j)=muVec[j];
			}
			hyperParams.mu0(mu0);
		}else if(inString.find("Tau0")==0){
			size_t pos = inString.find("=")+1;
			string tmpStr = inString.substr(pos,inString.size()-pos);
			vector<double> TauVec;
			while(tmpStr.find(" ")!=string::npos){
				pos = tmpStr.find(" ");
				if(pos==(tmpStr.size()-1)){
					string elem = tmpStr.substr(0,pos);
					TauVec.push_back((double)atof(elem.c_str()));
					tmpStr = tmpStr.substr(pos+1,tmpStr.size());
					break;
				}
				string elem = tmpStr.substr(0,pos);
				TauVec.push_back((double)atof(elem.c_str()));
				tmpStr = tmpStr.substr(pos+1,tmpStr.size()-pos-1);
			}
			unsigned int dim = (unsigned int)sqrt((double)TauVec.size());
			MatrixXd Tau0=MatrixXd::Zero(dim,dim);
			for(unsigned int j1=0;j1<dim;j1++){
				for(unsigned int j2=0;j2<dim;j2++){
					Tau0(j1,j2)=TauVec[j1*dim+j2];
				}
			}
			hyperParams.Tau0(Tau0);
		}else if (inString.find("Tau_Indep_0") == 0) {
			size_t pos = inString.find("=") + 1;
			string tmpStr = inString.substr(pos, inString.size() - pos);
			vector<double> TauVec_Indep;
			while (tmpStr.find(" ") != string::npos) {
				pos = tmpStr.find(" ");
				if (pos == (tmpStr.size() - 1)) {
					string elem = tmpStr.substr(0, pos);
					TauVec_Indep.push_back((double)atof(elem.c_str()));
					tmpStr = tmpStr.substr(pos + 1, tmpStr.size());
					break;
				}
				string elem = tmpStr.substr(0, pos);
				TauVec_Indep.push_back((double)atof(elem.c_str()));
				tmpStr = tmpStr.substr(pos + 1, tmpStr.size());
			}
			VectorXd Tau0_Indep = VectorXd::Zero(TauVec_Indep.size());
			for (unsigned int j = 0; j<TauVec_Indep.size(); j++) {
				Tau0_Indep(j) = TauVec_Indep[j];
			}
			hyperParams.Tau0_Indep(Tau0_Indep);

		}else if(inString.find("R0")==0){
			size_t pos = inString.find("=")+1;
			string tmpStr = inString.substr(pos,inString.size()-pos);
			vector<double> RVec;
			while(tmpStr.find(" ")!=string::npos){
				pos = tmpStr.find(" ");
				if(pos==(tmpStr.size()-1)){
					string elem = tmpStr.substr(0,pos);
					RVec.push_back((double)atof(elem.c_str()));
					tmpStr = tmpStr.substr(pos+1,tmpStr.size());
					break;
				}
				string elem = tmpStr.substr(0,pos);
				RVec.push_back((double)atof(elem.c_str()));
				tmpStr = tmpStr.substr(pos+1,tmpStr.size()-pos-1);
			}
			unsigned int dim = (unsigned int)sqrt((double)RVec.size());
			MatrixXd R0=MatrixXd::Zero(dim,dim);
			for(unsigned int j1=0;j1<dim;j1++){
				for(unsigned int j2=0;j2<dim;j2++){
					R0(j1,j2)=RVec[j1*dim+j2];
				}
			}
			hyperParams.R0(R0);
		}else if (inString.find("R_Indep_0") == 0) {
			size_t pos = inString.find("=") + 1;
			string tmpStr = inString.substr(pos, inString.size() - pos);
			vector<double> RVec_Indep;
			while (tmpStr.find(" ") != string::npos) {
				pos = tmpStr.find(" ");
				if (pos == (tmpStr.size() - 1)) {
					string elem = tmpStr.substr(0, pos);
					RVec_Indep.push_back((double)atof(elem.c_str()));
					tmpStr = tmpStr.substr(pos + 1, tmpStr.size());
					break;
				}
				string elem = tmpStr.substr(0, pos);
				RVec_Indep.push_back((double)atof(elem.c_str()));
				tmpStr = tmpStr.substr(pos + 1, tmpStr.size());
			}
			VectorXd R0_Indep = VectorXd::Zero(RVec_Indep.size());
			for (unsigned int j = 0; j<RVec_Indep.size(); j++) {
				R0_Indep(j) = RVec_Indep[j];
			}
			hyperParams.R0_Indep(R0_Indep);

		}else if(inString.find("kappa0")==0){
			size_t pos = inString.find("=")+1;
			string tmpStr = inString.substr(pos,inString.size()-pos);
			double kappa0 = (double)atof(tmpStr.c_str());
			hyperParams.kappa0(kappa0);
		}
		else if (inString.find("kappa1") == 0) {
			size_t pos = inString.find("=") + 1;
			string tmpStr = inString.substr(pos, inString.size() - pos);
			double kappa1 = (double)atof(tmpStr.c_str());
			hyperParams.kappa1(kappa1);
		}else if(inString.find("nu0")==0){
			size_t pos = inString.find("=")+1;
			string tmpStr = inString.substr(pos,inString.size()-pos);
			double nu0 = (double)atof(tmpStr.c_str());
			hyperParams.nu0(nu0);
		}else if(inString.find("muTheta")==0){
			size_t pos = inString.find("=")+1;
			string tmpStr = inString.substr(pos,inString.size()-pos);
			double muTheta = (double)atof(tmpStr.c_str());
			hyperParams.muTheta(muTheta);
		}else if(inString.find("sigmaTheta")==0){
			size_t pos = inString.find("=")+1;
			string tmpStr = inString.substr(pos,inString.size()-pos);
			double sigmaTheta = (double)atof(tmpStr.c_str());
			hyperParams.sigmaTheta(sigmaTheta);
		}else if(inString.find("dofTheta")==0){
			size_t pos = inString.find("=")+1;
			string tmpStr = inString.substr(pos,inString.size()-pos);
			unsigned int dofTheta = (unsigned int)atoi(tmpStr.c_str());
			hyperParams.dofTheta(dofTheta);
		}else if(inString.find("muBeta")==0){
			size_t pos = inString.find("=")+1;
			string tmpStr = inString.substr(pos,inString.size()-pos);
			double muBeta = (double)atof(tmpStr.c_str());
			hyperParams.muBeta(muBeta);
		}else if(inString.find("sigmaBeta")==0){
			size_t pos = inString.find("=")+1;
			string tmpStr = inString.substr(pos,inString.size()-pos);
			double sigmaBeta = (double)atof(tmpStr.c_str());
			hyperParams.sigmaBeta(sigmaBeta);
		}else if(inString.find("dofBeta")==0){
			size_t pos = inString.find("=")+1;
			string tmpStr = inString.substr(pos,inString.size()-pos);
			unsigned int dofBeta = (unsigned int)atoi(tmpStr.c_str());
			hyperParams.dofBeta(dofBeta);
		}else if(inString.find("shapeTauEpsilon")==0){
			size_t pos = inString.find("=")+1;
			string tmpStr = inString.substr(pos,inString.size()-pos);
			double shapeTauEpsilon = (double)atof(tmpStr.c_str());
			hyperParams.shapeTauEpsilon(shapeTauEpsilon);
		}else if(inString.find("rateTauEpsilon")==0){
			size_t pos = inString.find("=")+1;
			string tmpStr = inString.substr(pos,inString.size()-pos);
			double rateTauEpsilon = (double)atof(tmpStr.c_str());
			hyperParams.rateTauEpsilon(rateTauEpsilon);
		}else if(inString.find("aRho")==0){
			size_t pos = inString.find("=")+1;
			string tmpStr = inString.substr(pos,inString.size()-pos);
			double aRho = (double)atof(tmpStr.c_str());
			hyperParams.aRho(aRho);
		}else if(inString.find("bRho")==0){
			size_t pos = inString.find("=")+1;
			string tmpStr = inString.substr(pos,inString.size()-pos);
			double bRho = (double)atof(tmpStr.c_str());
			hyperParams.bRho(bRho);
		}else if(inString.find("atomRho")==0){
			size_t pos = inString.find("=")+1;
			string tmpStr = inString.substr(pos,inString.size()-pos);
			double atomRho = (double)atof(tmpStr.c_str());
			hyperParams.atomRho(atomRho);
			if(hyperParams.atomRho()<=0 || hyperParams.atomRho()>1){
				// Illegal atomRho value entered - it must be in (0,1] where 1 corresponds to the non-sparsity inducing var selection
				wasError=true;
				break;
			}
		}else if(inString.find("shapeSigmaSqY")==0){
			size_t pos = inString.find("=")+1;
			string tmpStr = inString.substr(pos,inString.size()-pos);
			double shapeSigmaSqY = (double)atof(tmpStr.c_str());
			hyperParams.shapeSigmaSqY(shapeSigmaSqY);
		}else if(inString.find("scaleSigmaSqY")==0){
			size_t pos = inString.find("=")+1;
			string tmpStr = inString.substr(pos,inString.size()-pos);
			double scaleSigmaSqY = (double)atof(tmpStr.c_str());
			hyperParams.scaleSigmaSqY(scaleSigmaSqY);
		}else if(inString.find("pQuantile")==0){
			size_t pos = inString.find("=")+1;
			string tmpStr = inString.substr(pos,inString.size()-pos);
			double pQuantile = (double)atof(tmpStr.c_str());
			hyperParams.pQuantile(pQuantile);
		}else if(inString.find("shapeNu")==0){
			size_t pos = inString.find("=")+1;
			string tmpStr = inString.substr(pos,inString.size()-pos);
			double shapeNu = (double)atof(tmpStr.c_str());
			hyperParams.shapeNu(shapeNu);
		}else if(inString.find("scaleNu")==0){
			size_t pos = inString.find("=")+1;
			string tmpStr = inString.substr(pos,inString.size()-pos);
			double scaleNu = (double)atof(tmpStr.c_str());
			hyperParams.scaleNu(scaleNu);
		}else if(inString.find("rSlice")==0){
			size_t pos = inString.find("=")+1;
			string tmpStr = inString.substr(pos,inString.size()-pos);
			double rSlice = (double)atof(tmpStr.c_str());
			hyperParams.rSlice(rSlice);
		}else if(inString.find("truncationEps")==0){
			size_t pos = inString.find("=")+1;
			string tmpStr = inString.substr(pos,inString.size()-pos);
			double truncationEps = (double)atof(tmpStr.c_str());
			hyperParams.truncationEps(truncationEps);
		}else if(inString.find("shapeTauCAR")==0){
			size_t pos = inString.find("=")+1;
			string tmpStr = inString.substr(pos,inString.size()-pos);
			double shapeTauCAR = (double)atof(tmpStr.c_str());
			hyperParams.shapeTauCAR(shapeTauCAR);
		}else if(inString.find("rateTauCAR")==0){
			size_t pos = inString.find("=")+1;
			string tmpStr = inString.substr(pos,inString.size()-pos);
			double rateTauCAR = (double)atof(tmpStr.c_str());
			hyperParams.rateTauCAR(rateTauCAR);
		}else if(inString.find("initAlloc")==0){
			size_t pos = inString.find("=")+1;
			string tmpStr = inString.substr(pos,inString.size()-pos);
			vector<double> initAl;
			while(tmpStr.find(" ")!=string::npos){
				pos = tmpStr.find(" ");
				if(pos==(tmpStr.size()-1)){
					string elem = tmpStr.substr(0,pos);
					initAl.push_back((double)atof(elem.c_str()));
					tmpStr = tmpStr.substr(pos+1,tmpStr.size());
					break;
				}
				string elem = tmpStr.substr(0,pos);
				initAl.push_back((double)atof(elem.c_str()));
				tmpStr = tmpStr.substr(pos+1,tmpStr.size());
			}
			hyperParams.initAlloc(initAl);
		}

	}

	// Return if there was an error
	if(wasError){
		Rprintf("There is a mistake in the arguments provided in profRegr.\n");
		Rprintf("The code will be run with default values.\n");
	//	Rprintf("Please use:\n");
	//	Rprintf("\t profileRegression --help\n");
	//	Rprintf("to get help on correct usage.\n");
	//	exit(-1);
	}
}

// Initialise the PReMiuM object (needed in this file as it calls
// function to read hyper parameters)
void initialisePReMiuM(baseGeneratorType& rndGenerator,
				const mcmcModel<pReMiuMParams,pReMiuMOptions,pReMiuMData>& model,
				pReMiuMParams& params){

	const pReMiuMData& dataset = model.dataset();
	const pReMiuMOptions& options = model.options();
	pReMiuMHyperParams& hyperParams = params.hyperParams();

	unsigned int nSubjects=dataset.nSubjects();
	unsigned int nCovariates=dataset.nCovariates();
	unsigned int nDiscreteCovs=dataset.nDiscreteCovs();
	unsigned int nContinuousCovs=dataset.nContinuousCovs();
	unsigned int nFixedEffects=dataset.nFixedEffects();
	unsigned int nCategoriesY=dataset.nCategoriesY();
	unsigned int nPredictSubjects=dataset.nPredictSubjects();
	unsigned int nClusInit = options.nClusInit();
	string covariateType = options.covariateType();
	string outcomeType = options.outcomeType();
	string hyperParamFileName = options.hyperParamFileName();
	string varSelectType = options.varSelectType();
	string samplerType = options.samplerType();
	bool includeResponse = options.includeResponse();
	bool responseExtraVar = options.responseExtraVar();
	bool includeCAR=options.includeCAR();
	bool includeuCARinit=options.includeuCARinit();
	string predictType = options.predictType();
	bool weibullFixedShape = options.weibullFixedShape();
	string uCARinitFileName = options.uCARinitFileName();
	bool useHyperpriorR1 = options.useHyperpriorR1();
	bool useIndependentNormal = options.useIndependentNormal();
	bool useSeparationPrior = options.useSeparationPrior();

	bool wasError=false;

	vector<unsigned int> nCategories;
	nCategories = dataset.nCategories();

	// Set the hyper parameters to their default values
	hyperParams.setSizes(nCovariates,nDiscreteCovs,
			nContinuousCovs,covariateType,useIndependentNormal, useHyperpriorR1, useSeparationPrior);
	hyperParams.setDefaults(dataset,options);
	// Read the parameters from file if file provided
	if(hyperParamFileName.compare("")!=0){
		readHyperParamsFromFile(hyperParamFileName,hyperParams);
	}

	// Allocate the right sizes for each of the parameter variables
	// This also switches "on" all variable indicators (gamma)
	// This gets changed below if variable selection is being done
	params.setSizes(nSubjects,nCovariates,nDiscreteCovs,nContinuousCovs,nFixedEffects,nCategoriesY,
		nPredictSubjects,nCategories,nClusInit,covariateType,weibullFixedShape,useHyperpriorR1,useIndependentNormal, useSeparationPrior);
	unsigned int maxNClusters=params.maxNClusters();


	// Define a uniform random number generator
	randomUniform unifRand(0,1);

	if(nClusInit==0){
		nClusInit=50+(unsigned int)11*unifRand(rndGenerator);
	}

	// Fix the number of clusters if we are using the truncated sampler
	if(samplerType.compare("Truncated")==0){
		maxNClusters=20;
		if((nClusInit+10)>maxNClusters){
			maxNClusters=nClusInit+10;
		}
		// Now compute the bound recommended in Ishwaran and James 2001
		double multiplier=0.0;
		if(options.fixedAlpha()>-1){
			multiplier=options.fixedAlpha();
		}else{
			// Use the expected value of alpha as the multiplier
			multiplier=hyperParams.shapeAlpha()/hyperParams.rateAlpha();
		}
		double computedBound=1+multiplier*(log(4.0*nSubjects)-log(hyperParams.truncationEps()));
		if(computedBound>maxNClusters){
			maxNClusters=computedBound;
		}
		params.maxNClusters(maxNClusters,covariateType, useIndependentNormal, useSeparationPrior);
	}

	// Copy the dataset X matrix to a working object in params
	params.workDiscreteX(dataset.discreteX());
	params.workContinuousX(dataset.continuousX());

	// Now initialise the actual parameters
	randomGamma gammaRand(hyperParams.shapeAlpha(),1.0/hyperParams.rateAlpha());

	double alpha=gammaRand(rndGenerator);
	if(options.fixedAlpha()>-1){
		alpha=options.fixedAlpha();
	}
	params.alpha(alpha);	

	double dPitmanYor = options.dPitmanYor();
	params.dPitmanYor(dPitmanYor);	

	vector<unsigned int> nXInCluster(maxNClusters,0);
	unsigned int maxZ=0;
	params.workNClusInit(nClusInit);
	if (hyperParams.initAlloc().empty()){
		for(unsigned int i=0;i<nSubjects+nPredictSubjects;i++){
			int c=(int) nClusInit*unifRand(rndGenerator);
			params.z(i,c,covariateType,useIndependentNormal);
			if(c>(int)maxZ){
				maxZ=c;
			}
			if(i<nSubjects){
				nXInCluster[c]++;
			}
		}
	} else {
		for(unsigned int i=0;i<nSubjects+nPredictSubjects;i++){
			int c = hyperParams.initAlloc(i);
			params.z(i,c,covariateType,useIndependentNormal);
			if(c>(int)maxZ){
				maxZ=c;
			}
			if(i<nSubjects){
				nXInCluster[c]++;
			}
		}
	}

	params.workNXInCluster(nXInCluster);
	params.workMaxZi(maxZ);


	// Sample v (for logPsi)
	// This is sampled from the posterior given the z vector above
	// Prior comes from the conjugacy of the dirichlet and multinomial
	// See Ishwaran and James 2001

	// Sample active V
	vector<unsigned int> sumCPlus1ToMaxMembers(maxZ+1,0);
	for(int c=maxZ-1;c>=0;c--){
		sumCPlus1ToMaxMembers[c]=sumCPlus1ToMaxMembers[c+1]+params.workNXInCluster(c+1);
	}

	double tmp=0.0;

	for(unsigned int c=0;c<=maxZ;c++){
		double vVal = betaRand(rndGenerator,1.0+params.workNXInCluster(c)-dPitmanYor,alpha+sumCPlus1ToMaxMembers[c]+dPitmanYor*(c+1));
		params.v(c,vVal);
		// Set logPsi
		params.logPsi(c,tmp+log(vVal));
		tmp += log(1-vVal);
	}

	if(samplerType.compare("Truncated")==0){
		// Just sample the remaining V from the prior
		vector<double> vNew=params.v();
		vector<double> logPsiNew=params.logPsi();

		for(unsigned int c=maxZ+1;c<maxNClusters;c++){
			double v=betaRand(rndGenerator,1.0-dPitmanYor,alpha+dPitmanYor*c);
			double logPsi=log(v)+log(1-vNew[c-1])-log(vNew[c-1])+logPsiNew[c-1];
			if(c>=vNew.size()){
				vNew.push_back(v);
				logPsiNew.push_back(logPsi);
			}else{
				vNew[c]=v;
				logPsiNew[c]=logPsi;
			}
		}
		params.v(vNew);
		params.logPsi(logPsiNew);

	}else{

		// Sample u (auxilliary variables). This will determine the maximum number of clusters
		double minU=1.0;
		for(unsigned int i=0;i<nSubjects+nPredictSubjects;i++){
			int zi=params.z(i);
			double ui=0.0;
			if(samplerType.compare("SliceDependent")==0){
				ui = exp(params.logPsi(zi))*unifRand(rndGenerator);
			}else if(samplerType.compare("SliceIndependent")==0){
				ui = hyperParams.workXiSlice(zi)*unifRand(rndGenerator);
			}
			if(ui<minU){
				minU=ui;
			}
			params.u(i,ui);
		}
		params.workMinUi(minU);

		// Sample V
		vector<double> cumPsi(maxZ+1,0.0);
		cumPsi[0] = exp(params.logPsi(0));
		for(unsigned int c=1;c<=maxZ;c++){
			cumPsi[c]=cumPsi[c-1]+exp(params.logPsi(c));
		}

		vector<double> vNew=params.v();
		vector<double> logPsiNew=params.logPsi();

		maxNClusters = maxZ+1;
		if(samplerType.compare("SliceIndependent")==0){

			maxNClusters=2+(int)((log(params.workMinUi())-log(1.0-hyperParams.rSlice()))/log(hyperParams.rSlice()));
		}

		bool continueLoop=true;
		unsigned int c=maxZ;
		while(continueLoop){
			if(samplerType.compare("SliceDependent")==0&&cumPsi[c]>1-minU){
				// We can stop
				maxNClusters=c+1;
				continueLoop=false;
			}else if(samplerType.compare("SliceIndependent")==0&&c>=maxNClusters){
				continueLoop=false;
			}else{
				c++;
				// We need a new sampled value of v
				double v=betaRand(rndGenerator,1.0-dPitmanYor,alpha+dPitmanYor*c);
				double logPsi=log(v)+log(1-vNew[c-1])-log(vNew[c-1])+logPsiNew[c-1];
				if(c>=vNew.size()){
					vNew.push_back(v);
					logPsiNew.push_back(logPsi);
				}else{
					vNew[c]=v;
					logPsiNew[c]=logPsi;
				}
				cumPsi.push_back(cumPsi[c-1]+exp(logPsi));
			}
		}

		params.maxNClusters(maxNClusters,covariateType, useIndependentNormal, useSeparationPrior);
		params.v(vNew);
		params.logPsi(logPsiNew);
	}

	if(covariateType.compare("Discrete")==0){
		// Sample logPhi
		// Need to count the number of X[i][j]==p for each covariate and category of p
		vector<vector<unsigned int> > nXpMembers(nCovariates);
		for(unsigned int j=0;j<nCovariates;j++){
			nXpMembers[j].resize(nCategories[j]);
			for(unsigned int p=0;p<nCategories[j];p++){
				nXpMembers[j][p]=0;
				for(unsigned int i=0;i<nSubjects;i++){
					if(dataset.discreteX(i,j)==(int)p&&!dataset.missingX(i,j)){
						nXpMembers[j][p]++;
					}
				}
			}
		}

		// Now we can sample. We don't use the priors, but instead look at the number
		// of people in each category and do a dirichlet sample that takes account of that
		boost::math::normal_distribution<double> norm01(0.0,1.0);
		for(unsigned int c=0;c<maxNClusters;c++){
			for(unsigned int j=0;j<nCovariates;j++){
				vector<double> dirichParams(nCategories[j]);
				for(unsigned int p=0;p<nCategories[j];p++){
					dirichParams[p]=(double)nXpMembers[j][p]+hyperParams.aPhi(j);
				}
				vector<double> logDirichSample(nCategories[j]);
				vector<double> dirichSample(nCategories[j]);
				dirichSample=dirichletRand(rndGenerator,dirichParams);
				for(unsigned int p=0;p<nCategories[j];p++){
					logDirichSample[p]=log(dirichSample[p]);
				}
				params.logPhi(c,j,logDirichSample);
			}
		}
		// Initialise the null parameters for the variable selection case
		// In all cases, initialise it at the value it will be fixed at for
		// the continuous indicator case
		if(varSelectType.compare("None")!=0){
			for(unsigned int j=0;j<nCovariates;j++){
				double sumVec=0.0;
				vector<double> probVec(nCategories[j],0.0000001);
				for(unsigned int p=0;p<nCategories[j];p++){
					probVec[p]+=(double)(nXpMembers[j][p]);
					sumVec+=(double)nXpMembers[j][p];
				}
				vector<double> logProbVec(nCategories[j]);
				for(unsigned int p=0;p<nCategories[j];p++){
					logProbVec[p]=log(probVec[p]/sumVec);
				}
				params.logNullPhi(j,logProbVec);
			}
		}

	}else if(covariateType.compare("Normal")==0){

		// In the following it is useful to have the rows of X as
		// Eigen dynamic vectors
		vector<VectorXd> xi(nSubjects);
		for(unsigned int i=0;i<nSubjects;i++){
			xi[i].setZero(nCovariates);
			for(unsigned int j=0;j<nCovariates;j++){
				xi[i](j)=dataset.continuousX(i,j);
			}
		}


		// Now we can sample from the conditionals (using Sigma_c=Sigma_0 and
		// mu_c=mu_0 for all c) to get mu_c and Sigma_c for each cluster

		// First we sample mu_c for each cluster

		// We begin by computing the mean X for individuals in each cluster
		vector<VectorXd> meanX(maxNClusters);
		for(unsigned int c=0;c<maxNClusters;c++){
			meanX[c].setZero(nCovariates);
		}
		for(unsigned int i=0;i<nSubjects;i++){
			meanX[params.z(i)]=meanX[params.z(i)]+xi[i];
		}

		for(unsigned int c=0;c<maxNClusters;c++){
			// Having computed this we can calcuate the posterior mean
			// and posterior covariance for each mu_c
			if(params.workNXInCluster(c)>0){
				meanX[c]=meanX[c]/(double)params.workNXInCluster(c);
			}else{
				meanX[c].setZero(nCovariates);
			}

			if (useIndependentNormal) {
				VectorXd mu(nCovariates);
				VectorXd mu0(nCovariates);
				mu0 = hyperParams.mu0();
				VectorXd Tau0(nCovariates);
				Tau0 = hyperParams.Tau0_Indep();
				for (unsigned int j = 0; j < nCovariates; j++) {
					double denom = 0;
					denom = (double)params.workNXInCluster(c)*(1 / Tau0(j)) + (1 / Tau0(j));
					double meanNum = 0;
					meanNum = (double)params.workNXInCluster(c)*(1 / Tau0(j))*meanX[c](j) + (1 / Tau0(j))*mu0(j);
					double variance = (1 / Tau0(j))*(1 / Tau0(j)) / denom;
					double mean = meanNum / denom;

					//Sample from the posterior
					mu(j) = NormalRand(rndGenerator, mean, variance);
				}

				// We store our sample
				params.mu(c, mu, useIndependentNormal);

			}else if (useHyperpriorR1) {
				VectorXd mu00 = hyperParams.mu0();
				params.mu00(mu00);

				params.Tau00(hyperParams.Tau0());

				MatrixXd covMat(nCovariates, nCovariates);
				covMat = (params.Tau00() + params.workNXInCluster(c)*params.Tau00()).inverse();
				VectorXd meanVec(nCovariates);
				meanVec = params.Tau00()*params.mu00() + params.workNXInCluster(c)*params.Tau00()*meanX[c];
				meanVec = covMat*meanVec;

				VectorXd mu(nCovariates);
				// We sample from this posterior
				mu = multivarNormalRand(rndGenerator, meanVec, covMat);

				// We store our sample
				params.mu(c, mu, useIndependentNormal);

			}
			else if (useSeparationPrior) {
				VectorXd mu00 = hyperParams.mu0();
				params.mu00(mu00);

				MatrixXd covMat(nCovariates, nCovariates);
				covMat = (hyperParams.Tau00() + params.workNXInCluster(c)*hyperParams.Tau00()).inverse();
				VectorXd meanVec(nCovariates);
				meanVec = hyperParams.Tau00()*params.mu00() + params.workNXInCluster(c)*hyperParams.Tau00()*meanX[c];
				meanVec = covMat*meanVec;

				VectorXd mu(nCovariates);
				// We sample from this posterior
				mu = multivarNormalRand(rndGenerator, meanVec, covMat);

				// We store our sample
				params.mu(c, mu, useIndependentNormal);

			}
			else {
				MatrixXd covMat(nCovariates, nCovariates);
				covMat = (hyperParams.Tau0() + params.workNXInCluster(c)*hyperParams.Tau0()).inverse();
				VectorXd meanVec(nCovariates);
				meanVec = hyperParams.Tau0()*hyperParams.mu0() + params.workNXInCluster(c)*hyperParams.Tau0()*meanX[c];
				meanVec = covMat*meanVec;

				VectorXd mu(nCovariates);
				// We sample from this posterior
				mu = multivarNormalRand(rndGenerator, meanVec, covMat);

				// We store our sample
				params.mu(c, mu, useIndependentNormal);
			}

		}


		//initialise mu00 and Tau00
		/* if (useHyperpriorR1 || useSeparationPrior) {
			VectorXd SumMu;
			SumMu.setZero(nCovariates);
			unsigned int workNactive = 0;
			for (unsigned int c = 0; c <= maxZ; c++) {
                // SumMu += hyperParams.mu0();
				SumMu += params.mu(c);
				workNactive += 1;
			}

			MatrixXd sigmaMu00(nCovariates, nCovariates);
			VectorXd meanMu00(nCovariates);

			sigmaMu00 = (workNactive*hyperParams.Tau0() + hyperParams.Tau0()).inverse();
			meanMu00 = sigmaMu00*(hyperParams.Tau0()*SumMu + hyperParams.Tau0()*hyperParams.mu0());


			VectorXd mu00(nCovariates);
			mu00 = multivarNormalRand(rndGenerator, meanMu00, sigmaMu00);

			params.mu00(mu00);

			if (useHyperpriorR1) {
				params.Tau00(hyperParams.Tau0());
			}  
		}   */



		//initialise R1 and kappa1 when using multivariate normal likelihood
		if (options.useHyperpriorR1()){

			//initialise kappa1
			params.kappa11(nCovariates + 1);

			//initialise R1
			MatrixXd SumTau = MatrixXd::Zero(nCovariates,nCovariates);
			unsigned int workNactive=0;
			for (unsigned int c=0; c<maxZ+1;c++){
				SumTau += hyperParams.Tau0();
				workNactive += 1;
			}
			SumTau += hyperParams.R0();
			MatrixXd R0Star = SumTau.inverse();
		
			MatrixXd inverseR1 =  wishartRand(rndGenerator,R0Star,workNactive*params.kappa11()+hyperParams.kappa0());
			MatrixXd R1 = inverseR1.inverse();
			params.R1(R1);

		}

		//initialise R1 when using independent normal likelihood
		if (useIndependentNormal) {
			VectorXd R1(nCovariates);
			unsigned int workNactive = 0;
			VectorXd sumTau;
			sumTau.setZero(nCovariates);
			for (unsigned int c = 0; c < maxZ + 1; c++) {
				sumTau += hyperParams.Tau0_Indep();
				workNactive = workNactive + 1;
			}

			VectorXd R0_Indep = hyperParams.R0_Indep();
			for (unsigned int j = 0; j < nCovariates; j++) {
				double kappaNew = workNactive*hyperParams.kappa1() + hyperParams.kappa0();
				double rNew = sumTau(j) + R0_Indep(j);

				randomGamma gammaRand(kappaNew, 1.0 / rNew);
				R1(j) = gammaRand(rndGenerator);

			}
			params.R1_Indep(R1);
		}


		// Now we can sample Tau_c for each cluster
		if (useIndependentNormal) {

			vector<VectorXd> XiMinusMuSq(nSubjects);
			for (unsigned int i = 0; i < nSubjects; i++) {
				XiMinusMuSq[i].setZero(nCovariates);
				VectorXd mu = params.mu(params.z(i));
				for (unsigned int j = 0; j < nCovariates; j++) {
					XiMinusMuSq[i](j) = (xi[i](j) - mu(j))*(xi[i](j) - mu(j));
				}
			}

			vector<VectorXd> sumXiMinusMuSq(maxNClusters);
			for (unsigned int c = 0; c <maxNClusters; c++) {
				sumXiMinusMuSq[c].setZero(nCovariates);
			}

			for (unsigned int i = 0; i<nSubjects; i++) {
				sumXiMinusMuSq[params.z(i)] = sumXiMinusMuSq[params.z(i)] + XiMinusMuSq[i];
			}

			for (unsigned int c = 0; c < maxNClusters; c++) {
				VectorXd tau(nCovariates);
				int nXInC = params.workNXInCluster(c);
				for (unsigned int j = 0; j < nCovariates; j++) {
					double kappaNew = (double)nXInC / 2 + hyperParams.kappa1();
					double rNew = (sumXiMinusMuSq[c](j) + 2 * params.R1_Indep(j)) / 2;

					randomGamma gammaRand(kappaNew, 1.0 / rNew);
					tau(j) = gammaRand(rndGenerator);
				}
				params.Tau_Indep(c, tau);
			}

		}
		else if (useSeparationPrior) {

			params.kappa11(nCovariates + 1);

			//set TauR to be indentity to begin with
			MatrixXd TauR;
			TauR.setZero(nCovariates, nCovariates);

			for (unsigned int j = 0; j < nCovariates; j++) {
				TauR(j, j) = 1;
			}
			for (unsigned int c = 0; c<maxNClusters; c++) {
				params.TauR(c, TauR);
			}

			//compute muc for each c for the initial allocation
			vector<VectorXd> muc(maxNClusters);
			for (unsigned int c = 0; c<maxNClusters; c++) {
				muc[c].setZero(nCovariates);
			}

			for (unsigned int i = 0; i<nSubjects; i++) {
				unsigned int zi = params.z(i);
				muc[zi] = muc[zi] + xi[i];
			} 

			for (unsigned int c = 0; c < maxNClusters; c++) {
				if (params.workNXInCluster(c) > 0) {
					muc[c] = muc[c] / (double)params.workNXInCluster(c);
				}
				else {
					muc[c].setZero(nCovariates);
				}
			}

			vector<MatrixXd> Sigmas(maxNClusters);
			for (unsigned int c = 0; c < maxNClusters; c++) {
				Sigmas[c].setZero(nCovariates, nCovariates);
			}

			for (unsigned int i = 0; i<nSubjects; i++) {
				unsigned int zi = params.z(i);
				Sigmas[zi] += (xi[i] - muc[zi])*((xi[i] - muc[zi]).transpose());
			}


			for (unsigned int c = 0; c < maxNClusters; c++) {
				if (params.workNXInCluster(c) > 1) {
					Sigmas[c] = Sigmas[c] / (double)params.workNXInCluster(c);
				}
				else {
					Sigmas[c].setZero(nCovariates, nCovariates);
					VectorXd beta_taus = hyperParams.beta_taus0();
					double alpha_taus = hyperParams.alpha_taus();
					for (unsigned int j = 0; j < nCovariates; j++) {
						randomGamma gammaRand(alpha_taus, 1.0 / beta_taus(j));
						double tausj = gammaRand(rndGenerator);
						Sigmas[c](j, j) = (1.0 / tausj)*(1.0 / tausj);
					}					
				}			
			}

	
			for (unsigned int c = 0; c<maxNClusters; c++) {
				for (unsigned int j = 0; j < nCovariates; j++) {
					double Tausj = 1.0 / sqrt(Sigmas[c](j, j));
					params.TauS(c, j, Tausj);
				}
			}

		}
		else {

			vector<MatrixXd> Rc(maxNClusters);
			for (unsigned int c = 0; c<maxNClusters; c++) {
				Rc[c].setZero(nCovariates, nCovariates);
			}

			for (unsigned int i = 0; i<nSubjects; i++) {
				unsigned int zi = params.z(i);
				Rc[zi] = Rc[zi] + (xi[i] - params.mu(zi))*((xi[i] - params.mu(zi)).transpose());
			}

			if (options.useHyperpriorR1()) {
				for (unsigned int c = 0; c<maxNClusters; c++) {
					Rc[c] = (params.R1().inverse() + Rc[c]).inverse();
					MatrixXd Tau = wishartRand(rndGenerator, Rc[c], params.workNXInCluster(c) + params.kappa11());
					params.Tau(c, Tau);
				}
			}
			else {
				for (unsigned int c = 0; c<maxNClusters; c++) {
					Rc[c] = (hyperParams.R0().inverse() + Rc[c]).inverse();
					MatrixXd Tau = wishartRand(rndGenerator, Rc[c], params.workNXInCluster(c) + hyperParams.kappa0());
					params.Tau(c, Tau);
				}

			}

		}

		if (useSeparationPrior) {
			VectorXd sumTau = VectorXd::Zero(nCovariates);
			unsigned int workNactive = 0;
			MatrixXd Tau0 = hyperParams.Tau0();
			for (unsigned int c = 0; c < maxZ + 1; c++) {
				for (unsigned int j = 0; j < nCovariates; j++) {
					sumTau(j) = sumTau(j) + params.TauS(c, j);
				}
				workNactive += 1;
			}

			VectorXd betas_new(nCovariates);
			VectorXd beta_taus0 = hyperParams.beta_taus0();
			for (unsigned int j = 0; j < nCovariates; j++) {
				double alpha_new = workNactive * hyperParams.alpha_taus() + hyperParams.alpha_taus0();
				double beta_new = sumTau(j) + beta_taus0(j);

				randomGamma gammaRand(alpha_new, 1.0 / beta_new);
				betas_new(j) = gammaRand(rndGenerator);

			}
			params.beta_taus(betas_new);

		}
		

		// Now do the null mu for variable selection
		// In all cases, initialise it at the value it will be fixed at for
		// the continuous indicator case
		if(varSelectType.compare("None")!=0){
			vector<double> meanXVec(nCovariates,0.0);
			vector<unsigned int> countXVec(nCovariates,0);
			for(unsigned int i=0;i<nSubjects;i++){
				for(unsigned int j=0;j<nCovariates;j++){
					if(!dataset.missingX(i,j)){
						meanXVec[j]+=dataset.continuousX(i,j);
						countXVec[j]+=1;
					}
				}
			}
			VectorXd nullMu=VectorXd::Zero(nCovariates);
			for(unsigned int j=0;j<nCovariates;j++){
				nullMu(j)=meanXVec[j]/(double)countXVec[j];
			}
			params.nullMu(nullMu,useIndependentNormal);
		}

	}else if(covariateType.compare("Mixed")==0){
		// Sample logPhi
		// Need to count the number of X[i][j]==p for each covariate and category of p
		vector<vector<unsigned int> > nXpMembers(nDiscreteCovs);
		for(unsigned int j=0;j<nDiscreteCovs;j++){
			nXpMembers[j].resize(nCategories[j]);
			for(unsigned int p=0;p<nCategories[j];p++){
				nXpMembers[j][p]=0;
				for(unsigned int i=0;i<nSubjects;i++){
					if(dataset.discreteX(i,j)==(int)p&&!dataset.missingX(i,j)){
						nXpMembers[j][p]++;
					}
				}
			}
		}

		// Now we can sample. We don't use the priors, but instead look at the number
		// of people in each category and do a dirichlet sample that takes account of that
		boost::math::normal_distribution<double> norm01(0.0,1.0);
		for(unsigned int c=0;c<maxNClusters;c++){
			for(unsigned int j=0;j<nDiscreteCovs;j++){
				vector<double> dirichParams(nCategories[j]);
				for(unsigned int p=0;p<nCategories[j];p++){
					dirichParams[p]=(double)nXpMembers[j][p]+hyperParams.aPhi(j);
				}
				vector<double> logDirichSample(nCategories[j]);
				vector<double> dirichSample(nCategories[j]);
				dirichSample=dirichletRand(rndGenerator,dirichParams);
				for(unsigned int p=0;p<nCategories[j];p++){
					logDirichSample[p]=log(dirichSample[p]);
				}
				params.logPhi(c,j,logDirichSample);
			}
		}
		// Initialise the null parameters for the variable selection case
		// In all cases, initialise it at the value it will be fixed at for
		// the continuous indicator case
		if(varSelectType.compare("None")!=0){
			for(unsigned int j=0;j<nDiscreteCovs;j++){
				double sumVec=0.0;
				vector<double> probVec(nCategories[j],0.0000001);
				for(unsigned int p=0;p<nCategories[j];p++){
					probVec[p]+=(double)(nXpMembers[j][p]);
					sumVec+=(double)nXpMembers[j][p];
				}
				vector<double> logProbVec(nCategories[j]);
				for(unsigned int p=0;p<nCategories[j];p++){
					logProbVec[p]=log(probVec[p]/sumVec);
				}
				params.logNullPhi(j,logProbVec);
			}
		}

		// In the following it is useful to have the rows of X as
		// Eigen dynamic vectors
		vector<VectorXd> xi(nSubjects);
		for(unsigned int i=0;i<nSubjects;i++){
			xi[i].setZero(nContinuousCovs);
			for(unsigned int j=0;j<nContinuousCovs;j++){
				xi[i](j)=dataset.continuousX(i,j);
			}
		}

		//initialise mu00
		/*if (useHyperpriorR1) {
			VectorXd SumMu;
			SumMu.setZero(nContinuousCovs);
			unsigned int workNactive = 0;
			for (unsigned int c = 0; c <= maxZ; c++) {
				SumMu += hyperParams.mu0();
				workNactive += 1;
			}

			MatrixXd sigmaMu00(nContinuousCovs, nContinuousCovs);
			VectorXd meanMu00(nContinuousCovs);

			sigmaMu00 = (workNactive*hyperParams.Tau0() + hyperParams.Tau0()).inverse();
			meanMu00 = sigmaMu00*(hyperParams.Tau0()*SumMu + hyperParams.Tau0()*hyperParams.mu0());

			VectorXd mu00(nContinuousCovs);
			mu00 = multivarNormalRand(rndGenerator, meanMu00, sigmaMu00);

			params.mu00(mu00);

			VectorXd SumMuDiff;
			SumMuDiff.setZero(nContinuousCovs);
			for (unsigned int c = 0; c <= maxZ; c++) {
				SumMuDiff = (params.mu00() - hyperParams.mu0())*((params.mu00() - hyperParams.mu0()).transpose());
				workNactive += 1;
			}

			SumMuDiff += hyperParams.R00().inverse();

			MatrixXd RUpadated(nContinuousCovs, nContinuousCovs);
			RUpadated = SumMuDiff.inverse();

			MatrixXd Tau00(nContinuousCovs, nContinuousCovs);
			Tau00 = wishartRand(rndGenerator, RUpadated, workNactive + hyperParams.kappa00());

			params.Tau00(Tau00);

		} */

		// Now we can sample from the conditionals (using Sigma_c=Sigma_0 and
		// mu_c=mu_0 for all c) to get mu_c and Sigma_c for each cluster

		// First we sample mu_c for each cluster

		// We begin by computing the mean X for individuals in each cluster
		vector<VectorXd> meanX(maxNClusters);
		for(unsigned int c=0;c<maxNClusters;c++){
			meanX[c].setZero(nContinuousCovs);
		}
		for(unsigned int i=0;i<nSubjects;i++){
			meanX[params.z(i)]=meanX[params.z(i)]+xi[i];
		}

		for(unsigned int c=0;c<maxNClusters;c++){
			// Having computed this we can calcuate the posterior mean
			// and posterior covariance for each mu_c
			if(params.workNXInCluster(c)>0){
				meanX[c]=meanX[c]/(double)params.workNXInCluster(c);
			}else{
				meanX[c].setZero(nContinuousCovs);
			}

			if (useIndependentNormal) {
				VectorXd mu(nContinuousCovs);
				VectorXd mu0(nContinuousCovs);
				mu0 = hyperParams.mu0();
				VectorXd Tau0(nContinuousCovs);
				Tau0 = hyperParams.Tau0_Indep();
				for (unsigned int j = 0; j < nContinuousCovs; j++) {
					double denom = 0;
					denom = (double)params.workNXInCluster(c)*(1 / Tau0(j)) + (1 / Tau0(j));
					double meanNum = 0;
					meanNum = (double)params.workNXInCluster(c)*(1 / Tau0(j))*meanX[c](j) + (1 / Tau0(j))*mu0(j);
					double variance = (1 / Tau0(j))*(1 / Tau0(j)) / denom;
					double mean = meanNum / denom;

					//Sample from the posterior
					mu(j) = NormalRand(rndGenerator, mean, variance);
				}
				// We store our sample
				params.mu(c, mu, useIndependentNormal);
			}
			else if (useHyperpriorR1) {

				VectorXd mu00 = hyperParams.mu0();
				params.mu00(mu00);

				params.Tau00(hyperParams.Tau0());

				MatrixXd covMat(nContinuousCovs, nContinuousCovs);
				covMat = (params.Tau00() + params.workNXInCluster(c)*params.Tau00()).inverse();
				VectorXd meanVec(nContinuousCovs);
				meanVec = params.Tau00()*params.mu00() + params.workNXInCluster(c)*params.Tau00()*meanX[c];
				meanVec = covMat*meanVec;

				VectorXd mu(nContinuousCovs);
				// We sample from this posterior
				mu = multivarNormalRand(rndGenerator, meanVec, covMat);

				// We store our sample
				params.mu(c, mu, useIndependentNormal);

			}else if (useSeparationPrior) {
				VectorXd mu00 = hyperParams.mu0();
				params.mu00(mu00);

				MatrixXd covMat(nContinuousCovs, nContinuousCovs);
				covMat = (hyperParams.Tau00() + params.workNXInCluster(c)*hyperParams.Tau00()).inverse();
				VectorXd meanVec(nContinuousCovs);
				meanVec = hyperParams.Tau00()*params.mu00() + params.workNXInCluster(c)*hyperParams.Tau00()*meanX[c];
				meanVec = covMat*meanVec;

				VectorXd mu(nContinuousCovs);
				// We sample from this posterior
				mu = multivarNormalRand(rndGenerator, meanVec, covMat);

				// We store our sample
				params.mu(c, mu, useIndependentNormal);

			}else {

				MatrixXd covMat(nContinuousCovs, nContinuousCovs);
				covMat = (hyperParams.Tau0() + params.workNXInCluster(c)*hyperParams.Tau0()).inverse();
				VectorXd meanVec(nContinuousCovs);
				meanVec = hyperParams.Tau0()*hyperParams.mu0() + params.workNXInCluster(c)*hyperParams.Tau0()*meanX[c];
				meanVec = covMat*meanVec;

				VectorXd mu(nContinuousCovs);
				// We sample from this posterior
				mu = multivarNormalRand(rndGenerator, meanVec, covMat);

				// We store our sample
				params.mu(c, mu, useIndependentNormal);
			}

		}

		//initialise R1 and kappa1 when using multivariate normal likelihood
		if (options.useHyperpriorR1()){

			params.kappa11(nContinuousCovs + 1);

			MatrixXd SumTau = MatrixXd::Zero(nContinuousCovs,nContinuousCovs);
			unsigned int workNactive=0;
			for (unsigned int c=0; c<maxZ+1;c++){
				SumTau += hyperParams.Tau0();
				workNactive += 1;
			}
			SumTau += hyperParams.R0();
			MatrixXd R0Star = SumTau.inverse();
		
			MatrixXd inverseR1 =  wishartRand(rndGenerator,R0Star,workNactive*params.kappa11()+hyperParams.kappa0());
			MatrixXd R1 = inverseR1.inverse();
			params.R1(R1);

		}

		//initialise R1 when using independent normal likelihood
		if (useIndependentNormal) {
			VectorXd R1(nContinuousCovs);
			unsigned int workNactive = 0;
			VectorXd sumTau;
			sumTau.setZero(nContinuousCovs);
			for (unsigned int c = 0; c < maxZ + 1; c++) {
				sumTau = hyperParams.Tau0_Indep();
				workNactive = workNactive + 1;
			}

			VectorXd R0_Indep = hyperParams.R0_Indep();
			for (unsigned int j = 0; j < nContinuousCovs; j++) {
				double kappaNew = workNactive*hyperParams.kappa1() + hyperParams.kappa0();
				double rNew = sumTau(j) + R0_Indep(j);

				randomGamma gammaRand(kappaNew, 1.0 / rNew);
				R1(j) = gammaRand(rndGenerator);
			}
			params.R1_Indep(R1);
		}


		// Now we can sample Tau_c for each cluster
		if (useIndependentNormal) {

			vector<VectorXd> XiMinusMuSq(nSubjects);
			for (unsigned int i = 0; i < nSubjects; i++) {
				XiMinusMuSq[i].setZero(nContinuousCovs);
				unsigned int zi = params.z(i);
				VectorXd mu = params.mu(zi);
				for (unsigned int j = 0; j < nContinuousCovs; j++) {
					XiMinusMuSq[i](j) = (xi[i](j) - mu(j))*(xi[i](j) - mu(j));
				}
			}

			vector<VectorXd> sumXiMinusMuSq(maxZ + 1);
			for (unsigned int c = 0; c <= maxZ; c++) {
				sumXiMinusMuSq[c].setZero(nContinuousCovs);
			}

			for (unsigned int i = 0; i<nSubjects; i++) {
				sumXiMinusMuSq[params.z(i)] = sumXiMinusMuSq[params.z(i)] + XiMinusMuSq[i];
			}

			for (unsigned int c = 0; c <= maxZ; c++) {
				VectorXd tau(nContinuousCovs);
				int nXInC = params.workNXInCluster(c);
				for (unsigned int j = 0; j < nContinuousCovs; j++) {
					double kappaNew = (double)nXInC / 2 + hyperParams.kappa1();
					double rNew = (sumXiMinusMuSq[c](j) + 2 * params.R1_Indep(j)) / 2;

					randomGamma gammaRand(kappaNew, 1.0 / rNew);
					tau(j) = gammaRand(rndGenerator);
				}
				params.Tau_Indep(c, tau);
			}

		}
		else if (useSeparationPrior) {

			params.kappa11(nContinuousCovs + 1);

			//set TauR to be indentity to begin with
			MatrixXd TauR;
			TauR.setZero(nContinuousCovs, nContinuousCovs);

			for (unsigned int j = 0; j < nContinuousCovs; j++) {
				TauR(j, j) = 1;
			}
			for (unsigned int c = 0; c<maxNClusters; c++) {
				params.TauR(c, TauR);
			}

			//compute muc for each c for the initial allocation
			vector<VectorXd> muc(maxNClusters);
			for (unsigned int c = 0; c<maxNClusters; c++) {
				muc[c].setZero(nContinuousCovs);
			}

			for (unsigned int i = 0; i<nSubjects; i++) {
				unsigned int zi = params.z(i);
				muc[zi] = muc[zi] + xi[i];
			}

			for (unsigned int c = 0; c < maxNClusters; c++) {
				if (params.workNXInCluster(c) > 0) {
					muc[c] = muc[c] / (double)params.workNXInCluster(c);
				}
				else {
					muc[c].setZero(nContinuousCovs);
				}
			}

			vector<MatrixXd> Sigmas(maxNClusters);
			for (unsigned int c = 0; c < maxNClusters; c++) {
				Sigmas[c].setZero(nContinuousCovs, nContinuousCovs);
			}

			for (unsigned int i = 0; i<nSubjects; i++) {
				unsigned int zi = params.z(i);
				Sigmas[zi] += (xi[i] - muc[zi])*((xi[i] - muc[zi]).transpose());
			}


			for (unsigned int c = 0; c < maxNClusters; c++) {
				if (params.workNXInCluster(c) > 1) {
					Sigmas[c] = Sigmas[c] / (double)params.workNXInCluster(c);
				}
				else {
					Sigmas[c].setZero(nContinuousCovs, nContinuousCovs);
					VectorXd beta_taus = hyperParams.beta_taus0();
					double alpha_taus = hyperParams.alpha_taus();
					for (unsigned int j = 0; j < nContinuousCovs; j++) {
						randomGamma gammaRand(alpha_taus, 1.0 / beta_taus(j));
						double tausj = gammaRand(rndGenerator);
						Sigmas[c](j, j) = (1.0 / tausj)*(1.0 / tausj);
					}
				}
			}


			for (unsigned int c = 0; c<maxNClusters; c++) {
				for (unsigned int j = 0; j < nContinuousCovs; j++) {
					double Tausj = 1.0 / sqrt(Sigmas[c](j, j));
					params.TauS(c, j, Tausj);
				}
			}

		}
		else {

			vector<MatrixXd> Rc(maxNClusters);
			for (unsigned int c = 0; c<maxNClusters; c++) {
				Rc[c].setZero(nContinuousCovs, nContinuousCovs);
			}

			for (unsigned int i = 0; i<nSubjects; i++) {
				unsigned int zi = params.z(i);
				Rc[zi] = Rc[zi] + (xi[i] - params.mu(zi))*((xi[i] - params.mu(zi)).transpose());
			}

			if (options.useHyperpriorR1()) {
				for (unsigned int c = 0; c<maxNClusters; c++) {
					Rc[c] = (params.R1().inverse() + Rc[c]).inverse();
					MatrixXd Tau = wishartRand(rndGenerator, Rc[c], params.workNXInCluster(c) + params.kappa11());
					params.Tau(c, Tau);
				}
			}
			else {
				for (unsigned int c = 0; c<maxNClusters; c++) {
					Rc[c] = (hyperParams.R0().inverse() + Rc[c]).inverse();
					MatrixXd Tau = wishartRand(rndGenerator, Rc[c], params.workNXInCluster(c) + hyperParams.kappa0());
					params.Tau(c, Tau);
				}

			}

		}

		if (useSeparationPrior) {
			VectorXd sumTau = VectorXd::Zero(nContinuousCovs);
			unsigned int workNactive = 0;
			MatrixXd Tau0 = hyperParams.Tau0();
			for (unsigned int c = 0; c < maxZ + 1; c++) {
				for (unsigned int j = 0; j < nContinuousCovs; j++) {
					sumTau(j) = sumTau(j) + params.TauS(c, j);
				}
				workNactive += 1;
			}

			VectorXd betas_new(nContinuousCovs);
			VectorXd beta_taus0 = hyperParams.beta_taus0();
			for (unsigned int j = 0; j < nContinuousCovs; j++) {
				double alpha_new = workNactive * hyperParams.alpha_taus() + hyperParams.alpha_taus0();
				double beta_new = sumTau(j) + beta_taus0(j);

				randomGamma gammaRand(alpha_new, 1.0 / beta_new);
				betas_new(j) = gammaRand(rndGenerator);

			}
			params.beta_taus(betas_new);

		}


		// Now do the null mu for variable selection
		// In all cases, initialise it at the value it will be fixed at for
		// the continuous indicator case
		if(varSelectType.compare("None")!=0){
			vector<double> meanXVec(nContinuousCovs,0.0);
			vector<unsigned int> countXVec(nContinuousCovs,0);
			for(unsigned int i=0;i<nSubjects;i++){
				for(unsigned int j=0;j<nContinuousCovs;j++){
					if(!dataset.missingX(i,nDiscreteCovs+j)){
						meanXVec[j]+=dataset.continuousX(i,j);
						countXVec[j]+=1;
					}
				}
			}
			VectorXd nullMu=VectorXd::Zero(nContinuousCovs);
			for(unsigned int j=0;j<nContinuousCovs;j++){
				nullMu(j)=meanXVec[j]/(double)countXVec[j];
			}
			params.nullMu(nullMu, useIndependentNormal);
		}


	}

	// Initialise the variable selection variables if appropriate
	// Bias towards having variables in
	if(varSelectType.compare("None")!=0){
		vector<vector<double> > gamma(maxNClusters);
		for(unsigned int c=0;c<maxNClusters;c++){
			gamma[c].resize(nCovariates);
		}
		vector<unsigned int> omega(nCovariates);
		vector<double> rho(nCovariates);
		for(unsigned int j=0;j<nCovariates;j++){
			if((unifRand(rndGenerator)<0.01) && (hyperParams.atomRho()!=1)){
				// We are in the point mass at 0 case - variable is switched off
				omega[j]=0;
				rho[j]=0;
				if(varSelectType.compare("BinaryCluster")==0){
					for(unsigned int c=0;c<maxNClusters;c++){
						gamma[c][j]=0;
					}
				}else{
					gamma[0][j]=0;
				}
			}else{
				omega[j]=1;
				rho[j]=0.75+0.25*unifRand(rndGenerator);
				if(varSelectType.compare("BinaryCluster")==0){
					for(unsigned int c=0;c<maxNClusters;c++){
						if(unifRand(rndGenerator)<rho[j]){
							gamma[c][j]=1;
						}else{
							gamma[c][j]=0;
						}
					}
				}
			}

			params.omega(j,omega[j]);
			params.rho(j,rho[j],covariateType,varSelectType,useIndependentNormal);
			if(varSelectType.compare("BinaryCluster")==0){
				for(unsigned int c=0;c<maxNClusters;c++){
					params.gamma(c,j,gamma[c][j],covariateType,useIndependentNormal);
				}
			}
			// Note in the case of the continuous variable selection indicators
			// gamma is deterministically equal to rho, and so is set in the method
			// for rho so we do nothing here.
		}

	}

	if(includeResponse){
		// Finally we sample the theta and beta values from uniform distributions
		for(unsigned int c=0;c<maxNClusters;c++){
			for (unsigned int k=0;k<nCategoriesY;k++){
				// Thetas are randomly between -2 and 2
				params.theta(c,k,-2.0+4.0*unifRand(rndGenerator));
			}
		}

		for(unsigned int j=0;j<nFixedEffects;j++){
			for (unsigned int k=0;k<nCategoriesY;k++){
				// Betas are randomly between -2 and 2
				params.beta(j,k,-2.0+4.0*unifRand(rndGenerator));
			}
		}

		if(outcomeType.compare("Normal")==0||outcomeType.compare("Quantile")==0){
			randomGamma gammaRand(hyperParams.shapeSigmaSqY(),1.0/hyperParams.scaleSigmaSqY());
			double sigmaSqY=1.0/(gammaRand(rndGenerator));
			params.sigmaSqY(sigmaSqY);
		}

		if(outcomeType.compare("Survival")==0){
			randomGamma gammaRand(hyperParams.shapeNu(),hyperParams.scaleNu());
			if(weibullFixedShape){
				double nu=gammaRand(rndGenerator);
				params.nu(0,nu);

			} else {
				for (unsigned int c=0;c<maxNClusters;c++){
					double nu=gammaRand(rndGenerator);
					params.nu(c,nu);
				}
			}
		}

		// And also the extra variation values if necessary
		if(responseExtraVar){
			// Shape and rate parameters
			double a=5,b=2;
			// Boost parameterised in terms of shape and scale
			randomGamma gammaRand(a,1.0/b);
			// Tau is now a Gamma(a,b)
			double tau = gammaRand(rndGenerator);
			params.tauEpsilon(tau);

			randomNormal normalRand(0,1.0/sqrt(tau));
			for(unsigned int i=0;i<nSubjects;i++){
				double eps = normalRand(rndGenerator);
				int zi = params.z(i);
				double meanVal = params.theta(zi,0);
				if(outcomeType.compare("Categorical")==0){
					for(unsigned int j=0;j<nFixedEffects;j++){
						meanVal+=params.beta(j,dataset.discreteY(i))*dataset.W(i,j);
					}
				} else {
					for(unsigned int j=0;j<nFixedEffects;j++){
						meanVal+=params.beta(j,0)*dataset.W(i,j);
					}
				}
				if(outcomeType.compare("Poisson")==0){
					meanVal+=dataset.logOffset(i);
				}
				params.lambda(i,meanVal+eps);

			}
		}
	        // And also _uCAR and _TauCAR if includeCAR==TRUE
	        if (includeCAR){
			// Boost parameterised in terms of shape and scale
			randomGamma gammaRand(5.0,0.5);
			// Tau is now a Gamma(shape,rate)
			double tau = gammaRand(rndGenerator);
			params.TauCAR(tau);

			if (includeuCARinit) {
			        ifstream uCARfile;
			        uCARfile.open(uCARinitFileName.c_str());	
			        if (!uCARfile.is_open()){
					Rprintf("File with initialisation values for uCAR not found\n");
					wasError = true;
        			}
				if (uCARfile.good()){
					for(unsigned int i=0;i<nSubjects;i++){
						string line;
						getline(uCARfile, line);
						stringstream streamline(line);
						double eps;
						streamline >> eps;
						params.uCAR(i,eps);
					}
				}
				uCARfile.close();
			} else {	
				double mean_w=0;
				for (unsigned int i=0; i<nSubjects; i++ ) mean_w+=dataset.nNeighbours(i);
				mean_w /= nSubjects;
				randomNormal normalRand(0,sqrt(mean_w/tau));
				for(unsigned int i=0;i<nSubjects;i++){
					double eps = normalRand(rndGenerator);
					params.uCAR(i,eps);
				}
			}

	        }

	}

	if(wasError){
		Rprintf("There is a mistake in the initialisation of PReMiuM.\n");
	}

}

// Write the sampler output
void writePReMiuMOutput(mcmcSampler<pReMiuMParams,pReMiuMOptions,pReMiuMPropParams,pReMiuMData>& sampler,
								const unsigned int& sweep){


	bool reportBurnIn = sampler.reportBurnIn();
	unsigned int nBurn = sampler.nBurn();
	unsigned int nFilter = sampler.nFilter();
	vector<ofstream*>& outFiles = sampler.outFiles();

	// Check if we need to do anything
	if((reportBurnIn||((!reportBurnIn)&&sweep>nBurn))&&(sweep%nFilter==0)){
		const pReMiuMParams& params = sampler.chain().currentState().parameters();

		unsigned int nSubjects = params.nSubjects();
		unsigned int nPredictSubjects = params.nPredictSubjects();
		unsigned int maxNClusters = params.maxNClusters();
		unsigned int nCovariates = params.nCovariates();
		unsigned int nDiscreteCovs=params.nDiscreteCovs();
		unsigned int nContinuousCovs=params.nContinuousCovs();
		unsigned int nCategoriesY = params.nCategoriesY();
		string covariateType = sampler.model().dataset().covariateType();
		bool includeResponse = sampler.model().options().includeResponse();
		bool includeCAR = sampler.model().options().includeCAR();
		bool responseExtraVar = sampler.model().options().responseExtraVar();
		double fixedAlpha = sampler.model().options().fixedAlpha();
		string outcomeType = sampler.model().options().outcomeType();
		bool computeEntropy = sampler.model().options().computeEntropy();
		unsigned int nFixedEffects = params.nFixedEffects(outcomeType);
		string varSelectType = sampler.model().options().varSelectType();
		string predictType = sampler.model().options().predictType();
		bool weibullFixedShape = sampler.model().options().weibullFixedShape();
		bool useHyperpriorR1 = sampler.model().options().useHyperpriorR1();
		bool useIndependentNormal= sampler.model().options().useIndependentNormal();
		bool useSeparationPrior = sampler.model().options().useSeparationPrior();

		const pReMiuMData& dataset = sampler.model().dataset();
		pReMiuMPropParams& proposalParams = sampler.proposalParams();

		vector<unsigned int> nCategories;
		if(covariateType.compare("Discrete")==0||covariateType.compare("Mixed")==0){
			nCategories = params.nCategories();
		}

		// Check if the files are already open
		if(outFiles.size()==0){
			string fileStem =sampler.outFileStem();
			string fileName = fileStem + "_nClusters.txt";
			outFiles.push_back(new ofstream(fileName.c_str()));
			fileName = fileStem + "_psi.txt";
			outFiles.push_back(new ofstream(fileName.c_str()));
			if(covariateType.compare("Discrete")==0){
				fileName = fileStem + "_phi.txt";
				outFiles.push_back(new ofstream(fileName.c_str()));
			}else if(covariateType.compare("Normal")==0){
				fileName = fileStem + "_mu.txt";
				outFiles.push_back(new ofstream(fileName.c_str()));
				fileName = fileStem + "_Sigma.txt";
				outFiles.push_back(new ofstream(fileName.c_str()));
				if (useHyperpriorR1||useIndependentNormal) {
					fileName = fileStem + "_R1.txt";	
					outFiles.push_back(new ofstream(fileName.c_str()));
				}
				if (useHyperpriorR1 ) {
					fileName = fileStem + "_kappa1.txt";
					outFiles.push_back(new ofstream(fileName.c_str()));

					fileName = fileStem + "_mu00.txt";
					outFiles.push_back(new ofstream(fileName.c_str()));

					fileName = fileStem + "_Sigma00.txt";
					outFiles.push_back(new ofstream(fileName.c_str()));

					fileName = fileStem + "_kappa1Prop.txt";
					outFiles.push_back(new ofstream(fileName.c_str()));
				} 

				if (useSeparationPrior) {
					fileName = fileStem + "_mu00.txt";
					outFiles.push_back(new ofstream(fileName.c_str()));

					fileName = fileStem + "_SigmaR.txt";
					outFiles.push_back(new ofstream(fileName.c_str()));

					fileName = fileStem + "_SigmaS.txt";
					outFiles.push_back(new ofstream(fileName.c_str()));

					fileName = fileStem + "_SigmaSProp.txt";
					outFiles.push_back(new ofstream(fileName.c_str()));

					fileName = fileStem + "_kappa1.txt";
					outFiles.push_back(new ofstream(fileName.c_str()));

					fileName = fileStem + "_kappa1Prop.txt";
					outFiles.push_back(new ofstream(fileName.c_str()));

				}

			}else if(covariateType.compare("Mixed")==0){
				fileName = fileStem + "_phi.txt";
				outFiles.push_back(new ofstream(fileName.c_str()));
				fileName = fileStem + "_mu.txt";
				outFiles.push_back(new ofstream(fileName.c_str()));
				fileName = fileStem + "_Sigma.txt";
				outFiles.push_back(new ofstream(fileName.c_str()));
				if (useHyperpriorR1|| useIndependentNormal) {
					fileName = fileStem + "_R1.txt";	
					outFiles.push_back(new ofstream(fileName.c_str()));
				}
				if (useHyperpriorR1) {
					fileName = fileStem + "_kappa1.txt";
					outFiles.push_back(new ofstream(fileName.c_str()));

					fileName = fileStem + "_mu00.txt";
					outFiles.push_back(new ofstream(fileName.c_str()));

					fileName = fileStem + "_Sigma00.txt";
					outFiles.push_back(new ofstream(fileName.c_str()));

					fileName = fileStem + "_kappa1Prop.txt";
					outFiles.push_back(new ofstream(fileName.c_str()));
				} 

				if (useSeparationPrior) {
					fileName = fileStem + "_mu00.txt";
					outFiles.push_back(new ofstream(fileName.c_str()));

					fileName = fileStem + "_SigmaR.txt";
					outFiles.push_back(new ofstream(fileName.c_str()));

					fileName = fileStem + "_SigmaS.txt";
					outFiles.push_back(new ofstream(fileName.c_str()));

					fileName = fileStem + "_SigmaSProp.txt";
					outFiles.push_back(new ofstream(fileName.c_str()));

					fileName = fileStem + "_kappa1.txt";
					outFiles.push_back(new ofstream(fileName.c_str()));

					fileName = fileStem + "_kappa1Prop.txt";
					outFiles.push_back(new ofstream(fileName.c_str()));

				}
			}
			fileName = fileStem + "_z.txt";
			outFiles.push_back(new ofstream(fileName.c_str()));
			fileName = fileStem + "_entropy.txt";
			outFiles.push_back(new ofstream(fileName.c_str()));
			fileName = fileStem + "_alpha.txt";
			outFiles.push_back(new ofstream(fileName.c_str()));
			fileName = fileStem + "_logPost.txt";
			outFiles.push_back(new ofstream(fileName.c_str()));
			fileName = fileStem + "_nMembers.txt";
			outFiles.push_back(new ofstream(fileName.c_str()));
			if(fixedAlpha<=-1){
				fileName = fileStem + "_alphaProp.txt";
				outFiles.push_back(new ofstream(fileName.c_str()));
			}
			if(includeResponse){
				fileName = fileStem + "_theta.txt";
				outFiles.push_back(new ofstream(fileName.c_str()));
				fileName = fileStem + "_beta.txt";
				outFiles.push_back(new ofstream(fileName.c_str()));
				fileName = fileStem + "_thetaProp.txt";
				outFiles.push_back(new ofstream(fileName.c_str()));
				fileName = fileStem + "_betaProp.txt";
				outFiles.push_back(new ofstream(fileName.c_str()));
				if(outcomeType.compare("Normal")==0||outcomeType.compare("Quantile")==0){
					fileName = fileStem + "_sigmaSqY.txt";
					outFiles.push_back(new ofstream(fileName.c_str()));
				}
				if(outcomeType.compare("Survival")==0){
					fileName = fileStem + "_nu.txt";
					outFiles.push_back(new ofstream(fileName.c_str()));
				}
				if(responseExtraVar){
					fileName = fileStem + "_epsilon.txt";
					outFiles.push_back(new ofstream(fileName.c_str()));
					fileName = fileStem + "_sigmaEpsilon.txt";
					outFiles.push_back(new ofstream(fileName.c_str()));
					fileName = fileStem + "_epsilonProp.txt";
					outFiles.push_back(new ofstream(fileName.c_str()));
				}
				if(nPredictSubjects>0){
					fileName = fileStem + "_predictThetaRaoBlackwell.txt";
					outFiles.push_back(new ofstream(fileName.c_str()));
				}
				if (includeCAR){
					fileName = fileStem + "_TauCAR.txt";
					outFiles.push_back(new ofstream(fileName.c_str()));
					fileName = fileStem + "_uCAR.txt";
					outFiles.push_back(new ofstream(fileName.c_str()));
				}
			}
			if(varSelectType.compare("None")!=0){
				fileName = fileStem + "_omega.txt";
				outFiles.push_back(new ofstream(fileName.c_str()));
				fileName = fileStem + "_rho.txt";
				outFiles.push_back(new ofstream(fileName.c_str()));
				fileName = fileStem + "_rhoOmegaProp.txt";
				outFiles.push_back(new ofstream(fileName.c_str()));
				if(varSelectType.compare("Continuous")!=0){
					fileName = fileStem + "_gamma.txt";
					outFiles.push_back(new ofstream(fileName.c_str()));
				}
				if(covariateType.compare("Discrete")==0){
					fileName = fileStem + "_nullPhi.txt";
					outFiles.push_back(new ofstream(fileName.c_str()));
				}else if(covariateType.compare("Normal")==0){
					fileName = fileStem + "_nullMu.txt";
					outFiles.push_back(new ofstream(fileName.c_str()));
				}else if(covariateType.compare("Mixed")==0){
					fileName = fileStem + "_nullPhi.txt";
					outFiles.push_back(new ofstream(fileName.c_str()));
					fileName = fileStem + "_nullMu.txt";
					outFiles.push_back(new ofstream(fileName.c_str()));
				}
			}
		}

		// File indices
		int nClustersInd=-1,psiInd=-1,phiInd=-1,muInd=-1,SigmaInd=-1,R1Ind=-1,zInd=-1,entropyInd=-1,alphaInd=-1,kappa1Ind=-1,mu00Ind=-1, Sigma00Ind = -1;
		int logPostInd = -1, nMembersInd = -1, alphaPropInd = -1, kappa1PropInd = -1, SigmaRInd = -1, SigmaSInd = -1, SigmaSPropInd = -1;
		int thetaInd=-1,betaInd=-1,thetaPropInd=-1,betaPropInd=-1,sigmaSqYInd=-1,nuInd=-1,epsilonInd=-1;
		int sigmaEpsilonInd=-1,epsilonPropInd=-1,omegaInd=-1,rhoInd=-1;
		int rhoOmegaPropInd=-1,gammaInd=-1,nullPhiInd=-1,nullMuInd=-1;
		int predictThetaRaoBlackwellInd=-1;
		int TauCARInd=-1,uCARInd=-1;

		int r=0;
		nClustersInd=r++;
		psiInd=r++;
		if(covariateType.compare("Discrete")==0){
			phiInd=r++;
		}else if(covariateType.compare("Normal")==0){
			muInd=r++;
			SigmaInd=r++;
			if (useHyperpriorR1||useIndependentNormal) R1Ind=r++;
			if (useHyperpriorR1) {
				kappa1Ind = r++;
				mu00Ind = r++;
				Sigma00Ind = r++;
				kappa1PropInd = r++;
			} 
			if (useSeparationPrior) {
				mu00Ind = r++;
				SigmaRInd = r++;
				SigmaSInd = r++;
				SigmaSPropInd = r++;
				kappa1Ind = r++;
				kappa1PropInd = r++;
			}
		}else if(covariateType.compare("Mixed")==0){
			phiInd=r++;
			muInd=r++;
			SigmaInd=r++;
			if (useHyperpriorR1 || useIndependentNormal) R1Ind=r++;
			if (useHyperpriorR1) {
				kappa1Ind = r++;
				mu00Ind = r++;
				Sigma00Ind = r++;
				kappa1PropInd = r++;
			} 
			if (useSeparationPrior) {
				mu00Ind = r++;
				SigmaRInd = r++;
				SigmaSInd = r++;
				SigmaSPropInd = r++;
				kappa1Ind = r++;
				kappa1PropInd = r++;
			}
		}
		zInd=r++;
		entropyInd=r++;
		alphaInd=r++;
		logPostInd=r++;
		nMembersInd=r++;
		if(fixedAlpha<=-1){
			alphaPropInd=r++;
		}

		if(includeResponse){
			thetaInd=r++;
			betaInd=r++;
			thetaPropInd=r++;
			betaPropInd=r++;
			if(outcomeType.compare("Normal")==0||outcomeType.compare("Quantile")==0){
				sigmaSqYInd=r++;
			}
			if(outcomeType.compare("Survival")==0){
				nuInd=r++;
			}
			if(responseExtraVar){
				epsilonInd=r++;
				sigmaEpsilonInd=r++;
				epsilonPropInd=r++;
			}
			if(nPredictSubjects>0){
				predictThetaRaoBlackwellInd=r++;
			}
			if (includeCAR){
				TauCARInd=r++;
				uCARInd=r++;
			}
		}

		if(varSelectType.compare("None")!=0){
			omegaInd=r++;
			rhoInd=r++;
			rhoOmegaPropInd=r++;
			if(varSelectType.compare("Continuous")!=0){
				gammaInd=r++;
			}
			if(covariateType.compare("Discrete")==0){
				nullPhiInd=r++;
			}else if (covariateType.compare("Normal")==0){
				nullMuInd=r++;
			}else if (covariateType.compare("Mixed")==0){
				nullPhiInd=r++;
				nullMuInd=r++;
			}
		}


		*(outFiles[nClustersInd]) << maxNClusters << endl;

		unsigned int sumMembers=0;
		for(unsigned int c=0;c<maxNClusters;c++){
			// Print logPsi
			*(outFiles[psiInd]) << exp(params.logPsi(c));
			if(includeResponse){
				// Print theta
				if(outcomeType.compare("Categorical")==0){
					for (unsigned int k=0;k<nCategoriesY;k++){
						*(outFiles[thetaInd]) << params.theta(c,k);
						if (k<(nCategoriesY-1)) {
							*(outFiles[thetaInd]) <<" ";
						}
					}
				} else {
					*(outFiles[thetaInd]) << params.theta(c,0);
				}
			}
			// Print number of members of each cluster
			unsigned int nXinC = params.workNXInCluster(c);
			*(outFiles[nMembersInd]) << nXinC;
			sumMembers +=nXinC;
			if(c<maxNClusters-1){
				*(outFiles[psiInd]) << " ";
				if(includeResponse){
					*(outFiles[thetaInd]) << " ";
				}
				*(outFiles[nMembersInd]) << " ";
			}else{
				*(outFiles[psiInd]) << endl;
				if(includeResponse){
					*(outFiles[thetaInd]) << endl;
				}
				*(outFiles[nMembersInd]) << " " << sumMembers << endl;
			}
		}


		unsigned int maxNCategories=0;

		if(covariateType.compare("Discrete")==0){
			for(unsigned int j=0;j<nCovariates;j++){
				if(nCategories[j]>maxNCategories){
					maxNCategories=nCategories[j];
				}
			}

			for(unsigned int j=0;j<nCovariates;j++){
				for(unsigned int p=0;p<maxNCategories;p++){
					for(unsigned int c=0;c<maxNClusters;c++){
						// Print Phi
						if(p<nCategories[j]){
							*(outFiles[phiInd]) << exp(params.logPhi(c,j,p));
						}else{
							// pad the output with dummy variables
							// to make reading in R easier
							*(outFiles[phiInd]) << -999;
						}
						if(c<(maxNClusters-1)||p<(maxNCategories-1)||j<(nCovariates-1)){
							*(outFiles[phiInd]) << " ";
						}

					}
				}
			}
			*(outFiles[phiInd]) << endl;
		}else if(covariateType.compare("Normal")==0){
			// To make the output comparable with discrete, we will write the
			// output grouped by covariate (for each cluster)
			for(unsigned int j=0;j<nCovariates;j++){
				for(unsigned int c=0;c<maxNClusters;c++){
					*(outFiles[muInd]) << params.mu(c,j);
					if(c<(maxNClusters-1)||j<(nCovariates-1)){
						*(outFiles[muInd]) << " ";
					}
				}
			}
			*(outFiles[muInd]) << endl;


			// For the covariance matrices we write by covariate x covariate (for each cluster)

			if (useIndependentNormal) {
				for (unsigned int j = 0; j<nCovariates; j++) {
					for (unsigned int c = 0; c<maxNClusters; c++) {
						*(outFiles[SigmaInd]) << params.Sigma_Indep(c, j);
						if (c<(maxNClusters - 1) || j<(nCovariates - 1)) {
							*(outFiles[SigmaInd]) << " ";
						}
					}
				}
				*(outFiles[SigmaInd]) << endl;

				for (unsigned int j = 0; j<nCovariates; j++) {
					*(outFiles[R1Ind]) << params.R1_Indep(j);
						if (j<(nCovariates - 1)) {
							*(outFiles[R1Ind]) << " ";
						}
				}
				*(outFiles[R1Ind]) << endl;


			}
			else {
				for (unsigned int j1 = 0; j1<nCovariates; j1++) {
					for (unsigned int j2 = 0; j2<nCovariates; j2++) {
						for (unsigned int c = 0; c<maxNClusters; c++) {
							*(outFiles[SigmaInd]) << params.Sigma(c, j1, j2);
							if (c<(maxNClusters - 1) || j1<(nCovariates - 1) || j2<(nCovariates - 1)) {
								*(outFiles[SigmaInd]) << " ";
							}
						}
					}
				}
				*(outFiles[SigmaInd]) << endl;
			}
				
			
			if (useHyperpriorR1){
				for(unsigned int j1=0;j1<nCovariates;j1++){
					for(unsigned int j2=0;j2<nCovariates;j2++){
						*(outFiles[R1Ind]) << params.R1(j1,j2);
						if(j1<(nCovariates-1)||j2<(nCovariates-1)){
							*(outFiles[R1Ind]) << " ";
						}
					}
				}
				*(outFiles[R1Ind]) << endl;

				*(outFiles[kappa1Ind]) << params.kappa11() << endl;

				bool anyUpdates = proposalParams.kappa1AnyUpdates();
				if (anyUpdates) {
					*(outFiles[kappa1PropInd]) << sampler.proposalParams().kappa1AcceptRate() <<
						" " << sampler.proposalParams().kappa1StdDev() << endl;
					proposalParams.kappa1AnyUpdates(false);
				}

				for (unsigned int j = 0; j<nCovariates; j++) {
					*(outFiles[mu00Ind]) << params.mu00(j);
					if (j<(nCovariates - 1)) {
						*(outFiles[mu00Ind]) << " ";
					}
				}
				*(outFiles[mu00Ind]) << endl;

				for (unsigned int j1 = 0; j1<nCovariates; j1++) {
					for (unsigned int j2 = 0; j2<nCovariates; j2++) {
						*(outFiles[Sigma00Ind]) << params.Sigma00(j1, j2);
						if (j1<(nCovariates - 1) || j2<(nCovariates - 1)) {
							*(outFiles[Sigma00Ind]) << " ";
						}
					}
				}

				*(outFiles[Sigma00Ind]) << endl;
			}

			if (useSeparationPrior) {

				for (unsigned int j = 0; j<nCovariates; j++) {
					*(outFiles[mu00Ind]) << params.mu00(j);
					if (j<(nCovariates - 1)) {
						*(outFiles[mu00Ind]) << " ";
					}
				}
				*(outFiles[mu00Ind]) << endl;

				for (unsigned int j1 = 0; j1<nCovariates; j1++) {
					for (unsigned int j2 = 0; j2<nCovariates; j2++) {
						for (unsigned int c = 0; c<maxNClusters; c++) {
							*(outFiles[SigmaRInd]) << params.SigmaR(c, j1, j2);
							if (c<(maxNClusters - 1) || j1<(nCovariates - 1) || j2<(nCovariates - 1)) {
								*(outFiles[SigmaRInd]) << " ";
							}
						}
					}
				}
				*(outFiles[SigmaRInd]) << endl;

				for (unsigned int j = 0; j<nCovariates; j++) {
					for (unsigned int c = 0; c<maxNClusters; c++) {
						*(outFiles[SigmaSInd]) << params.SigmaS(c, j);
						if (c<(maxNClusters - 1) || j<(nCovariates - 1)) {
							*(outFiles[SigmaSInd]) << " ";
						}
					}
				}
				*(outFiles[SigmaSInd]) << endl;

				// Print the acceptance rates for SigmaS
				bool anyUpdates = proposalParams.TauSAnyUpdates();
				if (anyUpdates) {
					for (unsigned int j = 0; j < nCovariates; j++) {
						*(outFiles[SigmaSPropInd]) << sampler.proposalParams().TauSAcceptRate(j) <<
							" " << sampler.proposalParams().TauSStdDev(j);
						if (j < (nCovariates - 1)) {
							*(outFiles[SigmaSPropInd]) << endl;
						}
						*(outFiles[SigmaSPropInd]) << endl;
						proposalParams.TauSAnyUpdates(false);
					}
				}


				*(outFiles[kappa1Ind]) << params.kappa11() << endl;

				anyUpdates = proposalParams.kappa1AnyUpdates();
				if (anyUpdates) {
					*(outFiles[kappa1PropInd]) << sampler.proposalParams().kappa1AcceptRate() <<
						" " << sampler.proposalParams().kappa1StdDev() << endl;
					proposalParams.kappa1AnyUpdates(false);
				}

			}


		}else if(covariateType.compare("Mixed")==0){
			for(unsigned int j=0;j<nDiscreteCovs;j++){
				if(nCategories[j]>maxNCategories){
					maxNCategories=nCategories[j];
				}
			}

			for(unsigned int j=0;j<nDiscreteCovs;j++){
				for(unsigned int p=0;p<maxNCategories;p++){
					for(unsigned int c=0;c<maxNClusters;c++){
						// Print Phi
						if(p<nCategories[j]){
							*(outFiles[phiInd]) << exp(params.logPhi(c,j,p));
						}else{
							// pad the output with dummy variables
							// to make reading in R easier
							*(outFiles[phiInd]) << -999;
						}
						if(c<(maxNClusters-1)||p<(maxNCategories-1)||j<(nDiscreteCovs-1)){
							*(outFiles[phiInd]) << " ";
						}

					}
				}
			}
			*(outFiles[phiInd]) << endl;

			// To make the output comparable with discrete, we will write the
			// output grouped by covariate (for each cluster)
			for(unsigned int j=0;j<nContinuousCovs;j++){
				for(unsigned int c=0;c<maxNClusters;c++){
					*(outFiles[muInd]) << params.mu(c,j);
					if(c<(maxNClusters-1)||j<(nContinuousCovs-1)){
						*(outFiles[muInd]) << " ";
					}
				}
			}
			*(outFiles[muInd]) << endl;


			// For the covariance matrices we write by covariate x covariate (for each cluster)

			if (useIndependentNormal) {
				for (unsigned int j = 0; j<nContinuousCovs; j++) {
					for (unsigned int c = 0; c<maxNClusters; c++) {
						*(outFiles[SigmaInd]) << params.Sigma_Indep(c, j);
						if (c<(maxNClusters - 1) || j<(nContinuousCovs - 1)) {
							*(outFiles[SigmaInd]) << " ";
						}
					}
				}
				*(outFiles[SigmaInd]) << endl;

				for (unsigned int j = 0; j<nContinuousCovs; j++) {
					*(outFiles[R1Ind]) << params.R1_Indep(j);
					if (j<(nContinuousCovs - 1)) {
						*(outFiles[R1Ind]) << " ";
					}
				}
				*(outFiles[R1Ind]) << endl;

			}
			else {
				for (unsigned int j1 = 0; j1<nContinuousCovs; j1++) {
					for (unsigned int j2 = 0; j2<nContinuousCovs; j2++) {
						for (unsigned int c = 0; c<maxNClusters; c++) {
							*(outFiles[SigmaInd]) << params.Sigma(c, j1, j2);
							if (c<(maxNClusters - 1) || j1<(nContinuousCovs - 1) || j2<(nContinuousCovs - 1)) {
								*(outFiles[SigmaInd]) << " ";
							}
						}
					}
				}
				*(outFiles[SigmaInd]) << endl;

			}

			if (useHyperpriorR1){
				for(unsigned int j1=0;j1<nContinuousCovs;j1++){
					for(unsigned int j2=0;j2<nContinuousCovs;j2++){
						*(outFiles[R1Ind]) << params.R1(j1,j2);
						if(j1<(nContinuousCovs-1)||j2<(nContinuousCovs-1)){
							*(outFiles[R1Ind]) << " ";
						}
					}
				}
			*(outFiles[R1Ind]) << endl;
		
			*(outFiles[kappa1Ind]) << params.kappa11() << endl;

			bool anyUpdates = proposalParams.kappa1AnyUpdates();
			if (anyUpdates) {
				*(outFiles[kappa1PropInd]) << sampler.proposalParams().kappa1AcceptRate() <<
					" " << sampler.proposalParams().kappa1StdDev() << endl;
				proposalParams.kappa1AnyUpdates(false);
			}

			for (unsigned int j = 0; j<nContinuousCovs; j++) {
				*(outFiles[mu00Ind]) << params.mu00(j);
				if (j<(nContinuousCovs - 1)) {
					*(outFiles[mu00Ind]) << " ";
				}
			}
			*(outFiles[mu00Ind]) << endl;

			for (unsigned int j1 = 0; j1<nContinuousCovs; j1++) {
				for (unsigned int j2 = 0; j2<nContinuousCovs; j2++) {
					*(outFiles[Sigma00Ind]) << params.Sigma00(j1, j2);
					if (j1<(nContinuousCovs - 1) || j2<(nContinuousCovs - 1)) {
						*(outFiles[Sigma00Ind]) << " ";
					}
				}
			}
			*(outFiles[Sigma00Ind]) << endl;
				
			}


			if (useSeparationPrior) {

				for (unsigned int j = 0; j<nContinuousCovs; j++) {
					*(outFiles[mu00Ind]) << params.mu00(j);
					if (j<(nContinuousCovs - 1)) {
						*(outFiles[mu00Ind]) << " ";
					}
				}
				*(outFiles[mu00Ind]) << endl;

				for (unsigned int j1 = 0; j1<nContinuousCovs; j1++) {
					for (unsigned int j2 = 0; j2<nContinuousCovs; j2++) {
						for (unsigned int c = 0; c<maxNClusters; c++) {
							*(outFiles[SigmaRInd]) << params.SigmaR(c, j1, j2);
							if (c<(maxNClusters - 1) || j1<(nContinuousCovs - 1) || j2<(nContinuousCovs - 1)) {
								*(outFiles[SigmaRInd]) << " ";
							}
						}
					}
				}
				*(outFiles[SigmaRInd]) << endl;

				for (unsigned int j = 0; j<nContinuousCovs; j++) {
					for (unsigned int c = 0; c<maxNClusters; c++) {
						*(outFiles[SigmaSInd]) << params.SigmaS(c, j);
						if (c<(maxNClusters - 1) || j<(nContinuousCovs - 1)) {
							*(outFiles[SigmaSInd]) << " ";
						}
					}
				}
				*(outFiles[SigmaSInd]) << endl;

				// Print the acceptance rates for SigmaS
				bool anyUpdates = proposalParams.TauSAnyUpdates();
				if (anyUpdates) {
					for (unsigned int j = 0; j < nContinuousCovs; j++) {
						*(outFiles[SigmaSPropInd]) << sampler.proposalParams().TauSAcceptRate(j) <<
							" " << sampler.proposalParams().TauSStdDev(j);
						if (j < (nContinuousCovs - 1)) {
							*(outFiles[SigmaSPropInd]) << endl;
						}
						*(outFiles[SigmaSPropInd]) << endl;
						proposalParams.TauSAnyUpdates(false);
					}
				}


				*(outFiles[kappa1Ind]) << params.kappa11() << endl;

				anyUpdates = proposalParams.kappa1AnyUpdates();
				if (anyUpdates) {
					*(outFiles[kappa1PropInd]) << sampler.proposalParams().kappa1AcceptRate() <<
						" " << sampler.proposalParams().kappa1StdDev() << endl;
					proposalParams.kappa1AnyUpdates(false);
				}
			}

		}


		if(includeResponse){
			if(outcomeType.compare("Categorical")==0){
				// Print beta
				for(unsigned int j=0;j<nFixedEffects;j++){
					for (unsigned int k=0;k<nCategoriesY;k++){
						*(outFiles[betaInd]) << params.beta(j,k) <<" ";
					}
					if(j==(nFixedEffects-1)){
						*(outFiles[betaInd]) << endl;
					}
				}
			} else {
				if(outcomeType.compare("Normal")==0||outcomeType.compare("Quantile")==0){
					*(outFiles[sigmaSqYInd]) << params.sigmaSqY() << endl;
				}
				if(outcomeType.compare("Survival")==0){
				// Print parameter nu for each cluster
					if (weibullFixedShape){
						*(outFiles[nuInd]) << params.nu(0) << endl;
					} else {
						for(unsigned int c=0;c<maxNClusters;c++){
							*(outFiles[nuInd]) << params.nu(c);
							if(c<maxNClusters-1){
								*(outFiles[nuInd]) << " ";
							}else{
								*(outFiles[nuInd]) << endl;
							}
						}				
					}
				}
				// Print beta
				for(unsigned int j=0;j<nFixedEffects;j++){
					*(outFiles[betaInd]) << params.beta(j,0);
					if(j<nFixedEffects-1){
						*(outFiles[betaInd]) << " ";
					}else{
						*(outFiles[betaInd]) << endl;
					}
				}
			}
			if (includeCAR){
				for(unsigned int i=0;i<nSubjects;i++){
					double uCARi = params.uCAR(i);
					*(outFiles[uCARInd]) << uCARi;
					if(i<nSubjects-1){
						*(outFiles[uCARInd]) << " ";
					}else{
						*(outFiles[uCARInd]) << endl;
					}
				}
				double tau=params.TauCAR();
				*(outFiles[TauCARInd]) << tau;
				*(outFiles[TauCARInd]) << endl;
			}
		}

		// Print the allocations
		for(unsigned int i=0;i<nSubjects+nPredictSubjects;i++){
			*(outFiles[zInd]) << params.z(i);
			if(computeEntropy){
				*(outFiles[entropyInd]) << params.workEntropy(i);
			}
			if(i<nSubjects+nPredictSubjects-1){
				if(computeEntropy){
					*(outFiles[entropyInd]) << " ";
				}
				*(outFiles[zInd]) << " ";
			}else{
				if(computeEntropy){
					*(outFiles[entropyInd]) << endl;
				}
				*(outFiles[zInd]) << endl;
			}
			// And print the expected theta for the prediction subjects
			if(i>=nSubjects){
				if(includeResponse){
					for (unsigned int k=0;k<nCategoriesY;k++){
						*(outFiles[predictThetaRaoBlackwellInd]) << params.workPredictExpectedTheta(i-nSubjects,k)<<" ";
					}
					if(i<nSubjects+nPredictSubjects-1){
						*(outFiles[predictThetaRaoBlackwellInd]) << " ";
					}else{
						*(outFiles[predictThetaRaoBlackwellInd]) << endl;
					}
				}
			}
		}

		// Print alpha
		if(fixedAlpha<=-1||sweep==0){
			*(outFiles[alphaInd]) << params.alpha() << endl;
		}

		// Print the log posterior
		*(outFiles[logPostInd]) << sampler.chain().currentState().logPosterior() << " ";
		*(outFiles[logPostInd]) << sampler.chain().currentState().logLikelihood() << " ";
		*(outFiles[logPostInd]) << sampler.chain().currentState().logPrior() << endl;

		bool anyUpdates;
		if(includeResponse){
			// Print the acceptance rates for theta
			anyUpdates = proposalParams.thetaAnyUpdates();
			if(anyUpdates){
				*(outFiles[thetaPropInd]) << sampler.proposalParams().thetaAcceptRate() <<
					" " << sampler.proposalParams().thetaStdDev() << endl;
				proposalParams.thetaAnyUpdates(false);
			}


			// Print the acceptance rates for beta
			anyUpdates = proposalParams.betaAnyUpdates();
			if(anyUpdates){
				for(unsigned int j=0;j<nFixedEffects;j++){
					*(outFiles[betaPropInd]) << sampler.proposalParams().betaAcceptRate(j) <<
							" " << sampler.proposalParams().betaStdDev(j);
					if(j<(nFixedEffects-1)){
						*(outFiles[betaPropInd]) << endl;
					}
					*(outFiles[betaPropInd]) << endl;
					proposalParams.betaAnyUpdates(false);
				}
				if(responseExtraVar){
					vector<double> meanVec(nSubjects,0.0);
					if(outcomeType.compare("Poisson")==0){
						meanVec=dataset.logOffset();
					}
					for(unsigned int i=0;i<nSubjects;i++){
						int zi = params.z(i);
						double meanVal = meanVec[i]+params.theta(zi,0);
						for(unsigned int j=0;j<nFixedEffects;j++){
							meanVal+=params.beta(j,0)*dataset.W(i,j);
						}
						double eps=params.lambda(i)-meanVal;
						*(outFiles[epsilonInd]) << eps;
						if(i<nSubjects-1){
							*(outFiles[epsilonInd]) << " ";
						}else{
							*(outFiles[epsilonInd]) << endl;
						}
					}
					anyUpdates = proposalParams.lambdaAnyUpdates();
					if(anyUpdates){
						*(outFiles[epsilonPropInd]) << sampler.proposalParams().lambdaAcceptRate() <<
								" " << sampler.proposalParams().lambdaStdDev() << endl;
						proposalParams.lambdaAnyUpdates(false);
					}
					*(outFiles[sigmaEpsilonInd]) << 1.0/sqrt(params.tauEpsilon()) << endl;
				}
			}
		}

		// Print the acceptance rates for alpha
		if(fixedAlpha<=-1){
			anyUpdates = proposalParams.alphaAnyUpdates();
			if(anyUpdates){
				*(outFiles[alphaPropInd]) << sampler.proposalParams().alphaAcceptRate() <<
							" " << sampler.proposalParams().alphaStdDev() << endl;
				proposalParams.alphaAnyUpdates(false);
			}
		}

		if(varSelectType.compare("None")!=0){
			// Print variable selection related quantities
			for(unsigned int j=0;j<nCovariates;j++){
				*(outFiles[omegaInd]) << params.omega(j);
				*(outFiles[rhoInd]) << params.rho(j);
				if(j<nCovariates-1){
					*(outFiles[omegaInd]) << " ";
					*(outFiles[rhoInd]) << " ";

				}else{
					*(outFiles[omegaInd]) << endl;
					*(outFiles[rhoInd]) << endl;
				}
				if(sweep!=0){ // this was "==0", it might be worth double checking 
					if(covariateType.compare("Discrete")==0){
						for(unsigned int p=0;p<maxNCategories;p++){
							if(p<nCategories[j]){
								*(outFiles[nullPhiInd]) << exp(params.logNullPhi(j,p));
							}else{
								*(outFiles[nullPhiInd]) << -999;
							}

							if(p<(maxNCategories-1)||j<(nCovariates-1)){
								*(outFiles[nullPhiInd]) << " ";
							}else{
								*(outFiles[nullPhiInd]) << endl;
							}

						}
					}else if(covariateType.compare("Normal")==0){
						*(outFiles[nullMuInd]) << params.nullMu(j);
						if(j<nCovariates-1){
							*(outFiles[nullMuInd]) << " ";
						}else{
							*(outFiles[nullMuInd]) << endl;
						}

					}else if(covariateType.compare("Mixed")==0){
						if (j < nDiscreteCovs){
							for(unsigned int p=0;p<maxNCategories;p++){
								if(p<nCategories[j]){
									*(outFiles[nullPhiInd]) << exp(params.logNullPhi(j,p));
								}else{
									*(outFiles[nullPhiInd]) << -999;
								}
								if(p<(maxNCategories-1)||j<(nDiscreteCovs-1)){
									*(outFiles[nullPhiInd]) << " ";
								}else{
									*(outFiles[nullPhiInd]) << endl;
								}
							}
						} else {
							*(outFiles[nullMuInd]) << params.nullMu(j-nDiscreteCovs);
							if(j<nCovariates-1){
								*(outFiles[nullMuInd]) << " ";
							}else{
								*(outFiles[nullMuInd]) << endl;
							}
						}
					}
				}
				if(varSelectType.compare("BinaryCluster")==0){
					for(unsigned int c=0;c<maxNClusters;c++){
						*(outFiles[gammaInd]) << params.gamma(c,j);
						if(c<maxNClusters-1||j<nCovariates-1){
							*(outFiles[gammaInd]) << " ";
						}else{
							*(outFiles[gammaInd]) << endl;
						}
					}
				}
			}


			anyUpdates = proposalParams.rhoAnyUpdates();
			if(anyUpdates){
				for(unsigned int j=0;j<nCovariates;j++){
					*(outFiles[rhoOmegaPropInd]) << sampler.proposalParams().rhoAcceptRate(j) <<
						" " << sampler.proposalParams().rhoStdDev(j);
					if(j<(nCovariates-1)){
						*(outFiles[rhoOmegaPropInd]) << " ";
					}else{
						*(outFiles[rhoOmegaPropInd]) << endl;
					}
				}
				proposalParams.rhoAnyUpdates(false);

			}
		}


	}

}

string storeLogFileData(const pReMiuMOptions& options,
								const pReMiuMData& dataset,
								const pReMiuMHyperParams& hyperParams,
								const unsigned int& nClusInit,
								const unsigned int& maxNClusters,
								const double& timeInSecs){

	ostringstream tmpStr;
	tmpStr << "Number of subjects: " << dataset.nSubjects() << endl;
	tmpStr << "Number of prediction subjects: " << dataset.nPredictSubjects() << endl;
	tmpStr << "Prediction type: " << options.predictType() << endl;
	tmpStr << "Sampler type: " << options.samplerType();
	if(options.samplerType().compare("Truncated")==0){
		tmpStr << " " << maxNClusters << " clusters" << endl;
	}else{
		tmpStr << endl;
	}
	tmpStr << "Number of initial clusters: " << nClusInit;
	if(options.nClusInit()==0){
		tmpStr << " (Random, Unif[50,60])" << endl;
	}else{
		tmpStr << endl;
	}
	tmpStr << "Covariates: " << endl;
	if(options.covariateType().compare("Mixed")==0){
		tmpStr << "Number of discrete covariates: " << dataset.nDiscreteCovs() << endl;
		tmpStr << "Number of continuous covariates: " << dataset.nContinuousCovs() << endl;
	}
	for(unsigned int j=0;j<dataset.nCovariates();j++){
		tmpStr << "\t" << dataset.covariateNames(j);
		if(options.covariateType().compare("Discrete")==0){
			tmpStr << " (categorical)";
		} else if(options.covariateType().compare("Mixed")==0){
			if (j < dataset.nDiscreteCovs()){
				tmpStr << " (categorical)";
			}
		}
		tmpStr << endl;
	}
	if(dataset.nFixedEffects()>0){
		tmpStr << "FixedEffects: " << endl;
		for(unsigned int j=0;j<dataset.nFixedEffects();j++){
			tmpStr << "\t" << dataset.fixedEffectNames(j) << endl;
		}
	}else{
		tmpStr<< "No fixed effects" << endl;
	}
	if(dataset.nCategoriesY()>1){
		tmpStr << "NumberOfCategoriesY: " <<  dataset.nCategoriesY() << endl;
	}
	tmpStr << "Model for Y: " << options.outcomeType() << endl;
	if(options.responseExtraVar()){
		tmpStr << "Extra Y variance: True" << endl;
	}else{
		tmpStr << "Extra Y variance: False" << endl;
	}
	if(options.includeResponse()){
		tmpStr << "Include response: True" << endl;
	}else{
		tmpStr << "Include response: False" << endl;
	}
	if(options.fixedAlpha()<=-1){
		tmpStr << "Update alpha: True" << endl;
	}else{
		tmpStr << "Update alpha: False" << endl;
		tmpStr << "Fixed alpha: " << options.fixedAlpha() << endl;
		if(options.dPitmanYor()==0) {
			tmpStr << "Dirichlet process prior, so dPitmanYor: " << options.dPitmanYor() << endl;
		} else {
			tmpStr << "dPitmanYor: " << options.dPitmanYor() << endl;
		}
	}
	if(options.computeEntropy()){
		tmpStr << "Compute allocation entropy: True" << endl;
	}else{
		tmpStr << "Compute allocation entropy: False" << endl;
	}

	tmpStr << "Model for X: " << options.covariateType() << endl;
	tmpStr << "Variable selection: " << options.varSelectType() << endl;

	tmpStr << endl << "Hyperparameters:" << endl;
	if(options.fixedAlpha()<=-1){
		tmpStr << "shapeAlpha: " << hyperParams.shapeAlpha() << endl;
		tmpStr << "rateAlpha: " << hyperParams.rateAlpha() << endl;
	}
	if(options.covariateType().compare("Discrete")==0 ||options.covariateType().compare("Mixed")==0 ){
		tmpStr << "aPhi[j]: ";
		if(options.covariateType().compare("Discrete")==0){
			for(unsigned int j=0;j<dataset.nCovariates();j++){
				tmpStr << hyperParams.aPhi(j) << " ";
			}
		}
		if(options.covariateType().compare("Mixed")==0){
			for(unsigned int j=0;j<dataset.nDiscreteCovs();j++){
				tmpStr << hyperParams.aPhi(j) << " ";
			}
		}
		tmpStr << endl;
	}

	if(options.covariateType().compare("Normal")==0 ||options.covariateType().compare("Mixed")==0 ){
		tmpStr << "mu0: " << endl;
		tmpStr << hyperParams.mu0() << endl;
		tmpStr << "Tau0:" << endl;
		if (options.useIndependentNormal()) {
			tmpStr << hyperParams.Tau0_Indep() << endl;
		}
		else {
			tmpStr << hyperParams.Tau0() << endl;
		}
		tmpStr << "R0: " << endl;
		if (options.useIndependentNormal()) {
			tmpStr << hyperParams.R0_Indep() << endl;
		}
		else {
			tmpStr << hyperParams.R0() << endl;
		}
		tmpStr << "kappa0: " << hyperParams.kappa0() << endl;
		if (options.useHyperpriorR1()|| options.useIndependentNormal()) tmpStr << "kappa1: " << hyperParams.kappa1() << endl;
		tmpStr << "nu0: " << hyperParams.nu0() << endl;
	}
	

	tmpStr << "muTheta: " << hyperParams.muTheta() << endl;
	tmpStr << "sigmaTheta: " << hyperParams.sigmaTheta() << endl;
	tmpStr << "dofTheta: " << hyperParams.dofTheta() << endl;

	if(dataset.nFixedEffects()>0){
		tmpStr << "muBeta: " << hyperParams.muBeta() << endl;
		tmpStr << "sigmaBeta: " << hyperParams.sigmaBeta() << endl;
		tmpStr << "dofBeta: " << hyperParams.dofBeta() << endl;
	}

	if(options.responseExtraVar()){
		tmpStr << "shapetauEpsilon: " << hyperParams.shapeTauEpsilon() << endl;
		tmpStr << "ratetauEpsilon: " << hyperParams.rateTauEpsilon() << endl;
	}

	if(options.varSelectType().compare("None")!=0){
		tmpStr << "aRho: " << hyperParams.aRho() << endl;
		tmpStr << "bRho: " << hyperParams.bRho() << endl;
		tmpStr << "atomRho: " << hyperParams.atomRho() << endl;
	}

	if(dataset.outcomeType().compare("Normal")==0){
		tmpStr << "shapeSigmaSqY: " << hyperParams.shapeSigmaSqY() << endl;
		tmpStr << "scaleSigmaSqY: " << hyperParams.scaleSigmaSqY() << endl;
	}

	if(dataset.outcomeType().compare("Quantile")==0){
		tmpStr << "shapeSigmaSqY: " << hyperParams.shapeSigmaSqY() << endl;
		tmpStr << "scaleSigmaSqY: " << hyperParams.scaleSigmaSqY() << endl;
		tmpStr << "pQuantile: " << hyperParams.pQuantile() << endl;
	}

	if(dataset.outcomeType().compare("Survival")==0){
		tmpStr << "Weibull with fixed shape parameter: " << options.weibullFixedShape() << endl;
		tmpStr << "shapeNu: " << hyperParams.shapeNu() << endl;
		tmpStr << "scaleNu: " << hyperParams.scaleNu() << endl;
	}

	if(options.includeCAR()){
		tmpStr << "shapeTauCAR: " << endl;
		tmpStr << hyperParams.shapeTauCAR() << endl;
		tmpStr << "rateTauCAR:" << endl;
		tmpStr << hyperParams.rateTauCAR() << endl;
	}

	tmpStr << endl << options.nSweeps()+options.nBurn() << " sweeps done in " <<
				timeInSecs << " seconds" << endl;

	return tmpStr.str();
}

#endif // DIPBACIO_H_
