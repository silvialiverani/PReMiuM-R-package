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

#include<Eigen/Core>
#include<Eigen/Cholesky>
#include<Eigen/LU>

// Custom includes
#include<MCMC/chain.h>
#include<MCMC/model.h>
#include<MCMC/sampler.h>
#include<PReMiuMOptions.h>
#include<PReMiuMData.h>
#include<PReMiuMModel.h>

using namespace Eigen;

using std::vector;
using std::cout;
using std::endl;
using std::ostringstream;
using std::istringstream;
using std::string;


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
			cout << endl << "### PReMiuMpp Help Page. ###" << endl;
			cout << endl;
			cout << "Possible arguments (defaults in parentheses):" << endl << endl;
			cout << "--help" << endl << "\tShow this help page" << endl;
			cout << "--input=<string>"<< endl << "\tThe full file path to the input data (./input.txt)" << endl;
			cout << "--output=<string>" << endl << "\tThe file stem (including full file path) where the data" << endl << "\tshould be written (./output)" << endl;
			cout << "--hyper=<string>" << endl << "\tThe full file path to the file containing hyper" << endl << "\t parameters. (Hyper parameter file not used)" << endl;
			cout << "--predict=<string>" << endl << "\tThe full file path to the file containing prediction" << endl << "\t covariates. (Prediction file not used)" << endl;
			cout << "--nSweeps=<unsigned int>" << endl << "\tThe number of sweeps (after burn in) to run the" << endl << "\t sampler for (10000)" << endl;
			cout << "--nBurn=<unsigned int>" << endl << "\tThe number of sweeps in the burn in period (1000)" << endl;
			cout << "--reportBurnIn=<bool>" << endl << "\tIt enables reporting in the output files of the burn-in period (true)." << endl;
			cout << "--nProgress=<unsigned int>" << endl << "\tThe number of sweeps at which to print a" << endl << "progress update (500)" << endl;
			cout << "--nFilter=<unsigned int>" << endl << "\tThe frequency (in sweeps) with which to write" << endl << "\tthe output to file (1)" << endl;
			cout << "--nClusInit=<unsigned int>" << endl << "\tThe number of clusters individuals should be" << endl << "\tinitially randomly assigned to (Unif[50,60])" << endl;
			cout << "--seed=<unsigned int>" << endl << "\tThe value for the seed for the random number" << endl << "\tgenerator (current time)" << endl;
			cout << "--yModel=<string>" << endl << "\tThe model type for the outcome variable. Options are" << endl << "\tcurrently 'Bernoulli','Poisson','Binomial', 'Categorical' and 'Normal' (Bernoulli)" << endl;
			cout << "--xModel=<string>" << endl << "\tThe model type for the covariates. Options are" << endl << "\tcurrently 'Discrete', 'Normal' and 'Mixed' (Discrete)" << endl;
			cout << "--sampler=<string>" << endl << "\tThe sampler type to be used. Options are" << endl << "\tcurrently 'SliceDependent', 'SliceIndependent' and 'Truncated' (SliceDependent)" << endl;
			cout << "--alpha=<double>" << endl << "\tThe value to be used if alpha is to remain fixed." << endl << "\tIf a negative value is used then alpha is updated (-1)" << endl;
			cout << "--excludeY" << endl << "\tIf included only the covariate data X is modelled (not included)" << endl;
			cout << "--extraYVar" << endl << "\tIf included extra Gaussian variance is included in the" << endl << "\tresponse model (not included)." << endl;
			cout << "--varSelect=<string>" << endl << "\tThe type of variable selection to be used 'None'," << endl << "\t'BinaryCluster' or 'Continuous' (None)" << endl;
			cout << "--entropy" << endl << "\tIf included then we compute allocation entropy (not included)" << endl;
			exit(0);
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
					long rndSeed=(long)atoi(tmpStr.c_str());
					options.seed(rndSeed);
				}else if(inString.find("--yModel")!=string::npos){
					size_t pos = inString.find("=")+1;
					string outcomeType = inString.substr(pos,inString.size()-pos);
					if(outcomeType.compare("Poisson")!=0&&outcomeType.compare("Bernoulli")!=0&&
							outcomeType.compare("Categorical")!=0&&
							outcomeType.compare("Binomial")!=0&&outcomeType.compare("Normal")!=0){
						// Illegal outcome model entered
						wasError=true;
						break;
					}
					options.outcomeType(outcomeType);
					if(outcomeType.compare("Normal")==0&&options.responseExtraVar()){
						cout << "Response extra variation not permitted with Normal response" << endl;
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
				}else if(inString.find("--sampler")!=string::npos){
					size_t pos = inString.find("=")+1;
					string samplerType = inString.substr(pos,inString.size()-pos);
					if(samplerType.compare("SliceDependent")!=0&&samplerType.compare("SliceIndependent")!=0
							&&samplerType.compare("Truncated")!=0){
						// Illegal covariate type entered
						wasError=true;
						break;
					}
					options.samplerType(samplerType);
				}else if(inString.find("--alpha")!=string::npos){
					size_t pos = inString.find("=")+1;
					string tmpStr = inString.substr(pos,inString.size()-pos);
					double alpha=(double)atof(tmpStr.c_str());
					options.fixedAlpha(alpha);
				}else if(inString.find("--excludeY")!=string::npos){
					options.includeResponse(false);
				}else if(inString.find("--extraYVar")!=string::npos){
					if(options.outcomeType().compare("Normal")!=0){
						options.responseExtraVar(true);
					}else{
						cout << "Response extra variation not permitted with Normal response" << endl;
					}
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
				}else{
					cout << "Unknown command line option." << endl;
					wasError=true;
					break;
				}
				currArg++;
			}
		}
	}

	// Return if there was an error
	if(wasError){
		cout << "Please use:" << endl;
		cout << "\t profileRegression --help" << endl;
		cout << "to get help on correct usage." << endl;
		exit(-1);
	}

	return options;

}

// Read the PReMiuM data set
void importPReMiuMData(const string& fitFilename,const string& predictFilename,pReMiuMData& dataset){

	ifstream inputFile,predictFile;
	inputFile.open(fitFilename.c_str());
	if(!inputFile.is_open()){
		cout << "Input file not found" << endl;
		exit(-1);
	}
	if(predictFilename.compare("")!=0){
		predictFile.open(predictFilename.c_str());
		if(!predictFile.is_open()){
			cout << "Prediction covariate file not found" << endl;
			exit(-1);
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
				cout << "Illegal number of covariates, discrete covariates or continuous covariates" <<endl;
				// Illegal number of covariates, discrete covariates or continuous covariates
				wasError=true;
		}
		if(nDiscreteCovs==0 || nContinuousCovs==0){
				cout << "If xModel=Mixed a positive number of discrete and continuous covariates must be provided " <<endl;
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
	missingX.resize(nSubjects+nPredictSubjects);
	nContinuousCovariatesNotMissing.resize(nSubjects+nPredictSubjects);
	vector<double> meanX(nCovariates,0);
	vector<unsigned int> nXNotMissing(nCovariates,0);
	for(unsigned int i=0;i<nSubjects;i++){
		if(outcomeType.compare("Normal")==0){
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
						continuousX[i][j-nDiscreteCovs]=meanX[j];
					}
				}
			}
		}
	}

	inputFile.close();
	if(predictFile.is_open()){
		predictFile.close();
	}

	// Return if there was an error
	if(wasError){
		cout << "Please use:" << endl;
		cout << "\t profileRegression --help" << endl;
		cout << "to get help on correct usage." << endl;
		exit(-1);
	}

}

// Function to read the hyper parameters from file
void readHyperParamsFromFile(const string& filename,pReMiuMHyperParams& hyperParams){

	ifstream inputFile;
	inputFile.open(filename.c_str());
	if(!inputFile.is_open()){
		cout << "Parameter file not found" << endl;
		exit(-1);
	}

	string inString;

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
		}else if(inString.find("useReciprocalNCatsPhi")==0){
			size_t pos = inString.find("=")+1;
			string tmpStr = inString.substr(pos,inString.size()-pos);
			bool useRecip = false;
			if(tmpStr.compare("true")==0){
				useRecip = true;
			}
			hyperParams.useReciprocalNCatsPhi(useRecip);
		}else if(inString.find("aPhi")==0){
			size_t pos = inString.find("=")+1;
			string tmpStr = inString.substr(pos,inString.size()-pos);
			vector<double> aVec;
			while(tmpStr.find(" ")!=string::npos){
				pos = tmpStr.find(" ");
				if(pos==(tmpStr.size()-1)){
					break;
				}
				string elem = tmpStr.substr(0,pos-1);
				aVec.push_back((double)atof(elem.c_str()));
				tmpStr = tmpStr.substr(pos+1,tmpStr.size()-pos-1);
			}
			hyperParams.aPhi(aVec);
		}else if(inString.find("mu0")==0){
			size_t pos = inString.find("=")+1;
			string tmpStr = inString.substr(pos,inString.size()-pos);
			vector<double> muVec;
			while(tmpStr.find(" ")!=string::npos){
				pos = tmpStr.find(" ");
				if(pos==(tmpStr.size()-1)){
					break;
				}
				string elem = tmpStr.substr(0,pos-1);
				muVec.push_back((double)atof(elem.c_str()));
				tmpStr = tmpStr.substr(pos+1,tmpStr.size()-pos-1);
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
					break;
				}
				string elem = tmpStr.substr(0,pos-1);
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
		}else if(inString.find("R0")==0){
			size_t pos = inString.find("=")+1;
			string tmpStr = inString.substr(pos,inString.size()-pos);
			vector<double> RVec;
			while(tmpStr.find(" ")!=string::npos){
				pos = tmpStr.find(" ");
				if(pos==(tmpStr.size()-1)){
					break;
				}
				string elem = tmpStr.substr(0,pos-1);
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
		}else if(inString.find("kappa0")==0){
			size_t pos = inString.find("=")+1;
			string tmpStr = inString.substr(pos,inString.size()-pos);
			unsigned int kappa0 = (unsigned int)atoi(tmpStr.c_str());
			hyperParams.kappa0(kappa0);
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
		}
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

	vector<unsigned int> nCategories;
	nCategories = dataset.nCategories();

	// Set the hyper parameters to their default values
	hyperParams.setSizes(nCovariates,nDiscreteCovs,
			nContinuousCovs,covariateType);
	hyperParams.setDefaults(dataset,options);
	// Read the parameters from file if file provided
	if(hyperParamFileName.compare("")!=0){
		readHyperParamsFromFile(hyperParamFileName,hyperParams);
	}

	// Allocate the right sizes for each of the parameter variables
	// This also switches "on" all variable indicators (gamma)
	// This gets changed below if variable selection is being done
	params.setSizes(nSubjects,nCovariates,nDiscreteCovs,nContinuousCovs,nFixedEffects,nCategoriesY,nPredictSubjects,nCategories,nClusInit,covariateType);
	unsigned int maxNClusters=params.maxNClusters();

	// Fix the number of clusters if we are using the truncated sampler
	if(samplerType.compare("Truncated")==0){
		maxNClusters=20;
		if((nClusInit+10)>maxNClusters){
			maxNClusters=nClusInit+10;
		}
		// Now compute the bound recommended in Ishwaran and James 2001
		double multiplier=0.0;
		if(options.fixedAlpha()>0){
			multiplier=options.fixedAlpha();
		}else{
			// Use the expected value of alpha as the multiplier
			multiplier=hyperParams.shapeAlpha()/hyperParams.rateAlpha();
		}
		double computedBound=1+multiplier*(log(4.0*nSubjects)-log(hyperParams.truncationEps()));
		if(computedBound>maxNClusters){
			maxNClusters=computedBound;
		}
		params.maxNClusters(maxNClusters,covariateType);
	}

	// Copy the dataset X matrix to a working object in params
	params.workDiscreteX(dataset.discreteX());
	params.workContinuousX(dataset.continuousX());

	// Now initialise the actual parameters
	// Define a uniform random number generator
	randomUniform unifRand(0,1);

	randomGamma gammaRand(hyperParams.shapeAlpha(),1.0/hyperParams.rateAlpha());

	double alpha=gammaRand(rndGenerator);
	if(options.fixedAlpha()>0){
		alpha=options.fixedAlpha();
	}
	params.alpha(alpha);

	vector<unsigned int> nXInCluster(maxNClusters,0);
	unsigned int maxZ=0;
	if(nClusInit==0){
		nClusInit=50+(unsigned int)11*unifRand(rndGenerator);
	}
	params.workNClusInit(nClusInit);
	for(unsigned int i=0;i<nSubjects+nPredictSubjects;i++){
		int c=(int) nClusInit*unifRand(rndGenerator);
		params.z(i,c,covariateType);
		if(c>(int)maxZ){
			maxZ=c;
		}
		if(i<nSubjects){
			nXInCluster[c]++;
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
		double vVal = betaRand(rndGenerator,1.0+params.workNXInCluster(c),alpha+sumCPlus1ToMaxMembers[c]);
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
			double v=betaRand(rndGenerator,1.0,alpha);
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
			if(samplerType.compare("SliceDependent")){
				ui = exp(params.logPsi(zi))*unifRand(rndGenerator);
			}else if(samplerType.compare("SliceIndependent")){
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
				double v=betaRand(rndGenerator,1.0,alpha);
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

		params.maxNClusters(maxNClusters,covariateType);
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
		normal_distribution<double> norm01(0.0,1.0);
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
			MatrixXd covMat(nCovariates,nCovariates);
			covMat = (hyperParams.Tau0()+params.workNXInCluster(c)*hyperParams.Tau0()).inverse();
			VectorXd meanVec(nCovariates);
			meanVec = hyperParams.Tau0()*hyperParams.mu0()+params.workNXInCluster(c)*hyperParams.Tau0()*meanX[c];
			meanVec = covMat*meanVec;

			VectorXd mu(nCovariates);
			// We sample from this posterior
			mu = multivarNormalRand(rndGenerator,meanVec,covMat);

			// We store our sample
			params.mu(c,mu);

		}

		// Now we can sample Tau_c for each cluster
		vector<MatrixXd> Rc(maxNClusters);
		for(unsigned int c=0;c<maxNClusters;c++){
			Rc[c].setZero(nCovariates,nCovariates);
		}

		for(unsigned int i=0;i<nSubjects;i++){
			unsigned int zi = params.z(i);
			Rc[zi]=Rc[zi]+(xi[i]-params.mu(zi))*((xi[i]-params.mu(zi)).transpose());
		}

		for(unsigned int c=0;c<maxNClusters;c++){
			Rc[c]=(hyperParams.R0().inverse()+Rc[c]).inverse();
			MatrixXd Tau = wishartRand(rndGenerator,Rc[c],params.workNXInCluster(c)+hyperParams.kappa0());
			params.Tau(c,Tau);
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
			params.nullMu(nullMu);
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
		normal_distribution<double> norm01(0.0,1.0);
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
			MatrixXd covMat(nContinuousCovs,nContinuousCovs);
			covMat = (hyperParams.Tau0()+params.workNXInCluster(c)*hyperParams.Tau0()).inverse();
			VectorXd meanVec(nContinuousCovs);
			meanVec = hyperParams.Tau0()*hyperParams.mu0()+params.workNXInCluster(c)*hyperParams.Tau0()*meanX[c];
			meanVec = covMat*meanVec;

			VectorXd mu(nContinuousCovs);
			// We sample from this posterior
			mu = multivarNormalRand(rndGenerator,meanVec,covMat);

			// We store our sample
			params.mu(c,mu);

		}

		// Now we can sample Tau_c for each cluster
		vector<MatrixXd> Rc(maxNClusters);
		for(unsigned int c=0;c<maxNClusters;c++){
			Rc[c].setZero(nContinuousCovs,nContinuousCovs);
		}

		for(unsigned int i=0;i<nSubjects;i++){
			unsigned int zi = params.z(i);
			Rc[zi]=Rc[zi]+(xi[i]-params.mu(zi))*((xi[i]-params.mu(zi)).transpose());
		}

		for(unsigned int c=0;c<maxNClusters;c++){
			Rc[c]=(hyperParams.R0().inverse()+Rc[c]).inverse();
			MatrixXd Tau = wishartRand(rndGenerator,Rc[c],params.workNXInCluster(c)+hyperParams.kappa0());
			params.Tau(c,Tau);
		}

		// Now do the null mu for variable selection
		// In all cases, initialise it at the value it will be fixed at for
		// the continuous indicator case
		if(varSelectType.compare("None")!=0){
			vector<double> meanXVec(nContinuousCovs,0.0);
			vector<unsigned int> countXVec(nContinuousCovs,0);
			for(unsigned int i=0;i<nSubjects;i++){
				for(unsigned int j=0;j<nContinuousCovs;j++){
					if(!dataset.missingX(i,j)){
						meanXVec[j]+=dataset.continuousX(i,j);
						countXVec[j]+=1;
					}
				}
			}
			VectorXd nullMu=VectorXd::Zero(nContinuousCovs);
			for(unsigned int j=0;j<nContinuousCovs;j++){
				nullMu(j)=meanXVec[j]/(double)countXVec[j];
			}
			params.nullMu(nullMu);
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
			if(unifRand(rndGenerator)<0.01){
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
			params.rho(j,rho[j],covariateType,varSelectType);
			if(varSelectType.compare("BinaryCluster")==0){
				for(unsigned int c=0;c<maxNClusters;c++){
					params.gamma(c,j,gamma[c][j],covariateType);
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

		if(outcomeType.compare("Normal")==0){
			randomGamma gammaRand(hyperParams.shapeSigmaSqY(),1.0/hyperParams.scaleSigmaSqY());
			double sigmaSqY=1.0/(gammaRand(rndGenerator));
			params.sigmaSqY(sigmaSqY);

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
		bool responseExtraVar = sampler.model().options().responseExtraVar();
		double fixedAlpha = sampler.model().options().fixedAlpha();
		string outcomeType = sampler.model().options().outcomeType();
		bool computeEntropy = sampler.model().options().computeEntropy();
		unsigned int nFixedEffects = params.nFixedEffects(outcomeType);
		string varSelectType = sampler.model().options().varSelectType();

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
			}else if(covariateType.compare("Mixed")==0){
				fileName = fileStem + "_phi.txt";
				outFiles.push_back(new ofstream(fileName.c_str()));
				fileName = fileStem + "_mu.txt";
				outFiles.push_back(new ofstream(fileName.c_str()));
				fileName = fileStem + "_Sigma.txt";
				outFiles.push_back(new ofstream(fileName.c_str()));
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
			if(fixedAlpha<0){
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
				if(outcomeType.compare("Normal")==0){
					fileName = fileStem + "_sigmaSqY.txt";
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
		int nClustersInd=-1,psiInd=-1,phiInd=-1,muInd=-1,SigmaInd=-1,zInd=-1,entropyInd=-1,alphaInd=-1;
		int logPostInd=-1,nMembersInd=-1,alphaPropInd;
		int thetaInd=-1,betaInd=-1,thetaPropInd=-1,betaPropInd=-1,sigmaSqYInd=-1,epsilonInd=-1;
		int sigmaEpsilonInd=-1,epsilonPropInd=-1,omegaInd=-1,rhoInd=-1;
		int rhoOmegaPropInd=-1,gammaInd=-1,nullPhiInd=-1,nullMuInd=-1;
		int predictThetaRaoBlackwellInd=-1;

		int r=0;
		nClustersInd=r++;
		psiInd=r++;
		if(covariateType.compare("Discrete")==0){
			phiInd=r++;
		}else if(covariateType.compare("Normal")==0){
			muInd=r++;
			SigmaInd=r++;
		}else if(covariateType.compare("Mixed")==0){
			phiInd=r++;
			muInd=r++;
			SigmaInd=r++;
		}
		zInd=r++;
		entropyInd=r++;
		alphaInd=r++;
		logPostInd=r++;
		nMembersInd=r++;
		if(fixedAlpha<0){
			alphaPropInd=r++;
		}

		if(includeResponse){
			thetaInd=r++;
			betaInd=r++;
			thetaPropInd=r++;
			betaPropInd=r++;
			if(outcomeType.compare("Normal")==0){
				sigmaSqYInd=r++;
			}
			if(responseExtraVar){
				epsilonInd=r++;
				sigmaEpsilonInd=r++;
				epsilonPropInd=r++;
			}
			if(nPredictSubjects>0){
				predictThetaRaoBlackwellInd=r++;
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

			for(unsigned int j1=0;j1<nCovariates;j1++){
				for(unsigned int j2=0;j2<nCovariates;j2++){
					for(unsigned int c=0;c<maxNClusters;c++){
						*(outFiles[SigmaInd]) << params.Sigma(c,j1,j2);
						if(c<(maxNClusters-1)||j1<(nCovariates-1)||j2<(nCovariates-1)){
							*(outFiles[SigmaInd]) << " ";
						}
					}
				}
			}
			*(outFiles[SigmaInd]) << endl;
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

			for(unsigned int j1=0;j1<nContinuousCovs;j1++){
				for(unsigned int j2=0;j2<nContinuousCovs;j2++){
					for(unsigned int c=0;c<maxNClusters;c++){
						*(outFiles[SigmaInd]) << params.Sigma(c,j1,j2);
						if(c<(maxNClusters-1)||j1<(nContinuousCovs-1)||j2<(nContinuousCovs-1)){
							*(outFiles[SigmaInd]) << " ";
						}
					}
				}
			}
			*(outFiles[SigmaInd]) << endl;
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
				if(outcomeType.compare("Normal")==0){
					*(outFiles[sigmaSqYInd]) << params.sigmaSqY() << endl;
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
		if(fixedAlpha<0||sweep==0){
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
		if(fixedAlpha<0){
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
				if(sweep==0){
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
	tmpStr << "Sampler type: " << options.samplerType();
	if(options.samplerType().compare("Truncated")==0){
		tmpStr << " " << maxNClusters << " clusters" << endl;
	}else{
		tmpStr << endl;
	}
	tmpStr << "Number of initial clusters: " << nClusInit;
	if(options.nClusInit()==0){
		tmpStr << " (Random, Unif[5,15])" << endl;
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
	if(options.fixedAlpha()<0){
		tmpStr << "Update alpha: True" << endl;
	}else{
		tmpStr << "Update alpha: False" << endl;
		tmpStr << "Fixed alpha: " << options.fixedAlpha() << endl;
	}
	if(options.computeEntropy()){
		tmpStr << "Compute allocation entropy: True" << endl;
	}else{
		tmpStr << "Compute allocation entropy: False" << endl;
	}

	tmpStr << "Model for X: " << options.covariateType() << endl;
	tmpStr << "Variable selection: " << options.varSelectType() << endl;

	tmpStr << endl << "Hyperparameters:" << endl;
	if(options.fixedAlpha()<0){
		tmpStr << "shapeAlpha: " << hyperParams.shapeAlpha() << endl;
		tmpStr << "rateAlpha: " << hyperParams.rateAlpha() << endl;
	}
	if(options.covariateType().compare("Discrete")==0 ||options.covariateType().compare("Mixed")==0 ){
		if(hyperParams.useReciprocalNCatsPhi()){
			tmpStr << "aPhi[j]: 1/nCategories[j]" << endl;
		}else{
			tmpStr << "aPhi[j]: ";
			for(unsigned int j=0;j<dataset.nCovariates();j++){
				tmpStr << hyperParams.aPhi(j) << " ";
			}
			tmpStr << endl;
		}

	}

	if(options.covariateType().compare("Normal")==0 ||options.covariateType().compare("Mixed")==0 ){
		tmpStr << "mu0: " << endl;
		tmpStr << hyperParams.mu0() << endl;
		tmpStr << "Tau0:" << endl;
		tmpStr << hyperParams.Tau0() << endl;
		tmpStr << "R0: "  << endl;
		tmpStr << hyperParams.R0() << endl;
		tmpStr << "kappa0: " << hyperParams.kappa0() << endl;
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
		tmpStr << "bRho: " << hyperParams.aRho() << endl;
	}

	if(dataset.outcomeType().compare("Normal")==0){
		tmpStr << "shapeSigmaSqY: " << hyperParams.shapeSigmaSqY() << endl;
		tmpStr << "scaleSigmaSqY: " << hyperParams.scaleSigmaSqY() << endl;
	}

	tmpStr << endl << options.nSweeps()+options.nBurn() << " sweeps done in " <<
				timeInSecs << " seconds" << endl;

	return tmpStr.str();
}

#endif // DIPBACIO_H_
