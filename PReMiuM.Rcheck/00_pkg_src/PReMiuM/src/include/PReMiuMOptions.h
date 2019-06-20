/// \file PReMiuMOptions.h
/// \author David Hastie
/// \brief Header file for options class for PReMiuMpp

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


#ifndef DIPBACOPTIONS_H_
#define DIPBACOPTIONS_H_

// Standard includes
#include<vector>
#include<iostream>
#include<fstream>
#include<string>
#include<numeric>
#include<ctime>
#include<cstdint>

using std::vector;
using std::ifstream;
using std::ofstream;
using std::string;
using std::time;

/// \class pReMiuMOptions PReMiuMOptions.h "PReMiuMOptions.h"
/// \brief A class for PReMiuM options
class pReMiuMOptions{

	public:
		/// \brief Default constructor
		pReMiuMOptions() {
			/* ---------- Default Values -------------------*/
			// Default values for the command line arguments
			// File names
			_inFileName = "input.txt";
			_outFileStem  = "output";
			_hyperParamFileName = "";
			_predictFileName = "";

			// MCMC sampler variables
			_nSweeps=1000;
			_nBurn=1000;
			_reportBurnIn=false;
			_nProgress=500;
			_nFilter=1;
			// The default of 0 initial clusters, means a random number between
			// 50 and 60 is used at initialisation
			_nClusInit=0;

			// Random number seed
			_seed=(uint_fast32_t) time(0);

			// Profile regression variables
			_outcomeType="Bernoulli";
			_covariateType="Discrete";
			_includeResponse = true;
			_whichLabelSwitch = "123";
			_responseExtraVar = false;
			_doPrediction = false;
			_varSelectType ="None";
			_fixedAlpha=-2;
			_dPitmanYor=0;
			_samplerType="SliceDependent";
			_computeEntropy=false;
			_includeCAR=false;
			_includeuCARinit=false;
			_neighbourFileName="Neighbour.txt";
			_uCARinitFileName="uCARinit.txt";
			_PoissonCARadaptive=false;
			_predictType ="RaoBlackwell";
			_weibullFixedShape=false;
			_useNormInvWishPrior=false;
			_useHyperpriorR1=false;
			_useIndependentNormal = false;
			_useSeparationPrior = false;
		};

		/// \brief Default destructor
		~pReMiuMOptions(){};

		/// \brief Return the number of sweeps
		unsigned int nSweeps() const{
			return _nSweeps;
		}

		/// \brief Set the number of sweeps
		void nSweeps(const unsigned int& nSw){
			_nSweeps=nSw;
		}

		/// \brief Return the number of burn in sweeps
		unsigned int nBurn() const{
			return _nBurn;
		}

		/// \brief Set the number of burn in sweeps
		void nBurn(const unsigned int& nB){
			_nBurn=nB;
		}

		/// \brief Return whether the burn in period is reported
		bool reportBurnIn() const{
			return _reportBurnIn;
		}

		/// \brief Sets whether the burn in period is reported
		void reportBurnIn(const bool& rBI){
			_reportBurnIn=rBI;
		}

		/// \brief Return the number of filter sweeps
		unsigned int nFilter() const{
			return _nFilter;
		}

		/// \brief Set the number of filter sweeps
		void nFilter(const unsigned int& nF){
			_nFilter=nF;
		}

		/// \brief Return the number of progress iterations
		unsigned int nProgress() const{
			return _nProgress;
		}

		/// \brief Set the number of sweeps
		void nProgress(const unsigned int& nPr){
			_nProgress=nPr;
		}

		/// \brief Return the number of progress iterations
		unsigned int nClusInit() const{
			return _nClusInit;
		}

		/// \brief Set the number of sweeps
		void nClusInit(const unsigned int& nCl){
			_nClusInit=nCl;
		}

		/// \brief Return the random number seed
		uint_fast32_t seed() const{
			return _seed;
		}

		/// \brief Set the random number seed
		void seed(const uint_fast32_t& rndSeed){
			uint_fast32_t rndSeedNew = rndSeed;
			if(rndSeedNew==0){
				rndSeedNew=(uint_fast32_t) time(0);
			}
			_seed=rndSeedNew;
		}

		/// \brief Return the input file name
		string inFileName() const{
			return _inFileName;
		}

		/// \brief Set the input file name
		void inFileName(const string& inFile){
			_inFileName=inFile;
		}

		/// \brief Return the output file stem
		string outFileStem() const{
			return _outFileStem;
		}

		/// \brief Set the output file stem
		void outFileStem(const string& outStem){
			_outFileStem=outStem;
		}

		/// \brief Return the input file name
		string hyperParamFileName() const{
			return _hyperParamFileName;
		}

		/// \brief Set the input file name
		void hyperParamFileName(const string& hyperFile){
			_hyperParamFileName=hyperFile;
		}

		/// \brief Return the prediction file name
		string predictFileName() const{
			return _predictFileName;
		}

		/// \brief Set the prediction file name
		void predictFileName(const string& predictFile){
			_predictFileName=predictFile;
			_doPrediction=true;
		}

		/// \brief Return the outcome type
		string outcomeType() const{
			return _outcomeType;
		}

		/// \brief Set the outcome type
		void outcomeType(const string& outType){
			_outcomeType=outType;
		}

		/// \brief Return the covariate type
		string covariateType() const{
			return _covariateType;
		}

		/// \brief Set the covariate type
		void covariateType(const string& covType){
			_covariateType=covType;
		}


		/// \brief Return whether we are including response
		bool includeResponse() const{
			return _includeResponse;
		}

		/// \brief Set whether we are including the response
		void includeResponse(const bool& incResp){
			_includeResponse=incResp;
		}

		/// \brief Return which label switch moves are implemented
		string whichLabelSwitch() const{
			return _whichLabelSwitch;
		}

		/// \brief Set which label switch moves are implemented
		void whichLabelSwitch(const string& ls){
			_whichLabelSwitch=ls;
		}

		/// \brief Return the fixed value of alpha (negative value means not fixed)
		double fixedAlpha() const{
			return _fixedAlpha;
		}

		/// \brief Return the fixed value of alpha (negative value means not fixed)
		void fixedAlpha(const double& alphaVal){
			_fixedAlpha=alphaVal;
		}

		/// \brief Return the fixed value of the discount parameter of the Pitman Yor process prior
		double dPitmanYor() const{
			return _dPitmanYor;
		}

		/// \brief Return the fixed value of the discount parameter of the Pitman Yor process prior
		void dPitmanYor(const double& dPitmanYorVal){
			_dPitmanYor=dPitmanYorVal;
		}

		/// \brief Return the sampler method
		string samplerType() const{
			return _samplerType;
		}

		/// \brief Set the outcome type
		void samplerType(const string& sampType){
			_samplerType=sampType;
		}

		/// \brief Return whether we are including response
		bool doPrediction() const{
			return _doPrediction;
		}

		/// \brief Return whether there is response extra variation
		bool responseExtraVar() const{
			return _responseExtraVar;
		}

		/// \brief Set whether there is response extra variation
		void responseExtraVar(const bool& respExVar){
			_responseExtraVar=respExVar;
		}

		/// \brief Return the variable selection type
		string varSelectType() const{
			return _varSelectType;
		}

		/// \brief Set the variable selection type
		void varSelectType(const string& varSelType){
			_varSelectType=varSelType;
		}

		/// \brief Return whether we are computing the entropy
		bool computeEntropy() const{
			return _computeEntropy;
		}

		/// \brief Set whether we are computing the entropy
		void computeEntropy(const bool& compEntr){
			_computeEntropy=compEntr;
		}

		/// \brief Return whether we are including CAR random term
		bool includeCAR() const{
			return _includeCAR;
		}

		/// \brief Set whether we are including CAR random term
		void includeCAR(const bool& incCar){
			_includeCAR=incCar;
		}

		/// \brief Return whether we are including uCAR initialisation
		bool includeuCARinit() const{
			return _includeuCARinit;
		}

		/// \brief Set whether we are including uCAR initialisation
		void includeuCARinit(const bool& incCar){
			_includeuCARinit=incCar;
		}

		/// \brief Return the neighbour structure file name
		string neighbourFileName() const{
			return _neighbourFileName;
		}

		/// \brief Set whether we are including CAR random term
		void neighbourFileName(const string& neiFileName){
			_neighbourFileName=neiFileName;
		}

		/// \brief Return the uCAR initialisation
		string uCARinitFileName() const{
			return _uCARinitFileName;
		}

		/// \brief Set whether we are including uCAR initialisation
		void uCARinitFileName(const string& uCARFileName){
			_uCARinitFileName=uCARFileName;
		}

		/// \brief Return whether the adaptive rejection sampler is to be used for Poisson CAR
		bool PoissonCARadaptive() const{
			return _PoissonCARadaptive;
		}

		/// \brief Set whether the adaptive rejection sampler is to be used for Poisson CAR
		void PoissonCARadaptive(const bool& adaptive){
			_PoissonCARadaptive=adaptive;
		}


		/// \brief Return the prediction type
		string predictType() const{
			return _predictType;
		}

		/// \brief Set the prediction type
		void predictType(const string& predType){
			_predictType=predType;
		}

		/// \brief Return whether the shape parameter nu of the weibull for survival outcome is fixed or cluster specific
		bool weibullFixedShape() const{
			return _weibullFixedShape;
		}

		/// \brief Set whether the shape parameter nu of the weibull for survival outcome is fixed or cluster specific
		void weibullFixedShape(const bool& fixedShape){
			_weibullFixedShape=fixedShape;
		}

		/// \brief Return whether we are including response
		bool useNormInvWishPrior() const{
			return _useNormInvWishPrior;
		}

		/// \brief Set whether we are including the response
		void useNormInvWishPrior(const bool& useNIWP){
			_useNormInvWishPrior=useNIWP;
		}


		/// \brief Return whether we are including hyperprior for R1
		bool useHyperpriorR1() const{
			return _useHyperpriorR1;
		}

		/// \brief Set whether we are including hyperprior for R1
		void useHyperpriorR1(const bool& useR1){
			_useHyperpriorR1=useR1;
		}

		/// \brief Return whether we are using the Independent likelihood for continuous variables
		bool useIndependentNormal() const {
			return _useIndependentNormal;
		}

		/// \brief Set whether we are using the Independent likelihood for continuous variables
		void useIndependentNormal(const bool& useIndep) {
			_useIndependentNormal = useIndep;
		}

		bool useSeparationPrior() const {
			return _useSeparationPrior;
		}

		void useSeparationPrior(const bool& useSP) {
			_useSeparationPrior = useSP;
		}


		// Copy operator
		pReMiuMOptions& operator=(const pReMiuMOptions& options){

			_inFileName = options.inFileName();
			_outFileStem = options.outFileStem();
			_hyperParamFileName = options.hyperParamFileName();
			_predictFileName = options.predictFileName();
			_nSweeps = options.nSweeps();
			_nBurn = options.nBurn();
			_reportBurnIn = options.reportBurnIn();
			_nFilter=options.nFilter();
			_nProgress=options.nProgress();
			_nClusInit=options.nClusInit();
			_seed=options.seed();
			_outcomeType=options.outcomeType();
			_covariateType=options.covariateType();
			_includeResponse=options.includeResponse();
			_whichLabelSwitch=options.whichLabelSwitch();
			_fixedAlpha=options.fixedAlpha();
			_dPitmanYor=options.dPitmanYor();
			_samplerType=options.samplerType();
			_doPrediction=options.doPrediction();
			_responseExtraVar=options.responseExtraVar();
			_varSelectType=options.varSelectType();
			_computeEntropy=options.computeEntropy();
			_includeCAR=options.includeCAR();
			_includeuCARinit=options.includeuCARinit();
			_neighbourFileName=options.neighbourFileName();
			_uCARinitFileName=options.uCARinitFileName();
			_PoissonCARadaptive=options.PoissonCARadaptive();
			_predictType=options.predictType();
			_weibullFixedShape=options.weibullFixedShape();
			_useNormInvWishPrior=options.useNormInvWishPrior();
			_useHyperpriorR1=options.useHyperpriorR1();
			_useIndependentNormal =options.useIndependentNormal();
			_useSeparationPrior = options.useSeparationPrior();
			return *this;
		}


	private:
		// The input file name
		string _inFileName;
		// The output file stem
		string _outFileStem;
		// The hyper parameter file
		string _hyperParamFileName;
		// The predictions file
		string _predictFileName;
		// The number of sweeps
		unsigned int _nSweeps;
		// The length of burn in
		unsigned int _nBurn;
		// Whether to report the burn in period
		bool _reportBurnIn;
		// How often the output should be filtered for reporting
		unsigned int _nFilter;
		// The number of iterations to print how sampler is progressing
		unsigned int _nProgress;
		// The number of initial clusters to initialise into
		unsigned int _nClusInit;
		// The random number seed
		uint_fast32_t _seed;
		// The model for the outcome
		string _outcomeType;
		// The model for the covariates
		string _covariateType;
		// This notes whether we are including the response
		bool _includeResponse;
		// This notes which label switching moves are run
		string _whichLabelSwitch;
		// This has a fixed value of alpha (if negative we update alpha)
		double _fixedAlpha;
		// This has a fixed value of the discount parameter of the Pitman Yor process prior
		double _dPitmanYor;
		// The method used by the sampler
		string _samplerType;
		// This notes whether we are also doing predictions
		bool _doPrediction;
		// This notes whether we have extra variation in the response
		bool _responseExtraVar;
		// The type of extra variation in the response
		string _varSelectType;
		// This notes whether we are computing entropy
		bool _computeEntropy;
		// This notes whether we are including CAR random term
		bool _includeCAR;
		// This notes whether we are including initialisation values for uCAR
		bool _includeuCARinit;
		// Neighbours structure file name
		string _neighbourFileName;
		// uCAR initialisation file name
		string _uCARinitFileName;
		// For Poisson response and spatial CAR, to choose whether to use  the adaptive rejection sampler
		bool _PoissonCARadaptive;
		// The type of predictions (RaoBlackwell or random - which is only for yModel=Normal or yModel=Quantile)
		string _predictType;
		// For Survival response, whether the weibull shape parameter is fixed or cluster specific
		bool _weibullFixedShape;
		// Whether the conjugated Normal Inverse Wishart prior is used instead of the independent Normal and inverse Wisharts priors
		bool _useNormInvWishPrior;
		// Whether we include a hyperprior for R1
		bool _useHyperpriorR1;
		// Whether we use the independent conditional likelihood for continuous variables
		bool _useIndependentNormal;
		// Whether we use the separation prior for the Wishart distriution for xModel=Normal or xModel=Mixed
		bool _useSeparationPrior;
};

#endif // DIPBACOPTIONS_H_
