/// \file sampler.h
/// \author David Hastie
/// \date 30 Sep 2010
/// \brief Header file defining classes for Markov chain Monte Carlo sampler.

/// \note (C) Copyright David Hastie and Silvia Liverani, 2012.

/// DiPBaC++ is free software; you can redistribute it and/or modify it under the
/// terms of the GNU Lesser General Public License as published by the Free Software
/// Foundation; either version 3 of the License, or (at your option) any later
/// version.

/// DiPBaC++ is distributed in the hope that it will be useful, but WITHOUT ANY
/// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
/// PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

/// You should have received a copy of the GNU Lesser General Public License
/// along with DiPBaC++ in the documentation directory. If not, see
/// <http://www.gnu.org/licenses/>.

/// The external linear algebra library Eigen, parts of which are included  in the
/// lib directory is released under the LGPL3+ licence. See comments in file headers
/// for details.

/// The Boost C++ header library, parts of which are included in the  lib directory
/// is released under the Boost Software Licence, Version 1.0, a copy  of which is
/// included in the documentation directory.


#ifndef SAMPLER_H_
#define SAMPLER_H_

// Standard includes
#include<vector>
#include<ctime>
#include<string>
#include<cassert>
#include<iostream>
#include<cstdlib>
#include<fstream>
#include<sstream>

// Custom includes
#include<Math/random.h>
#include<MCMC/model.h>
#include<MCMC/chain.h>
#include<MCMC/proposal.h>

using std::string;
using std::vector;
using std::endl;
using std::cout;
using std::ostringstream;

/// \class mcmcSampler sampler.h "MCMC/sampler.h"
/// \brief Class to implement an MCMC sampler
/// \param[in] modelParamType The user defined type to contain model parameters
/// \param[in] dataType The user defined type to contain a data observation
///	(can be omitted for problems with no data)
template<class modelParamType,class optionType,class propParamType,class dataType = int> class mcmcSampler{

	public:
		/// \brief Default constructor
		mcmcSampler() : _chain(), _model(), _proposalVec(),_rndGenerator(),_outFiles() {
			_nSweeps = 0;
			_nBurn = 0;
			_nFilter = 1;
			_reportBurnIn = false;
			_outFileStem = "output";
		}

		/// \brief Explicit constructor
		/// \param[in] nSweeps Number of sweeps
		/// \param[in] nBurn Number of sweeps for burn in
		/// \param[in] nFilter Filter for writing output (optional, default =1)
		/// \note nFilter filters the output for writing. If nFilter = 1, the
		/// output is written at the end of each sweep
		mcmcSampler(const unsigned int& nSweeps,const unsigned int& nBurn,
					const unsigned int& nFilter=1)
				:  _chain(), _model(), _proposalVec(), _rndGenerator(),_outFiles() {
			_nSweeps = nSweeps;
			_nBurn = nBurn;
			_nFilter = nFilter;
			_reportBurnIn = false;
			_outFileStem = "output";
		}

		/// \brief Destructor
		~mcmcSampler(){};

		/// \brief Member function to set the number of sweeps for the sampler
		/// \param[in] nS The number of sweeps
		void nSweeps(const unsigned int& nS){
			_nSweeps = nS;
		}

		/// \brief Member function to set the number of burn in sweeps for the sampler
		/// \param[in] nB The number of burn in sweeps
		void nBurn(const unsigned int& nB){
			_nBurn = nB;
		}

		unsigned int nBurn() const{
			return _nBurn;
		}

		/// \brief Member function to set the number of sweeps at which to report progress update
		/// \param[in] nPr The number of sweeps
		void nProgress(const unsigned int& nPr){
			_nProgress = nPr;
		}

		/// \brief Member function to set the number of sweeps that the
		/// output is filtered for writing
		/// \param[in] nF The filter to apply
		void nFilter(const unsigned int& nF){
			_nFilter = nF;
		}

		unsigned int nFilter() const{
			return _nFilter;
		}

		/// \brief Member function to define whether we want to report output
		/// from the burn in sweeps
		/// \param[in] repBurnIn true if we want output from the burn in, false
		/// otherwise
		void reportBurnIn(const bool& repBurnIn){
			_reportBurnIn = repBurnIn;
		}

		bool reportBurnIn() const{
			return _reportBurnIn;
		}


		/// \brief Member function to set the model options
		/// \param[in] modelOpts An object of optionsType
		void options(const optionType& modelOpts){
			_model.options(modelOpts);
		}

		/// \brief Member function to set the proposal parameters
		/// \param[in] modelOpts An object of optionsType
		void proposalParams(const propParamType& proposalParams){
			_proposalParams = proposalParams;
		}

		const propParamType& proposalParams() const{
			return _proposalParams;
		}

		propParamType& proposalParams(){
			return _proposalParams;
		}

		/// \brief Member function to set the output file path
		/// \param[in] fileStem The output file path
		void initialiseOutputFiles(const string& fileStem){
			_outFileStem = fileStem;
			string logFileName = fileStem + "_log.txt";
			_logFile.open(logFileName.c_str());
		}

		/// \brief Member function to close the output files
		void closeOutputFiles(){
			_logFile.close();
			for(unsigned int i=0;i<_outFiles.size();i++){
				(*(_outFiles[i])).close();
				delete _outFiles[i];
			}
		}

		const string& outFileStem() const{
			return _outFileStem;
		}

		vector<std::ofstream*>& outFiles(){
			return _outFiles;
		}

		/// \brief Member function to set the user output function
		/// \param[in] f Pointer to user function for writing the output
		void userOutputFn(void (*f)(mcmcSampler<modelParamType,optionType,propParamType,dataType>&,
										const unsigned int&)){
			_writeOutput = f;
		}

		// Full comments with function definition below
		void model(void (*initialiseParamsFn)(baseGeneratorType&,const mcmcModel<modelParamType,optionType,dataType>&,modelParamType&),
						vector<double> (*logPostFn)(const modelParamType&,const mcmcModel<modelParamType,optionType,dataType>&));

		// Full comments with function definition below
		void model(void (*importDataFn)(const string&,const string&,dataType&),
						void (*initialiseParamsFn)(baseGeneratorType&,const mcmcModel<modelParamType,optionType,dataType>&,modelParamType&),
						vector<double> (*logPostFn)(const modelParamType&,const mcmcModel<modelParamType,optionType,dataType>&),
						const bool&);

		/// \brief Return the mcmc model
		const mcmcModel<modelParamType,optionType,dataType>& model() const{
			return _model;
		}

		/// \brief Return the mcmc model
		mcmcModel<modelParamType,optionType,dataType>& model(){
			return _model;
		}

		const mcmcChain<modelParamType>& chain() const{
			return _chain;
		}

		/// \brief Set the missing data function for the model
		void updateMissingDataFn(void (*f)(baseGeneratorType&,
										modelParamType&,
										const optionType&,
										dataType&)){
			_model.updateMissingDataFunction(f);
		}

		// Full comments with function definition below
		void addProposal(const string& propName,const double& weight,
						const unsigned int& frequency,
						const unsigned int& firstSweep,
						void (*updateFn)(mcmcChain<modelParamType>&,
										unsigned int&,unsigned int&,
										const mcmcModel<modelParamType,optionType,dataType>&,
										propParamType&, baseGeneratorType&));

		/// \brief Member function to seed the random number generator
		/// \param[in] seedValue The seed value to use (if 0 is used the clock is
		/// used to seed the generator)
		void seedGenerator(const unsigned int& seedValue = 0){
			if(seedValue==0){
				_seed = static_cast<unsigned int>(time(NULL));
			}else{
				_seed = seedValue;
			}
			cout << "Random number seed: " << _seed << endl;
			_rndGenerator.seed(_seed);
		}

		/// \brief Member function to initialise the MCMC chain
		void initialiseChain(){
			modelParamType tmpModelParams;
			_model.initialiseParams(_rndGenerator,tmpModelParams);
			vector<double> logPostVec;
			logPostVec = _model.logPosterior(tmpModelParams);
			mcmcState<modelParamType> tmpState(tmpModelParams,logPostVec);
			_chain.currentState(tmpState);
		}

		/// \brief Member function to initialise the proposal parameters
		void initialiseProposalParams(const propParamType& propParams){
			_proposalParams=propParams;
		}

		/// \brief Member function to import the data
		/// \param[in] dataFilePath The file path where the data is stored
 		void importData(const string& dataFilePath,const string& predictFilePath){
			_dataFilePath=dataFilePath;
			_predictFilePath=predictFilePath;
			_model.importData(dataFilePath,predictFilePath);
		}

		/// \brief Member function to updateMissingData the data
		void updateMissingData(){
			if(_model.hasMissingData()){
				_model.updateMissingData(_rndGenerator,_chain.currentState().parameters());
			}
		}


		// Full comments with function definition below
		void writeLogFile();

		// Full comments with function definition below
		void run();

		/// \brief Member function to get the history of the mcmcChain associated
		/// with the sampler
		vector<mcmcState<modelParamType> > getSamplerHistory(){
			return _chain.history();
		}

		/// \brief Member function to append to log file
		void appendToLogFile(const string& logString){
			_logFile << logString;
		}

	private:
		/// \brief The mcmcChain associated with the sampler
		mcmcChain<modelParamType> _chain;

		/// \brief The model that the sampler is for
		mcmcModel<modelParamType,optionType,dataType> _model;

		/// \brief The proposal parameters used in the sampler
		propParamType _proposalParams;

		/// \brief The vector of proposals that are used for the sampler
		vector<mcmcProposal<modelParamType,optionType,propParamType,dataType> > _proposalVec;

		/// \brief A boost random number generator
		baseGeneratorType _rndGenerator;

		/// \brief The seed for the random number generator
		unsigned int _seed;

		/// \var _nSweeps
		/// \brief The number of sweeps that we are running the sampler for

		/// \var _nBurn
		/// \brief The number of sweeps that we are discarding as burn in

		/// \var _nFilter
		/// \brief The number of sweeps that we are filtering when writing output

		/// \var _nProgress
		/// \brief The number of sweeps at which to print a progress update
		unsigned int _nSweeps,_nBurn,_nFilter,_nProgress;

		/// \brief The input file path
		string _dataFilePath;

		/// \brief The output file stem
		string _outFileStem;

		/// \brief The input file path
		string _predictFilePath;

		/// \var _reportBurnIn
		/// \brief Boolean to indicate whether we want to report back output
		/// during the burn in period
		bool _reportBurnIn;

		/// \brief Pointer to user function to process the output
		void (*_writeOutput)(mcmcSampler<modelParamType,optionType,propParamType,dataType>&,
								const unsigned int&);

		/// \var _logPostFile
		/// \brief File object for writing the log posterior, log likelihood and
		/// log prior values

		/// \var _acceptFile
		/// \brief File object for writing the acceptance probabilities of proposals

		/// \var _logFile
		/// \brief File object for writing the log of the run
		std::ofstream _logFile;

		vector<std::ofstream*> _outFiles;

		/// \brief Private member function for writing the output of the sampler
		/// \note This is a just a wrapper for the user supplied function
		void writeOutput(const unsigned int&);

		/// \brief Private member function for writing the acceptance rates of the sampler
		void writeAcceptanceRates();


};

/// \brief Member function to set the model for the sampler for problems
/// with no data
/// \param[in] initialiseParamsFn Pointer to function for initialising
/// model parameters
/// \param[in] targetFn Pointer to function for evaluating the log posterior
template<class modelParamType,class optionType,class propParamType,class dataType>
void mcmcSampler<modelParamType,optionType,propParamType,dataType>::
	model(void (*initialiseParamsFn)(baseGeneratorType&,const mcmcModel<modelParamType,optionType,dataType>&,modelParamType&),
				vector<double> (*logPostFn)(const modelParamType&,const mcmcModel<modelParamType,optionType,dataType>&)){

	// Setting up the model associated with the sampler for problems when there
	// is no data

	// First call the mcmcModel method to assign the function for initialising the model
	// parameters
	_model.initialiseParamsFunction(initialiseParamsFn);

	// Now call the mcmcModel method to assign the function to calculate the
	// log posterior of a modelParamType object
	_model.logPosteriorFunction(logPostFn);

	// Now set the indicator whether the model has missing data
	_model.hasMissingData(false);


}

/// \brief Member function to set the model for the sampler for problems
/// with data
/// \param[in] importDataFn Pointer to function for importing the data
/// \param[in] initialiseParamsFn Pointer to function for initialising
/// model parameters
/// \param[in] logPostFn Pointer to function for evaluating the log posterior
template<class modelParamType,class optionType,class propParamType,class dataType>
void mcmcSampler<modelParamType,optionType,propParamType,dataType>::
	model(void (*importDataFn)(const string&,const string&,dataType&),
				void (*initialiseParamsFn)(baseGeneratorType&,const mcmcModel<modelParamType,optionType,dataType>&,modelParamType&),
				vector<double> (*logPostFn)(const modelParamType&,const mcmcModel<modelParamType,optionType,dataType>&),
				const bool& hasMissingData){

	// Setting up the model associated with the sampler for problems when there
	// is data

	// First call the mcmcModel method to assign the function for importing the data
	_model.importDataFunction(importDataFn);

	// Now call the mcmcModel method to assign the function for initialising the model
	// parameters
	_model.initialiseParamsFunction(initialiseParamsFn);

	// Now call the mcmcModel method to assign the function to calculate the
	// log posterior of a modelParamType object
	_model.logPosteriorFunction(logPostFn);

	// Now set the indicator whether the model has missing data
	_model.hasMissingData(hasMissingData);

}

/// \brief Member function to add a proposal to the sampler
/// \param[in] proposalParameters Initial value of proposal parameters
/// \param[in] weight Probability proposal is used if it is attempted
/// \param[in] frequency Frequency in sweeps that this proposal is attempted
/// \param[in] updateFrequency Frequency in proposal attempts that the proposal
/// parameters are updated (for adaptive proposals)
/// \param[in] updateFn Pointer to the function for updating the chain
template<class modelParamType,class optionType,class propParamType,class dataType>
void mcmcSampler<modelParamType,optionType,propParamType,dataType>::
	addProposal(const string& propName,const double& weight,
				const unsigned int& frequency,
				const unsigned int& firstSweep,
				void (*updateFn)(mcmcChain<modelParamType>&,
								unsigned int&, unsigned int&,
								const mcmcModel<modelParamType,optionType,dataType>&,
								propParamType&, baseGeneratorType&)){
	// Adding a proposal to the sampler

	// First define a temporary object of type mcmcProposal
	mcmcProposal<modelParamType,optionType,propParamType,dataType> tmpProposal;

	// Now call the mcmcProposal method to initialise the proposal
	tmpProposal.initialiseProposal(propName,&_proposalParams,weight,frequency,firstSweep);

	// Now call the mcmcProposal method to set the function for updating model
	// parameters
	tmpProposal.updateFunction(updateFn);

	// Finally we add this proposal to the vector of proposals for this sampler
	_proposalVec.push_back(tmpProposal);

}

/// \brief Member function to write a log file for the sampler
template<class modelParamType,class optionType,class propParamType,class dataType>
void mcmcSampler<modelParamType,optionType,propParamType,dataType>::writeLogFile(){

	// Function for writing the log file for reproducing the run
	// First of all get the current time
	time_t now;
	struct tm *timeinfo;
	char theDate[20];
	time(&now);
	timeinfo = localtime(&now);
	strftime(theDate,20,"%Y-%m-%d %H:%M%S",timeinfo);

	// Write the date and time
	_logFile << "Date and time: " << theDate << "\n\n";
	// Write the file path
	_logFile << "Data file path: " << _dataFilePath << "\n\n";
	_logFile << "Output file path: " << _outFileStem << "\n\n";
	if(_predictFilePath.compare("")!=0){
		_logFile << "Prediction file path: " << _predictFilePath << "\n\n";
	}else{
		_logFile << "Predict file path: No predictions run.\n\n";
	}
	// Write the seed that was used
	_logFile << "Seed: " << _seed << "\n\n";
	// Write the user inputed options
	_logFile << "Number of sweeps: " << _nSweeps << "\n";
	_logFile << "Burn in sweeps: " << _nBurn << "\n";
	_logFile << "Output filter: " << _nFilter << "\n";

}

template<class modelParamType,class optionType,class propParamType,class dataType>
void mcmcSampler<modelParamType,optionType,propParamType,dataType>::writeOutput(const unsigned int& sweep){
	// Function for writing the output of the chain

	// Write the output
	(*_writeOutput)(*this,sweep);

}

template<class modelParamType,class optionType,class propParamType,class dataType>
void mcmcSampler<modelParamType,optionType,propParamType,dataType>::writeAcceptanceRates(){

	// Write the acceptance rate of each of the proposals into the
	// log file. The log file is opened externally in a member function
	ostringstream tmpStr;
	typename vector<mcmcProposal<modelParamType,optionType,propParamType,dataType> >::iterator it;
	for(it=_proposalVec.begin();it<_proposalVec.end();++it){
		tmpStr << "Proposal Type: " << it->proposalName() << ", Acceptance Rate: " << it->acceptanceRate()<< endl;
	}
	appendToLogFile(tmpStr.str());


}

/// \brief Member function to run the sampler
/// \note This is the main function that does all the work. The function
/// itself is commented with normal C++ style comments to explain what
/// it does
template<class modelParamType,class optionType,class propParamType,class dataType>
void mcmcSampler<modelParamType,optionType,propParamType,dataType>::run(){
	// This is the main function that runs the sampler

	// Define a uniform random number generator
	randomUniform unifRand(0,1);

	// This is the main sampler
	// We loop over the sweeps

	// Write the output of initialisation before sampler begins
	writeOutput(0);
	for(unsigned int sweep=1; sweep<=_nBurn+_nSweeps; sweep++){
		if(sweep==1||sweep%_nProgress==0){
			cout << "Sweep: " << sweep << endl;
		}
		// Update the missing data (this will only do anything if the
		// _model.hasMissingData flag is true)
		updateMissingData();

		// At each sweep we loop over the proposals
		typename vector<mcmcProposal<modelParamType,optionType,propParamType,dataType> >::iterator it;
		for(it=_proposalVec.begin(); it<_proposalVec.end(); it++){
			// Only use this proposal if it is due to be tried at this sweep
			if(sweep >= it->proposalFirstSweep() && sweep % it->proposalFrequency()==0){

				// Only try this proposal with probability as defined
				if(unifRand(_rndGenerator)<it->proposalWeight()){

					// Update the chain state
					it->updateParameters(_chain,_model,_rndGenerator);
				}
			}

		}

		// At the end of the sweep make sure the log posterior is up to date.
		_chain.currentState().logPosterior(_model.logPosterior(_chain.currentState().parameters()));

		// Now write the output (this is controlled by the user defined function
		writeOutput(sweep);

	}
	writeAcceptanceRates();

}

#endif /* SAMPLER_H_ */
