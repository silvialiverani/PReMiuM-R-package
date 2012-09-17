/// \file model.h
/// \author David Hastie
/// \date 30 Sep s2010
/// \brief Header file defining class for Markov chain Monte Carlo model.

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


#ifndef MODEL_H_
#define MODEL_H_

// Standard includes
#include<vector>

// Custom includes
#include<Math/random.h>
#include<MCMC/state.h>

using std::vector;
using std::string;
using std::cout;
using std::endl;

/// \class mcmcModel model.h "MCMC/model.h"
/// \brief Template class for Markov chain Monte Carlo models
/// \param[in] modelParamType The user defined type to contain model parameters
/// \param[in] dataType The user defined type to contain a data observation
///	(can be omitted for problems with no data)
template<class modelParamType,class optionType,class dataType = vector<int> > class mcmcModel{

	public:
		/// \brief Default constructor
		mcmcModel() : _dataset() {};

		/// \brief Destructor
		~mcmcModel(){};

		/// \brief Member function for accessing the model options
		const optionType& options() const{
			return _options;
		}

		/// \brief member function for accessing the dataset
		const dataType& dataset() const{
			return _dataset;
		}

		/// \brief member function for accessing the dataset
		dataType& dataset(){
			return _dataset;
		}


		/// \brief Member function for setting the model options
		void options(const optionType& modelOpts){
			_options = modelOpts;
		}

		/// \brief Member function for returning whether there is missing data
		bool hasMissingData(){
			return _hasMissingData;
		}

		/// \brief Member function for setting whether there is missing data
		void hasMissingData(const bool& missingDataInd){
			_hasMissingData = missingDataInd;
		}

		/// \brief Member function for assigning the (user defined) function for
		///	importing data
		/// \param[in] *f Pointer to a function for importing data
		void importDataFunction(void (*f)(const string&,const string&,dataType&)){
					_importData=f;
		}

		/// \brief Member function for assigning the (user defined) function for
		///	updating missing data
		/// \param[in] *f Pointer to a function for updating missing data
		void updateMissingDataFunction(void (*f)(baseGeneratorType&,
				modelParamType&,const optionType& options,dataType& dataset)){
			_updateMissingData=f;
		}

		/// \brief Member function for assigning the (user defined) function for
		///	initialising parameters
		/// \param[in] *f Pointer to a function for initialising the parameters
		void initialiseParamsFunction(void (*f)(baseGeneratorType&,
											const mcmcModel<modelParamType,optionType,dataType>&
											,modelParamType&)){
			_initialiseParams=f;
		}

		/// \brief Member function for assigning the (user defined) log posterior when
		///	there is data
		/// \param[in] *f Pointer to a function for calculating the log
		///	Posterior given the current model parameters and data
		void logPosteriorFunction(vector<double> (*f)(const modelParamType&,
													const mcmcModel<modelParamType,optionType,dataType>&)){
			_logPosterior=f;
		}

		/// \brief Member function as wrapper to import the data
		/// \param[in] fitFilePath Name of the filepath where fitting data should be imported from
		/// \param[in] predictFilePath Name of the filepath where predictive data should be imported from
		void importData(const string& fitFilePath,const string& predictFilePath){
			(*_importData)(fitFilePath,predictFilePath,_dataset);
		}

		/// \brief Member function as wrapper to update missing data
		/// \param[in] rndGenerator A boost random number generator
		/// \param[out] modelParams The initialised parameters
		void updateMissingData(baseGeneratorType& rndGenerator,modelParamType& modelParams){
			if(_hasMissingData){
				(*_updateMissingData)(rndGenerator,modelParams,_options,_dataset);
			}
		}

		/// \brief Member function as wrapper to initialise the parameters
		/// \param[in] rndGenerator A boost random number generator
		/// \param[out] modelParams The initialised parameters
		void initialiseParams(baseGeneratorType& rndGenerator,modelParamType& modelParams){
			(*_initialiseParams)(rndGenerator,*this,modelParams);

		}

		/// \brief Member function wrapper to return log posterior. The function also
		/// returns log-likelihood and log-prior in 2nd and 3rd elements
		/// \param[in] parameters Model parameter values at which log posterior is to be
		/// evaluated
		/// \return A vector containing the log posterior in the first element,
		/// the log likelihood in the second element and the log prior in the
		/// third element. For problems with no data the last two elements are 0
		vector<double> logPosterior(const modelParamType& parameters) const{
			return (*_logPosterior)(parameters,*this);
		}

	private:
		/// \brief The user provided dataset of (user defined) dataType objects
		dataType _dataset;

		/// \brief An indicator of whether there is missing data in the dataset
		bool _hasMissingData;

		/// \brief The options for the model
		optionType _options;

		/// \brief Pointer to function for importing the data given a filePath
		void (*_importData)(const string&,const string&,dataType&);

		/// \brief Pointer to function for updating missing data given random number
		/// generator, model parameters and model options
		void (*_updateMissingData)(baseGeneratorType& rndGenerator,modelParamType&,const optionType&,dataType& dataset);

		/// \brief Pointer to function for initialising the parameters given a
		/// random number generator and a dataset
		void (*_initialiseParams)(baseGeneratorType& rndGenerator,const mcmcModel<modelParamType,optionType,dataType>&,modelParamType&);


		/// \brief Pointer to function for calculating the log posterior when
		/// there is data
		vector<double> (*_logPosterior)(const modelParamType&, const mcmcModel<modelParamType,optionType,dataType>&);


};

#endif /* MODEL_H_ */
