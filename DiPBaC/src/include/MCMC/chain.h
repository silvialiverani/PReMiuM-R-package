/// \file chain.h
/// \author David Hastie
/// \date 30 Sep 2010
/// \brief Header file defining classes for Markov chain Monte Carlo chains.

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


#ifndef CHAIN_H_
#define CHAIN_H_

// Standard includes
#include<fstream>
#include<string>
#include<vector>

// Custom includes
#include<MCMC/state.h>

using std::vector;
using std::cout;
using std::ostream;

/// \class mcmcChain chain.h "MCMC/chain.h"
/// \brief Class for Markov chain Monte Carlo chains
template<class modelParamType> class mcmcChain{

	public:
		/// \brief Default constructor
		mcmcChain() {};

		/// \brief Explicit constructor
		/// \param[in] currentState The initial state of the chain
		mcmcChain(const mcmcState<modelParamType>& currentState) :
			_currentState(currentState) {};

		/// \brief Destructor
		~mcmcChain(){};

		/// \brief Member function for setting the current state of the chain
		/// \param[in] x The state to be made the current state
		void currentState(const mcmcState<modelParamType>& x){
			_currentState = x;
		}

		/// \brief Member function to get the current state of the chain
		/// \return The current state of the chain
		const mcmcState<modelParamType>& currentState() const{
			return _currentState;
		}

		/// \brief Member function to get the current state of the chain
		/// \return The current state of the chain
		mcmcState<modelParamType>& currentState(){
			return _currentState;
		}


	private:
		/// \brief Current value of the chain
		mcmcState<modelParamType> _currentState;

};


#endif /* CHAIN_H_ */
