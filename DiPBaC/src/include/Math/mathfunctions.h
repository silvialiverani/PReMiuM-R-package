/// \file mathfunctions.h
/// \author David Hastie
/// \date 19 Mar 2012
/// \brief Header file to define distributions

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


#ifndef MATHFUNCTIONS_H_
#define MATHFUNCTIONS_H_

#include<cmath>

#include<boost/math/special_functions/gamma.hpp>

using namespace boost::math::constants;

using boost::math::lgamma;

double logMultivarGammaFn(const double& x,const unsigned int& p){

	double out;
	out = 0.25*(double)(p*(p-1))*log(pi<double>());
	for(unsigned int i=1;i<=p;i++){
		out += lgamma(x+(1.0-(double)i)/2.0);
	}
	return out;
}

double logit(const double& lambda){
	return 1.0/(1.0+exp(-lambda));
}

#endif /*MATHFUNCTIONS_H_*/
