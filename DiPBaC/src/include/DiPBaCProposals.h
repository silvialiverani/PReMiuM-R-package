/// \file DiPBaCProposals.h
/// \author David Hastie
/// \brief Header file for model specification for DiPBaCpp

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

#ifndef DIPBACPROPOSALS_H_
#define DIPBACPROPOSALS_H_

// Standard includes
#include<cmath>
#include<vector>
#include<iostream>
#include<string>
#include<numeric>
#include<limits>

#include<boost/math/distributions/normal.hpp>
#include<boost/math/distributions/beta.hpp>
#include<boost/math/distributions/students_t.hpp>
#include<boost/math/special_functions/gamma.hpp>

#include<Eigen/Core>
#include<Eigen/Cholesky>
#include<Eigen/LU>

// Custom includes
#include<MCMC/chain.h>
#include<MCMC/model.h>
#include<Math/random.h>
#include<Math/distribution.h>
#include<DiPBaCOptions.h>
#include<DiPBaCModel.h>
#include<DiPBaCData.h>

using namespace Eigen;

using std::vector;
using std::accumulate;
using std::numeric_limits;
using boost::math::normal_distribution;
using boost::math::students_t_distribution;
using boost::math::lgamma;


class diPBaCPropParams{

	public:
		// Default constructor
		diPBaCPropParams() {};

		diPBaCPropParams(const unsigned int& nSweeps,const unsigned int& nCovariates,
								const unsigned int& nFixedEffects,const unsigned int& nCategoriesY){
				_thetaStdDev=1.0;
				_thetaStdDevLower=0.1;
				_thetaStdDevUpper=99.9;
				_nTryTheta=0;
				_nAcceptTheta=0;
				_nLocalAcceptTheta=0;
				_nResetTheta=0;
				_thetaAcceptTarget=0.44;
				_thetaUpdateFreq=25;
				_thetaAnyUpdates=true;

				_nTryBeta.resize(nFixedEffects);
				_nAcceptBeta.resize(nFixedEffects);
				_nLocalAcceptBeta.resize(nFixedEffects);
				_nResetBeta.resize(nFixedEffects);
				_betaStdDev.resize(nFixedEffects);
				_betaStdDevLower.resize(nFixedEffects);
				_betaStdDevUpper.resize(nFixedEffects);
				for(unsigned int j=0;j<nFixedEffects;j++){
					_betaStdDev[j]=1.0;
					_betaStdDevLower[j]=0.1;
					_betaStdDevUpper[j]=99.9;
					_nTryBeta[j]=0;
					_nAcceptBeta[j]=0;
					_nLocalAcceptBeta[j]=0;
					_nResetBeta[j]=0;
				}
				_betaAcceptTarget = 0.44;
				_betaUpdateFreq = 25;
				_betaAnyUpdates=true;

				_alphaStdDev=1.0;
				_alphaStdDevLower=0.1;
				_alphaStdDevUpper=99.9;
				_nTryAlpha=0;
				_nAcceptAlpha=0;
				_nLocalAcceptAlpha=0;
				_nResetAlpha=0;
				_alphaAcceptTarget = 0.44;
				_alphaUpdateFreq = 25;
				_alphaAnyUpdates=true;


				_nTryRho.resize(nCovariates);
				_nAcceptRho.resize(nCovariates);
				_nLocalAcceptRho.resize(nCovariates);
				_nResetRho.resize(nCovariates);
				_rhoStdDev.resize(nCovariates);
				_rhoStdDevLower.resize(nCovariates);
				_rhoStdDevUpper.resize(nCovariates);
				for(unsigned int j=0;j<nCovariates;j++){
					_rhoStdDev[j]=0.5;
					_rhoStdDevLower[j]=0.0001;
					_rhoStdDevUpper[j]=9.9999;
					_nTryRho[j]=0;
					_nAcceptRho[j]=0;
					_nLocalAcceptRho[j]=0;
					_nResetRho[j]=0;
				}
				_rhoAcceptTarget = 0.44;
				_rhoUpdateFreq = 10;
				_rhoAnyUpdates=true;


				_lambdaStdDev=1.0;
				_lambdaStdDevLower=0.1;
				_lambdaStdDevUpper=99.9;
				_nTryLambda=0;
				_nAcceptLambda=0;
				_nLocalAcceptLambda=0;
				_nResetLambda=0;
				_lambdaAcceptTarget = 0.44;
				_lambdaUpdateFreq = 500;
				_lambdaAnyUpdates=true;

		};

		~diPBaCPropParams(){};


		unsigned int nTryTheta() const{
			return _nTryTheta;
		}

		unsigned int nAcceptTheta() const{
			return _nAcceptTheta;
		}

		double thetaAcceptRate() const{
			if(_nTryTheta>0){
				return (double)_nAcceptTheta/(double)_nTryTheta;
			}else{
				return 0.0;
			}
		}

		unsigned int thetaUpdateFreq() const{
			return _thetaUpdateFreq;
		}

		unsigned int nLocalAcceptTheta() const{
			return _nLocalAcceptTheta;
		}

		double thetaLocalAcceptRate() const{
					return (double)_nLocalAcceptTheta/(double)_thetaUpdateFreq;
		}


		double thetaAcceptTarget() const{
			return _thetaAcceptTarget;
		}

		void thetaAddTry(){
			_nTryTheta++;
		}

		void thetaAddAccept(){
			_nAcceptTheta++;
			_nLocalAcceptTheta++;
		}

		void thetaLocalReset(){
			_nLocalAcceptTheta=0;
		}

		unsigned int nResetTheta() const{
			return _nResetTheta;
		}

		void thetaStdDevReset(){
			_thetaStdDev = 1.0;
			_nResetTheta++;
			_thetaStdDevLower = pow(10.0,-((double)_nResetTheta+1.0));
			_thetaStdDevUpper = 100.0-pow(10.0,-((double)_nResetTheta+1.0));
		}

		double& thetaStdDev(){
			return _thetaStdDev;
		}

		double thetaStdDev() const{
			return _thetaStdDev;
		}

		// Member function for setting the standard deviation for
		// proposal for theta
		void thetaStdDev(const double& sd){
			_thetaStdDev=sd;
		}

		double thetaStdDevLower() const{
			return _thetaStdDevLower;
		}

		double thetaStdDevUpper() const{
			return _thetaStdDevUpper;
		}

		bool thetaAnyUpdates() const{
			return _thetaAnyUpdates;
		}

		void thetaAnyUpdates(const bool& newStatus){
			_thetaAnyUpdates = newStatus;
		}


		vector<unsigned int> nTryBeta() const{
			return _nTryBeta;
		}

		unsigned int nTryBeta(const unsigned int& j) const{
			return _nTryBeta[j];
		}

		vector<unsigned int> nAcceptBeta() const{
			return _nAcceptBeta;
		}

		double betaAcceptRate(const unsigned int& j) const{
			if(_nTryBeta[j]>0){
				return (double)_nAcceptBeta[j]/(double)_nTryBeta[j];
			}else{
				return 0.0;
			}
		}

		unsigned int betaUpdateFreq() const{
			return _betaUpdateFreq;
		}

		vector<unsigned int> nLocalAcceptBeta() const{
			return _nLocalAcceptBeta;
		}

		double betaLocalAcceptRate(const unsigned int& j) const{
				return (double)_nLocalAcceptBeta[j]/(double)_betaUpdateFreq;
		}

		double betaAcceptTarget() const{
			return _betaAcceptTarget;
		}

		void betaAddTry(const unsigned int& j){
			_nTryBeta[j]++;
		}

		void betaAddAccept(const unsigned int& j){
			_nAcceptBeta[j]++;
			_nLocalAcceptBeta[j]++;
		}

		void betaLocalReset(const unsigned int& j){
			_nLocalAcceptBeta[j]=0;
		}

		vector<unsigned int> nResetBeta() const{
			return _nResetBeta;
		}

		void betaStdDevReset(const unsigned int& j){
			_betaStdDev[j] = 1.0;
			_nResetBeta[j]++;
			_betaStdDevLower[j] = pow(10.0,-((double)_nResetBeta[j]+1.0));
			_betaStdDevUpper[j] = 100.0-pow(10.0,-((double)_nResetBeta[j]+1.0));
		}

		vector<double> betaStdDev() const{
			return _betaStdDev;
		}

		vector<double>& betaStdDev(){
			return _betaStdDev;
		}

		double& betaStdDev(const unsigned int& j){
			return _betaStdDev[j];
		}

		const double& betaStdDev(const unsigned int& j) const{
			return _betaStdDev[j];
		}

		vector<double> betaStdDevLower() const{
			return _betaStdDevLower;
		}

		double betaStdDevLower(const unsigned int& j) const{
			return _betaStdDevLower[j];
		}

		vector<double> betaStdDevUpper() const{
			return _betaStdDevUpper;
		}

		double betaStdDevUpper(const unsigned int& j) const{
			return _betaStdDevUpper[j];
		}

		// Member function for setting the standard deviation for
		// proposal for beta for fixed effect j
		void betaStdDev(const unsigned int& j,const double& sd){
			_betaStdDev[j]=sd;
		}

		bool betaAnyUpdates() const{
			return _betaAnyUpdates;
		}

		void betaAnyUpdates(const bool& newStatus){
			_betaAnyUpdates = newStatus;
		}



		unsigned int nTryAlpha() const{
			return _nTryAlpha;
		}

		unsigned int nAcceptAlpha() const{
			return _nAcceptAlpha;
		}

		double alphaAcceptRate() const{
			if(_nTryAlpha>0){
				return (double)_nAcceptAlpha/(double)_nTryAlpha;
			}else{
				return 0.0;
			}
		}

		unsigned int alphaUpdateFreq() const{
			return _alphaUpdateFreq;
		}

		unsigned int nLocalAcceptAlpha() const{
			return _nLocalAcceptAlpha;
		}

		double alphaLocalAcceptRate() const{
			return (double)_nLocalAcceptAlpha/(double)_alphaUpdateFreq;
		}

		double alphaAcceptTarget() const{
			return _alphaAcceptTarget;
		}

		void alphaAddTry(){
			_nTryAlpha++;
		}

		void alphaAddAccept(){
			_nAcceptAlpha++;
			_nLocalAcceptAlpha++;
		}

		void alphaLocalReset(){
			_nLocalAcceptAlpha=0;
		}

		unsigned int nResetAlpha() const{
			return _nResetAlpha;
		}

		void alphaStdDevReset(){
			_alphaStdDev = 1.0;
			_nResetAlpha++;
			_alphaStdDevLower = pow(10.0,-((double)_nResetAlpha+1.0));
			_alphaStdDevUpper = 100.0-pow(10.0,-((double)_nResetAlpha+1.0));
		}


		const double alphaStdDev() const{
			return _alphaStdDev;
		}

		double& alphaStdDev(){
			return _alphaStdDev;
		}

		double alphaStdDevLower() const{
			return _alphaStdDevLower;
		}

		double alphaStdDevUpper() const{
			return _alphaStdDevUpper;
		}

		// Member function for setting the standard deviation for
		// proposal for beta for fixed effect j
		void alphaStdDev(const double& sd){
			_alphaStdDev=sd;
		}

		bool alphaAnyUpdates() const{
			return _alphaAnyUpdates;
		}

		void alphaAnyUpdates(const bool& newStatus){
			_alphaAnyUpdates = newStatus;
		}

		vector<unsigned int> nTryRho() const{
			return _nTryRho;
		}

		unsigned int nTryRho(const unsigned int& j) const{
			return _nTryRho[j];
		}

		vector<unsigned int> nAcceptRho() const{
			return _nAcceptRho;
		}

		double rhoAcceptRate(const unsigned int& j) const{
			if(_nTryRho[j]>0){
				return (double)_nAcceptRho[j]/(double)_nTryRho[j];
			}else{
				return 0.0;
			}
		}

		unsigned int rhoUpdateFreq() const{
			return _rhoUpdateFreq;
		}

		vector<unsigned int> nLocalAcceptRho() const{
			return _nLocalAcceptRho;
		}


		double rhoLocalAcceptRate(const unsigned int& j) const{
			return (double)_nLocalAcceptRho[j]/(double)_rhoUpdateFreq;
		}


		double rhoAcceptTarget() const{
			return _rhoAcceptTarget;
		}

		void rhoAddTry(const unsigned int& j){
			_nTryRho[j]++;
		}

		void rhoAddAccept(const unsigned int& j){
			_nAcceptRho[j]++;
			_nLocalAcceptRho[j]++;
		}

		void rhoLocalReset(const unsigned int& j){
			_nLocalAcceptRho[j]=0;
		}

		vector<unsigned int> nResetRho() const{
			return _nResetRho;
		}

		void rhoStdDevReset(const unsigned int& j){
			_rhoStdDev[j] = 0.5;
			_nResetRho[j]++;
			_rhoStdDevLower[j] = pow(10.0,-((double)_nResetRho[j]+4.0));
			_rhoStdDevUpper[j] = 10.0-pow(10.0,-((double)_nResetRho[j]+4.0));
		}

		vector<double>& rhoStdDev(){
			return _rhoStdDev;
		}

		vector<double> rhoStdDev() const{
			return _rhoStdDev;
		}

		double& rhoStdDev(const unsigned int& j){
			return _rhoStdDev[j];
		}

		const double& rhoStdDev(const unsigned int& j) const{
			return _rhoStdDev[j];
		}

		vector<double> rhoStdDevLower() const{
			return _rhoStdDevLower;
		}

		double rhoStdDevLower(const unsigned int& j) const{
			return _rhoStdDevLower[j];
		}

		vector<double> rhoStdDevUpper() const{
			return _rhoStdDevUpper;
		}

		double rhoStdDevUpper(const unsigned int& j) const{
			return _rhoStdDevUpper[j];
		}

		void rhoStdDev(const unsigned int& j,const double& sd){
			_rhoStdDev[j]=sd;
		}

		bool rhoAnyUpdates() const{
			return _rhoAnyUpdates;
		}

		void rhoAnyUpdates(const bool& newStatus){
			_rhoAnyUpdates = newStatus;
		}


				unsigned int nTryLambda() const{
			return _nTryLambda;
		}

		unsigned int nAcceptLambda() const{
			return _nAcceptLambda;
		}

		double lambdaAcceptRate() const{
			if(_nTryLambda>0){
				return (double)_nAcceptLambda/(double)_nTryLambda;
			}else{
				return 0.0;
			}
		}

		unsigned int nLocalAcceptLambda() const{
			return _nLocalAcceptLambda;
		}

		unsigned int lambdaUpdateFreq() const{
			return _lambdaUpdateFreq;
		}

		double lambdaLocalAcceptRate() const{
			return (double)_nLocalAcceptLambda/(double)_lambdaUpdateFreq;
		}

		double lambdaAcceptTarget() const{
			return _lambdaAcceptTarget;
		}

		void lambdaAddTry(){
			_nTryLambda++;
		}

		void lambdaAddAccept(){
			_nAcceptLambda++;
			_nLocalAcceptLambda++;
		}

		void lambdaLocalReset(){
			_nLocalAcceptLambda=0;
		}

		unsigned int nResetLambda() const{
			return _nResetLambda;
		}

		void lambdaStdDevReset(){
			_lambdaStdDev = 1.0;
			_nResetLambda++;
			_lambdaStdDevLower = pow(10.0,-((double)_nResetLambda+1.0));
			_lambdaStdDevUpper = 100.0-pow(10.0,-((double)_nResetLambda+1.0));
		}

		double lambdaStdDev() const{
			return _lambdaStdDev;
		}

		double& lambdaStdDev(){
			return _lambdaStdDev;
		}

		double lambdaStdDevLower() const{
			return _lambdaStdDevLower;
		}

		double lambdaStdDevUpper() const{
			return _lambdaStdDevUpper;
		}

		// Member function for setting the standard deviation for
		// proposal for beta for fixed effect j
		void lambdaStdDev(const double& sd){
			_lambdaStdDev=sd;
		}

		bool lambdaAnyUpdates() const{
			return _lambdaAnyUpdates;
		}

		void lambdaAnyUpdates(const bool& newStatus){
			_lambdaAnyUpdates = newStatus;
		}


		// Need to define a copy iterator
		diPBaCPropParams& operator=(const diPBaCPropParams& propParams){
			_nTryTheta=propParams.nTryTheta();
			_nAcceptTheta=propParams.nAcceptTheta();
			_nLocalAcceptTheta=propParams.nLocalAcceptTheta();
			_nResetTheta=propParams.nResetTheta();
			_thetaStdDev=propParams.thetaStdDev();
			_thetaStdDevLower=propParams.thetaStdDevLower();
			_thetaStdDevUpper=propParams.thetaStdDevUpper();
			_thetaAcceptTarget=propParams.thetaAcceptTarget();
			_thetaUpdateFreq=propParams.thetaUpdateFreq();
			_thetaAnyUpdates=propParams.thetaAnyUpdates();
			_nTryBeta=propParams.nTryBeta();
			_nAcceptBeta=propParams.nAcceptBeta();
			_nLocalAcceptBeta=propParams.nLocalAcceptBeta();
			_nResetBeta=propParams.nResetBeta();
			_betaStdDev=propParams.betaStdDev();
			_betaStdDevLower=propParams.betaStdDevLower();
			_betaStdDevUpper=propParams.betaStdDevUpper();
			_betaAcceptTarget=propParams.betaAcceptTarget();
			_betaUpdateFreq=propParams.betaUpdateFreq();
			_betaAnyUpdates=propParams.betaAnyUpdates();
			_nTryAlpha=propParams.nTryAlpha();
			_nAcceptAlpha=propParams.nAcceptAlpha();
			_nLocalAcceptAlpha=propParams.nLocalAcceptAlpha();
			_nResetAlpha=propParams.nResetAlpha();
			_alphaStdDev=propParams.alphaStdDev();
			_alphaStdDevLower=propParams.alphaStdDevLower();
			_alphaStdDevUpper=propParams.alphaStdDevUpper();
			_alphaAcceptTarget=propParams.alphaAcceptTarget();
			_alphaUpdateFreq=propParams.alphaUpdateFreq();
			_alphaAnyUpdates=propParams.alphaAnyUpdates();
			_nTryRho=propParams.nTryRho();
			_nAcceptRho=propParams.nAcceptRho();
			_nLocalAcceptRho=propParams.nLocalAcceptRho();
			_nResetRho=propParams.nResetRho();
			_rhoStdDev=propParams.rhoStdDev();
			_rhoStdDevLower=propParams.rhoStdDevLower();
			_rhoStdDevUpper=propParams.rhoStdDevUpper();
			_rhoAcceptTarget=propParams.rhoAcceptTarget();
			_rhoUpdateFreq=propParams.rhoUpdateFreq();
			_rhoAnyUpdates=propParams.rhoAnyUpdates();
			_nTryLambda=propParams.nTryLambda();
			_nAcceptLambda=propParams.nAcceptLambda();
			_nLocalAcceptLambda=propParams.nLocalAcceptLambda();
			_nResetLambda=propParams.nResetLambda();
			_lambdaStdDev=propParams.lambdaStdDev();
			_lambdaStdDevLower=propParams.lambdaStdDevLower();
			_lambdaStdDevUpper=propParams.lambdaStdDevUpper();
			_lambdaAcceptTarget=propParams.lambdaAcceptTarget();
			_lambdaUpdateFreq=propParams.lambdaUpdateFreq();
			_lambdaAnyUpdates=propParams.lambdaAnyUpdates();
			return *this;

		}

	private:
		unsigned int _nTryTheta;
		unsigned int _nAcceptTheta;
		unsigned int _nLocalAcceptTheta;
		unsigned int _nResetTheta;
		double _thetaStdDev;
		double _thetaStdDevLower;
		double _thetaStdDevUpper;
		double _thetaAcceptTarget;
		unsigned int _thetaUpdateFreq;
		bool _thetaAnyUpdates;
		vector<unsigned int> _nTryBeta;
		vector<unsigned int> _nAcceptBeta;
		vector<unsigned int> _nLocalAcceptBeta;
		vector<unsigned int> _nResetBeta;
		vector<double> _betaStdDev;
		vector<double> _betaStdDevLower;
		vector<double> _betaStdDevUpper;
		double _betaAcceptTarget;
		unsigned int _betaUpdateFreq;
		bool _betaAnyUpdates;
		unsigned int _nTryAlpha;
		unsigned int _nAcceptAlpha;
		unsigned int _nLocalAcceptAlpha;
		unsigned int _nResetAlpha;
		double _alphaStdDev;
		double _alphaStdDevLower;
		double _alphaStdDevUpper;
		double _alphaAcceptTarget;
		unsigned int _alphaUpdateFreq;
		bool _alphaAnyUpdates;
		vector<unsigned int> _nTryRho;
		vector<unsigned int> _nAcceptRho;
		vector<unsigned int> _nLocalAcceptRho;
		vector<unsigned int> _nResetRho;
		vector<double> _rhoStdDev;
		vector<double> _rhoStdDevLower;
		vector<double> _rhoStdDevUpper;
		double _rhoAcceptTarget;
		unsigned int _rhoUpdateFreq;
		bool _rhoAnyUpdates;
		unsigned int _nTryLambda;
		unsigned int _nAcceptLambda;
		unsigned int _nLocalAcceptLambda;
		unsigned int _nResetLambda;
		double _lambdaStdDev;
		double _lambdaStdDevLower;
		double _lambdaStdDevUpper;
		double _lambdaAcceptTarget;
		unsigned int _lambdaUpdateFreq;
		bool _lambdaAnyUpdates;

};

/*********** BLOCK 1 p(v^A,Theta^A,u|.) **********************************/
// A=Active, and Theta contains: phi, mu, Tau, gamma, theta
// We proceed by sampling from p(v^A,Theta^A|.) i.e. marginalising out the u
// and then sampling from p(u|v^A,Theta^A,.). The first of these steps, is achieved
// in a number of stages, updating the various components in turn, and then
// performing label switching.

// Gibbs move for updating the v (for psi) which are active i.e. v_c where c<=Z_max
// This is done by using Gibbs to sample from p(v^A|.), which is the conditional
// for this block with u integrated out (the marginal over u). This can be thought of as
// a MH sample from p(v^A,Theta^A|.), where we sample v^A from conditional, and
// leave Theta^A unchanged (thus having an acceptance rate of 1).
void gibbsForVActive(mcmcChain<diPBaCParams>& chain,
				unsigned int& nTry,unsigned int& nAccept,
				const mcmcModel<diPBaCParams,diPBaCOptions,diPBaCData>& model,
				diPBaCPropParams& propParams,
				baseGeneratorType& rndGenerator){

	mcmcState<diPBaCParams>& currentState = chain.currentState();
	diPBaCParams& currentParams = currentState.parameters();

	nTry++;
	nAccept++;

	// Find the active clusters
	unsigned int maxZ = currentParams.workMaxZi();

	// This is sampled from the posterior given the z vector above
	// Prior comes from the conjugacy of the dirichlet and multinomial
	vector<unsigned int> sumCPlus1ToMaxMembers(maxZ+1);

	sumCPlus1ToMaxMembers[maxZ]=0;
	for(int c=maxZ-1;c>=0;c--){
		sumCPlus1ToMaxMembers[c]=sumCPlus1ToMaxMembers[c+1]+currentParams.workNXInCluster(c+1);
	}

	double tmp=0.0;
	double alpha = currentParams.alpha();
	for(unsigned int c=0;c<=maxZ;c++){
		double vVal = betaRand(rndGenerator,1.0+currentParams.workNXInCluster(c),alpha+sumCPlus1ToMaxMembers[c]);
		currentParams.v(c,vVal);
		// Set psi
		currentParams.logPsi(c,tmp+log(vVal));
		tmp += log(1-vVal);
	}
}

// Moves for updating the Theta which are active i.e. Theta_c where c<=Z_max
// Several different moves, which each update specific elements of Theta

// Move for updating phi
// Gibbs used if no variable selection, or binary switch based variable selection
// Otherwise Metropolis Hastings is used
void updateForPhiActive(mcmcChain<diPBaCParams>& chain,
				unsigned int& nTry,unsigned int& nAccept,
				const mcmcModel<diPBaCParams,diPBaCOptions,diPBaCData>& model,
				diPBaCPropParams& propParams,
				baseGeneratorType& rndGenerator){

	mcmcState<diPBaCParams>& currentState = chain.currentState();
	diPBaCParams& currentParams = currentState.parameters();
	diPBaCHyperParams hyperParams = currentParams.hyperParams();

	const diPBaCData& dataset = model.dataset();
	string varSelectType = model.options().varSelectType();
	// Find the number of clusters
	unsigned int maxZ = currentParams.workMaxZi();
	// Find the number of covariates
	unsigned int nCovariates = currentParams.nCovariates();
	// Find the number of subjects
	unsigned int nSubjects = dataset.nSubjects();

	// Define a uniform random number generator
	randomUniform unifRand(0,1);

	double currentLogPost=0.0;

	for(unsigned int c=0;c<=maxZ;c++){
		// Loop over the covariates
		for(unsigned int j=0;j<nCovariates;j++){
			nTry++;

			if(varSelectType.compare("Continuous")==0){
				currentLogPost=logCondPostPhicj(currentParams,model,c,j);
			}

			unsigned int nCategories = currentParams.nCategories(j);
			// We are updating phis
			// First we must count how many individuals have Xij in each of the
			// possible categories for covariate j.
			vector<double> dirichParams(nCategories,hyperParams.aPhi(j));
			double gammacj = currentParams.gamma(c,j);
			for(unsigned int i=0;i<nSubjects;i++){
				int zi = currentParams.z(i);
				if(zi==(int)c){
					int Xij = dataset.discreteX(i,j);
					// When no variable selection will always add 1.
					// In the binary variable selection case
					// this will add a 1 only when the
					// variable is switched on (as required)
					// In the continuous case this seems like a sensible proposal
					dirichParams[Xij]=dirichParams[Xij]+gammacj;
				}
			}
			vector<double> currentLogPhi(nCategories);
			currentLogPhi=currentParams.logPhi(c,j);

			vector<double> proposedLogPhi(nCategories);
			proposedLogPhi=dirichletRand(rndGenerator,dirichParams);
			for(unsigned int p=0;p<nCategories;p++){
				proposedLogPhi[p]=log(proposedLogPhi[p]);
			}
			currentParams.logPhi(c,j,proposedLogPhi);
			// If no variable selection or binary variable selection
			// this is a sample from full conditional so no accept reject
			// step need. If it is continuous variable selection we
			// need an accept reject decision
			if(varSelectType.compare("Continuous")==0){
				double proposedLogPost = logCondPostPhicj(currentParams,model,c,j);
				double logAcceptRatio=0.0;
				logAcceptRatio=proposedLogPost-currentLogPost;
				logAcceptRatio+=logPdfDirichlet(currentLogPhi,dirichParams,true);
				logAcceptRatio-=logPdfDirichlet(proposedLogPhi,dirichParams,true);
				if(unifRand(rndGenerator)<exp(logAcceptRatio)){
					// Move accepted
					nAccept++;
				}else{
					// Move rejected
					// Reset phi
					currentParams.logPhi(c,j,currentLogPhi);
				}
			}else{
				nAccept++;
			}

		}

	}

}

// Gibbs update for mu in Normal covariate case
void gibbsForMuActive(mcmcChain<diPBaCParams>& chain,
				unsigned int& nTry,unsigned int& nAccept,
				const mcmcModel<diPBaCParams,diPBaCOptions,diPBaCData>& model,
				diPBaCPropParams& propParams,
				baseGeneratorType& rndGenerator){

	mcmcState<diPBaCParams>& currentState = chain.currentState();
	diPBaCParams& currentParams = currentState.parameters();
	diPBaCHyperParams hyperParams = currentParams.hyperParams();

	const diPBaCData& dataset = model.dataset();

	// Find the number of clusters
	unsigned int maxZ = currentParams.workMaxZi();
	// Find the number of covariates
	unsigned int nCovariates = currentParams.nCovariates();
	// Find the number of subjects
	unsigned int nSubjects = dataset.nSubjects();

	nTry++;
	nAccept++;

	// In the following it is useful to have the rows of X as
	// Eigen dynamic vectors
	vector<VectorXd> xi(nSubjects);
	for(unsigned int i=0;i<nSubjects;i++){
		xi[i].setZero(nCovariates);
		for(unsigned int j=0;j<nCovariates;j++){
			xi[i](j)=dataset.continuousX(i,j);
		}
	}

	// We begin by computing the mean X for individuals in each cluster
	vector<VectorXd> meanX(maxZ+1);
	for(unsigned int c=0;c<=maxZ;c++){
		meanX[c].setZero(nCovariates);
	}
	for(unsigned int i=0;i<nSubjects;i++){
		meanX[currentParams.z(i)]=meanX[currentParams.z(i)]+xi[i];
	}

	vector<MatrixXd> gammaMat(maxZ+1);
	vector<MatrixXd> oneMinusGammaMat(maxZ+1);
	for(unsigned int c=0;c<=maxZ;c++){
		gammaMat[c].setZero(nCovariates,nCovariates);
		oneMinusGammaMat[c].setZero(nCovariates,nCovariates);
		for(unsigned int j=0;j<nCovariates;j++){
			gammaMat[c](j,j)=currentParams.gamma(c,j);
			oneMinusGammaMat[c](j,j)=1-gammaMat[c](j,j);
		}
	}

	for(unsigned int c=0;c<=maxZ;c++){
		// Having computed this we can calcuate the posterior mean
		// and posterior covariance for each mu_c
		int nXInC = currentParams.workNXInCluster(c);
		if(nXInC>0){
			meanX[c]=meanX[c]/(double)nXInC;
		}else{
			meanX[c].setZero(nCovariates);
		}

		MatrixXd covMat(nCovariates,nCovariates);
		covMat = (hyperParams.Tau0()+nXInC*gammaMat[c]*currentParams.Tau(c)*gammaMat[c]).inverse();
		VectorXd meanVec(nCovariates);
        meanVec = hyperParams.Tau0()*hyperParams.mu0()+
					nXInC*gammaMat[c]*currentParams.Tau(c)*(meanX[c]-oneMinusGammaMat[c]*currentParams.nullMu());
		meanVec = covMat*meanVec;
		VectorXd mu(nCovariates);
		// We sample from this posterior
		mu = multivarNormalRand(rndGenerator,meanVec,covMat);

		// We store our sample
		currentParams.mu(c,mu);

	}

}

// Gibbs update for Tau in the Normal covariate case
void gibbsForTauActive(mcmcChain<diPBaCParams>& chain,
				unsigned int& nTry,unsigned int& nAccept,
				const mcmcModel<diPBaCParams,diPBaCOptions,diPBaCData>& model,
				diPBaCPropParams& propParams,
				baseGeneratorType& rndGenerator){

	mcmcState<diPBaCParams>& currentState = chain.currentState();
	diPBaCParams& currentParams = currentState.parameters();
	diPBaCHyperParams hyperParams = currentParams.hyperParams();

	const diPBaCData& dataset = model.dataset();

	// Find the number of clusters
	unsigned int maxZ = currentParams.workMaxZi();
	// Find the number of covariates
	unsigned int nCovariates = currentParams.nCovariates();
	// Find the number of subjects
	unsigned int nSubjects = dataset.nSubjects();

	nTry++;
	nAccept++;

	// In the following it is useful to have the rows of X as
	// Eigen dynamic vectors
	vector<VectorXd> xi(nSubjects);
	for(unsigned int i=0;i<nSubjects;i++){
		xi[i].setZero(nCovariates);
		for(unsigned int j=0;j<nCovariates;j++){
			xi[i](j)=dataset.continuousX(i,j);
		}
	}

	vector<MatrixXd> Rc(maxZ+1);
	for(unsigned int c=0;c<=maxZ;c++){
		Rc[c].setZero(nCovariates,nCovariates);
	}

	for(unsigned int i=0;i<nSubjects;i++){
		unsigned int zi = currentParams.z(i);
		Rc[zi]=Rc[zi]+(xi[i]-currentParams.workMuStar(zi))*((xi[i]-currentParams.workMuStar(zi)).transpose());
	}

	for(unsigned int c=0;c<=maxZ;c++){
		Rc[c]=(hyperParams.R0().inverse()+Rc[c]).inverse();
		MatrixXd Tau = wishartRand(rndGenerator,Rc[c],currentParams.workNXInCluster(c)+hyperParams.kappa0());

		currentParams.Tau(c,Tau);

	}


}


// Gibbs update for update of gamma (only used in the binary variable selection case)
void gibbsForGammaActive(mcmcChain<diPBaCParams>& chain,
					unsigned int& nTry,unsigned int& nAccept,
					const mcmcModel<diPBaCParams,
									diPBaCOptions,
									diPBaCData>& model,
					diPBaCPropParams& propParams,
					baseGeneratorType& rndGenerator){


	mcmcState<diPBaCParams>& currentState = chain.currentState();
	diPBaCParams& currentParams = currentState.parameters();
	diPBaCHyperParams hyperParams = currentParams.hyperParams();

	// Find the number of subjects
	unsigned int nCovariates = currentParams.nCovariates();
	unsigned int nSubjects = currentParams.nSubjects();
	unsigned int maxZ = currentParams.workMaxZi();
	string covariateType = model.options().covariateType();
	string varSelectType = model.options().varSelectType();

	// Define a uniform random number generator
	randomUniform unifRand(0,1);

	nTry++;
	nAccept++;

	for(unsigned int j=0;j<nCovariates;j++){

		for(unsigned int c=0;c<=maxZ;c++){

			vector<double> currentGamma=currentParams.gamma(c);
			// work in terms of prob of sticking with current value and switching value
			// Compute probability of sticking
			double logProbStick=0;
			double logProbSwitch=0;
			double probStick=0;
			if(currentParams.omega(j)==0){
				// Nothing to do - not allowed to change
				continue;
			}else{
				for(unsigned int i=0;i<nSubjects;i++){
					unsigned int zi = currentParams.z(i);
					if(zi==c){
						logProbStick+=currentParams.workLogPXiGivenZi(i);
					}
				}
				logProbStick+=(currentGamma[j]*log(currentParams.rho(j))+
							(1-currentGamma[j])*log(1-currentParams.rho(j)));

				// Now compute probability of switching
				currentGamma[j]=1-currentGamma[j];
				currentParams.gamma(c,j,currentGamma[j],covariateType);

				for(unsigned int i=0;i<nSubjects;i++){
					unsigned int zi = currentParams.z(i);
					if(zi==c){
						logProbSwitch+=currentParams.workLogPXiGivenZi(i);
					}
				}
				logProbSwitch+=(currentGamma[j]*log(currentParams.rho(j))+
						(1-currentGamma[j])*log(1-currentParams.rho(j)));

				double maxLogProb;
				if(logProbSwitch<logProbStick){
					maxLogProb=logProbStick;
				}else{
					maxLogProb=logProbSwitch;
				}


				probStick=exp(logProbStick-maxLogProb)/(exp(logProbStick-maxLogProb)+exp(logProbSwitch-maxLogProb));
			}
			if(unifRand(rndGenerator)<probStick){
				// Sticking (we actually need to revert back to what we were
				// before doing the calculations)
				currentGamma[j]=1-currentGamma[j];
				currentParams.gamma(c,j,currentGamma[j],covariateType);
			}
			// Otherwise switching but nothing to do here as we had already done the
			// switch in the calculations

		}
	}

}


// Adaptive Metropolis-Hastings for theta
void metropolisHastingsForThetaActive(mcmcChain<diPBaCParams>& chain,
								unsigned int& nTry,unsigned int& nAccept,
								const mcmcModel<diPBaCParams,
												diPBaCOptions,
												diPBaCData>& model,
								diPBaCPropParams& propParams,
								baseGeneratorType& rndGenerator){

	mcmcState<diPBaCParams>& currentState = chain.currentState();
	diPBaCParams& currentParams = currentState.parameters();
	const string outcomeType = model.dataset().outcomeType();
	unsigned int nCategoriesY = currentParams.nCategoriesY();

	// Find the number of clusters
	unsigned int maxZ = currentParams.workMaxZi();

	// Define a uniform random number generator
	randomUniform unifRand(0,1);
	// Define a normal random number generator
	randomNormal normRand(0,1);

	double thetaTargetRate = propParams.thetaAcceptTarget();
	unsigned int thetaUpdateFreq = propParams.thetaUpdateFreq();

	double currentCondLogPost = logCondPostThetaBeta(currentParams,model);

	for(unsigned int c=0;c<=maxZ;c++){
		for (unsigned int k=0;k<nCategoriesY;k++){
			nTry++;
			propParams.thetaAddTry();
			double& stdDev = propParams.thetaStdDev();
			double thetaOrig = currentParams.theta(c,k);
			double thetaProp = thetaOrig +stdDev*normRand(rndGenerator);
			currentParams.theta(c,k,thetaProp);
			double propCondLogPost = logCondPostThetaBeta(currentParams,model);
			double logAcceptRatio = propCondLogPost - currentCondLogPost;
			if(unifRand(rndGenerator)<exp(logAcceptRatio)){
				nAccept++;
				propParams.thetaAddAccept();
				currentCondLogPost = propCondLogPost;
				// Update the std dev of the proposal
				if(propParams.nTryTheta()%thetaUpdateFreq==0){
					stdDev += 10*(propParams.thetaLocalAcceptRate()-thetaTargetRate)/
								pow((double)(propParams.nTryTheta()/thetaUpdateFreq)+2.0,0.75);
					propParams.thetaAnyUpdates(true);
					if(stdDev>propParams.thetaStdDevUpper()||stdDev<propParams.thetaStdDevLower()){
						propParams.thetaStdDevReset();
					}
					propParams.thetaLocalReset();
				}
			}else{
				currentParams.theta(c,k,thetaOrig);
				// Update the std dev of the proposal
				if(propParams.nTryTheta()%thetaUpdateFreq==0){
					stdDev += 10*(propParams.thetaLocalAcceptRate()-thetaTargetRate)/
								pow((double)(propParams.nTryTheta()/thetaUpdateFreq)+2.0,0.75);
					propParams.thetaAnyUpdates(true);
					if(stdDev<propParams.thetaStdDevLower()||stdDev>propParams.thetaStdDevUpper()){
						propParams.thetaStdDevReset();
					}
					propParams.thetaLocalReset();
				}
			}
		}
	}


}


// Label switching moves (as recommended in Papaspiliopoulos and Roberts, 2008
void metropolisHastingsForLabels(mcmcChain<diPBaCParams>& chain,
				unsigned int& nTry,unsigned int& nAccept,
				const mcmcModel<diPBaCParams,diPBaCOptions,diPBaCData>& model,
				diPBaCPropParams& propParams,
				baseGeneratorType& rndGenerator){

	mcmcState<diPBaCParams>& currentState = chain.currentState();
	diPBaCParams& currentParams = currentState.parameters();

	unsigned int maxZ = currentParams.workMaxZi();
	if(maxZ==0){
		// If there is only one cluster with individuals in, don't do anything
		return;
	}
	string varSelectType = model.options().varSelectType();
	string covariateType = model.options().covariateType();

	randomUniform unifRand(0,1);

	// Move 1 - swap labels of 2 randomly selected non-empty clusters,
	//          leaving psi_c^prop = psi_c for all c

	// Compute how many non-empty clusters
	unsigned int nNotEmpty=0;
	vector<unsigned int> nonEmptyIndices;
	for(unsigned int c=0;c<=maxZ;c++){
		if(currentParams.workNXInCluster(c)>0){
			nNotEmpty++;
			nonEmptyIndices.push_back(c);
		}
	}

	// Select two non-empty clusters at random
	//nTry++;
	unsigned int i1=(unsigned int)nNotEmpty*unifRand(rndGenerator);
	unsigned int c1=nonEmptyIndices[i1];
	nonEmptyIndices.erase(nonEmptyIndices.begin()+i1);
	unsigned int i2=(unsigned int)(nNotEmpty-1)*unifRand(rndGenerator);
	unsigned int c2=nonEmptyIndices[i2];

	// Check whether we accept the move
	double logAcceptRatio = ((double)currentParams.workNXInCluster(c2)-(double)currentParams.workNXInCluster(c1))
								*(currentParams.logPsi(c1)-currentParams.logPsi(c2));

	if(unifRand(rndGenerator)<exp(logAcceptRatio)){
		//nAccept++;
		// Switch the labels
		currentParams.switchLabels(c1,c2,covariateType,varSelectType);
	}

	// Move 2 - swap labels of 2 randomly selected neighbouring clusters,
	//			also swapping the v at the same time
	//nTry++;
	c1=(unsigned int)maxZ*unifRand(rndGenerator);

	logAcceptRatio=(double)currentParams.workNXInCluster(c1)*log(1-currentParams.v(c1+1))
							- (double)currentParams.workNXInCluster(c1+1)*log(1-currentParams.v(c1));

	if(unifRand(rndGenerator)<exp(logAcceptRatio)){
		//nAccept++;
		// Switch the labels
		currentParams.switchLabels(c1,c1+1,covariateType,varSelectType);

		// Also switch the v's
		double v1=currentParams.v(c1);
		double v2=currentParams.v(c1+1);
		double logPsi1=currentParams.logPsi(c1);
		double logPsi2=currentParams.logPsi(c1+1);

		currentParams.logPsi(c1,log(v2)+logPsi1-log(v1));
		currentParams.logPsi(c1+1,log(v1)+log(1-v2)+logPsi2-log(v2)-log(1-v1));
		currentParams.v(c1,v2);
		currentParams.v(c1+1,v1);

		if(c1==maxZ-1){
			// Just check if the maximum Z has now changed to be one less
			if(currentParams.workNXInCluster(c1+1)==0){
				currentParams.workMaxZi(c1);
			}
		}


	}

	// Move 3

	nTry++;
	c1=(unsigned int)maxZ*unifRand(rndGenerator);
	// Compute the acceptance ratio
	unsigned int sumNAfterC1Plus1=0;
	for(unsigned int c=c1+2;c<=maxZ;c++){
		sumNAfterC1Plus1+=currentParams.workNXInCluster(c);
	}
	double const1=0.0,const2=0.0;
	double alpha=currentParams.alpha();
	const1=(1.0+alpha+(double)currentParams.workNXInCluster(c1+1)+(double)sumNAfterC1Plus1)/
			(alpha+(double)currentParams.workNXInCluster(c1+1)+(double)sumNAfterC1Plus1);
	const2=(alpha+(double)currentParams.workNXInCluster(c1)+(double)sumNAfterC1Plus1)/
			(1.0+alpha+(double)currentParams.workNXInCluster(c1)+(double)sumNAfterC1Plus1);
	logAcceptRatio=(double)(currentParams.workNXInCluster(c1)+currentParams.workNXInCluster(c1+1))*
					log(exp(currentParams.logPsi(c1))+exp(currentParams.logPsi(c1+1)));
	logAcceptRatio-=(double)(currentParams.workNXInCluster(c1)+currentParams.workNXInCluster(c1+1))*
					log(exp(currentParams.logPsi(c1))*const1+exp(currentParams.logPsi(c1+1))*const2);
	logAcceptRatio+=(double)(currentParams.workNXInCluster(c1+1))*log(const1);
	logAcceptRatio+=(double)(currentParams.workNXInCluster(c1))*log(const2);

	if(unifRand(rndGenerator)<exp(logAcceptRatio)){
		nAccept++;
		currentParams.switchLabels(c1,c1+1,covariateType,varSelectType);
		double currPsiC1 = exp(currentParams.logPsi(c1));
		double currPsiC1Plus1 = exp(currentParams.logPsi(c1+1));
		double sumCurrPsi = currPsiC1+currPsiC1Plus1;
		double normConst = sumCurrPsi/(const1*currPsiC1Plus1+const2*currPsiC1);
		double propPsiC1 = normConst*const1*currPsiC1Plus1;
		double propPsiC1Plus1 = normConst*const2*currPsiC1;

		double productPrev1MinusV = 1.0;
		if(c1>0){
			productPrev1MinusV = exp(currentParams.logPsi(c1-1))*(1-currentParams.v(c1-1))/
							currentParams.v(c1-1);
		}

		double propVC1=propPsiC1/productPrev1MinusV;
		double propVC1Plus1=propPsiC1Plus1/(productPrev1MinusV*(1-propVC1));

		currentParams.logPsi(c1,log(propPsiC1));
		currentParams.logPsi(c1+1,log(propPsiC1Plus1));
		currentParams.v(c1,propVC1);
		currentParams.v(c1+1,propVC1Plus1);

		if(c1==maxZ-1){
			// Just check if the maximum Z has now changed to be one less
			if(currentParams.workNXInCluster(c1+1)==0){
				currentParams.workMaxZi(c1);
			}
		}
	}

}


// Gibbs move for updating the auxiliary variables u
// This is the second part of block 1. The first part used the marginal
// distribution with u integrated out, we now use the conditional distribution
// for u, conditional on the v^A,Theta^A parameters generated above
// This is done by using Gibbs to sample from p(u|v^A.), which is the conditional

void gibbsForU(mcmcChain<diPBaCParams>& chain,
				unsigned int& nTry,unsigned int& nAccept,
				const mcmcModel<diPBaCParams,diPBaCOptions,diPBaCData>& model,
				diPBaCPropParams& propParams,
				baseGeneratorType& rndGenerator){

	mcmcState<diPBaCParams>& currentState = chain.currentState();
	diPBaCParams& currentParams = currentState.parameters();
	diPBaCHyperParams hyperParams = currentParams.hyperParams();
	string samplerType = model.options().samplerType();

	nTry++;
	nAccept++;

	unsigned int nSubjects = currentParams.nSubjects();
	unsigned int nPredictSubjects = currentParams.nPredictSubjects();

	// Define a uniform random number generator
	randomUniform unifRand(0,1);

	double minUi = 1.0;
	for(unsigned int i=0;i<nSubjects+nPredictSubjects;i++){
		int zi = currentParams.z(i);
		double ui=unifRand(rndGenerator);
		if(samplerType.compare("SliceDependent")==0){
			ui*=exp(currentParams.logPsi((unsigned int)zi));
		}else if(samplerType.compare("SliceIndependent")==0){
			ui*=hyperParams.workXiSlice((unsigned int)zi);
		}

		// This is to avoid numerical errors be
		if(ui<0.0000000001){
			ui=0.0000000001;
		}

		// We only take the minimum over the fitting subjects, because
		// we don't allow predictiction subjects to be allocated to clusters
		// where there are no fitting members, so no need to calc probabilities
		// for cluster which are potential only for predictions subjects
		if(ui<minUi&&i<nSubjects){
			minUi=ui;
		}
		currentParams.u(i,ui);
	}
	currentParams.workMinUi(minUi);

}

/*********** BLOCK 2 p(alpha,v^I|.) **********************************/
// I=Inactive
// We proceed by sampling alpha from p(alpha|.) i.e. marginalising out
// v^I. Then we sample from p(v^I|alpha,.)

// Adaptive Metropolis Hastings move for alpha
void metropolisHastingsForAlpha(mcmcChain<diPBaCParams>& chain,
								unsigned int& nTry,unsigned int& nAccept,
								const mcmcModel<diPBaCParams,
												diPBaCOptions,
												diPBaCData>& model,
								diPBaCPropParams& propParams,
								baseGeneratorType& rndGenerator){

	mcmcState<diPBaCParams>& currentState = chain.currentState();
	diPBaCParams& currentParams = currentState.parameters();
	diPBaCHyperParams hyperParams = currentParams.hyperParams();

	// Find the number of clusters
	unsigned int maxZ = currentParams.workMaxZi();

	// Define a uniform random number generator
	randomUniform unifRand(0,1);

	double& stdDev = propParams.alphaStdDev();
	double alphaCurrent = currentParams.alpha();

	double alphaProp;
	alphaProp = truncNormalRand(rndGenerator,alphaCurrent,stdDev,"L",0,0);

	double logAcceptRatio = 0.0;
	for(unsigned int c=0;c<=maxZ;c++){
		double v=currentParams.v(c);
		logAcceptRatio += logPdfBeta(v,1.0,alphaProp)-logPdfBeta(v,1.0,
																  alphaCurrent);
	}

	// Add in the prior contribution (no contribution if uniform)
	logAcceptRatio += logPdfGamma(alphaProp,hyperParams.shapeAlpha(),hyperParams.rateAlpha());
	logAcceptRatio -= logPdfGamma(alphaCurrent,hyperParams.shapeAlpha(),hyperParams.rateAlpha());

	// Add the proposal contribution
	logAcceptRatio += logPdfTruncatedNormal(alphaCurrent,alphaProp,stdDev,"L",0,0);
	logAcceptRatio -= logPdfTruncatedNormal(alphaProp,alphaCurrent,stdDev,"L",0,0);

	propParams.alphaAddTry();
	nTry++;
	if(unifRand(rndGenerator)<exp(logAcceptRatio)){
		nAccept++;
		propParams.alphaAddAccept();
		// If the move was accepted update the state
		currentParams.alpha(alphaProp);

		// Also update the proposal standard deviation
		if(propParams.nTryAlpha()%propParams.alphaUpdateFreq()==0){
			stdDev += 10*(propParams.alphaLocalAcceptRate()-propParams.alphaAcceptTarget())/
							pow((double)(propParams.nTryAlpha()/propParams.alphaUpdateFreq())+2.0,0.75);
			propParams.alphaAnyUpdates(true);
			if(stdDev>propParams.alphaStdDevUpper()||stdDev<propParams.alphaStdDevLower()){
				propParams.alphaStdDevReset();
			}
			propParams.alphaLocalReset();
		}
	}else{
		// Otherwise update the proposal standard deviation
		if(propParams.nTryAlpha()%propParams.alphaUpdateFreq()==0){
			stdDev += 10*(propParams.alphaLocalAcceptRate()-propParams.alphaAcceptTarget())/
							pow((double)(propParams.nTryAlpha()/propParams.alphaUpdateFreq())+2.0,0.75);
			propParams.alphaAnyUpdates(true);
			if(stdDev>propParams.alphaStdDevUpper()||stdDev<propParams.alphaStdDevLower()){
				propParams.alphaStdDevReset();
			}
			propParams.alphaLocalReset();
		}

	}

}

// Gibbs move for v which are inactive. Only update to maxNClusters, which
// needs to be computed here
void gibbsForVInActive(mcmcChain<diPBaCParams>& chain,
				unsigned int& nTry,unsigned int& nAccept,
				const mcmcModel<diPBaCParams,diPBaCOptions,diPBaCData>& model,
				diPBaCPropParams& propParams,
				baseGeneratorType& rndGenerator){

	mcmcState<diPBaCParams>& currentState = chain.currentState();
	diPBaCParams& currentParams = currentState.parameters();
	diPBaCHyperParams hyperParams = currentParams.hyperParams();
	string samplerType = model.options().samplerType();
	string covariateType = model.options().covariateType();

	nTry++;
	nAccept++;

	unsigned int maxZ = currentParams.workMaxZi();
	unsigned int maxNClusters = currentParams.maxNClusters();
	double minUi = currentParams.workMinUi();

	vector<double> vNew=currentParams.v();
	vector<double> logPsiNew=currentParams.logPsi();

	double alpha = currentParams.alpha();

	if(samplerType.compare("Truncated")==0){
		// Just sample from the prior

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
	}else{
		if(samplerType.compare("SliceIndependent")==0){
			maxNClusters=2+(int)((log(minUi)-log(1.0-hyperParams.rSlice()))/log(hyperParams.rSlice()));
		}

		// Sample V
		vector<double> cumPsi(maxZ+1,0.0);
		cumPsi[0] = exp(currentParams.logPsi(0));
		for(unsigned int c=1;c<=maxZ;c++){
			cumPsi[c]=cumPsi[c-1]+exp(currentParams.logPsi(c));
		}

		bool continueLoop=true;
		unsigned int c=maxZ;
		while(continueLoop){
			if(samplerType.compare("SliceDependent")==0&&cumPsi[c]>1-minUi){
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
		currentParams.maxNClusters(maxNClusters,covariateType);

	}

	currentParams.v(vNew);
	currentParams.logPsi(logPsiNew);


}

/*********** BLOCK 3 p(Theta^I|.) **********************************/
// I=Inactive. Sample the inactive cluster variables from the prior
// Theta contains phi, mu, Tau, gamma, theta. Only need to sample
// up to maxNClusters = max_i{Ci}. Several different routines here for
// each of the variables

// Gibbs move for updating phi
void gibbsForPhiInActive(mcmcChain<diPBaCParams>& chain,
				unsigned int& nTry,unsigned int& nAccept,
				const mcmcModel<diPBaCParams,diPBaCOptions,diPBaCData>& model,
				diPBaCPropParams& propParams,
				baseGeneratorType& rndGenerator){

	mcmcState<diPBaCParams>& currentState = chain.currentState();
	diPBaCParams& currentParams = currentState.parameters();
	diPBaCHyperParams hyperParams = currentParams.hyperParams();

	string varSelectType = model.options().varSelectType();
	// Find the number of clusters
	unsigned int maxZ = currentParams.workMaxZi();
	unsigned int maxNClusters = currentParams.maxNClusters();

	// Find the number of covariates
	unsigned int nCovariates = currentParams.nCovariates();

	nTry++;
	nAccept++;

	for(unsigned int c=maxZ+1;c<maxNClusters;c++){
		// Loop over the covariates
		for(unsigned int j=0;j<nCovariates;j++){
			unsigned int nCategories = currentParams.nCategories(j);
			vector<double> dirichParams(nCategories,hyperParams.aPhi(j));
			vector<double> proposedLogPhi(nCategories);
			proposedLogPhi=dirichletRand(rndGenerator,dirichParams);

			for(unsigned int p=0;p<nCategories;p++){
				proposedLogPhi[p]=log(proposedLogPhi[p]);
			}
			currentParams.logPhi(c,j,proposedLogPhi);
		}
	}


}


// Gibbs update for mu in Normal covariate case
void gibbsForMuInActive(mcmcChain<diPBaCParams>& chain,
				unsigned int& nTry,unsigned int& nAccept,
				const mcmcModel<diPBaCParams,diPBaCOptions,diPBaCData>& model,
				diPBaCPropParams& propParams,
				baseGeneratorType& rndGenerator){

	mcmcState<diPBaCParams>& currentState = chain.currentState();
	diPBaCParams& currentParams = currentState.parameters();
	diPBaCHyperParams hyperParams = currentParams.hyperParams();

	// Find the number of clusters
	unsigned int maxZ = currentParams.workMaxZi();
	unsigned int maxNClusters = currentParams.maxNClusters();
		// Find the number of covariates
	unsigned int nCovariates = currentParams.nCovariates();

	nTry++;
	nAccept++;

	MatrixXd covMat(nCovariates,nCovariates);
	covMat = hyperParams.Tau0().inverse();
	VectorXd meanVec(nCovariates);
	meanVec = hyperParams.mu0();

	for(unsigned int c=maxZ+1;c<maxNClusters;c++){
		VectorXd mu(nCovariates);
		// We sample from this posterior
		mu = multivarNormalRand(rndGenerator,meanVec,covMat);
		// We store our sample
		currentParams.mu(c,mu);
	}

}

// Gibbs update for Tau in the Normal covariate case
void gibbsForTauInActive(mcmcChain<diPBaCParams>& chain,
				unsigned int& nTry,unsigned int& nAccept,
				const mcmcModel<diPBaCParams,diPBaCOptions,diPBaCData>& model,
				diPBaCPropParams& propParams,
				baseGeneratorType& rndGenerator){

	mcmcState<diPBaCParams>& currentState = chain.currentState();
	diPBaCParams& currentParams = currentState.parameters();
	diPBaCHyperParams hyperParams = currentParams.hyperParams();

	// Find the number of clusters
	unsigned int maxZ = currentParams.workMaxZi();
	unsigned int maxNClusters = currentParams.maxNClusters();

	nTry++;
	nAccept++;

	for(unsigned int c=maxZ+1;c<maxNClusters;c++){
		MatrixXd Tau = wishartRand(rndGenerator,hyperParams.R0(),hyperParams.kappa0());
		currentParams.Tau(c,Tau);
	}

}

// Gibbs update for update of gamma (only used in the binary variable selection case)
void gibbsForGammaInActive(mcmcChain<diPBaCParams>& chain,
					unsigned int& nTry,unsigned int& nAccept,
					const mcmcModel<diPBaCParams,
									diPBaCOptions,
									diPBaCData>& model,
					diPBaCPropParams& propParams,
					baseGeneratorType& rndGenerator){


	mcmcState<diPBaCParams>& currentState = chain.currentState();
	diPBaCParams& currentParams = currentState.parameters();
	diPBaCHyperParams hyperParams = currentParams.hyperParams();

	// Find the number of subjects
	unsigned int nCovariates = currentParams.nCovariates();
	unsigned int maxZ = currentParams.workMaxZi();
	unsigned int maxNClusters = currentParams.maxNClusters();
	string covariateType = model.options().covariateType();
	string varSelectType = model.options().varSelectType();

	// Define a uniform random number generator
	randomUniform unifRand(0,1);

	nTry++;
	nAccept++;

	for(unsigned int j=0;j<nCovariates;j++){

		for(unsigned int c=maxZ+1;c<maxNClusters;c++){

			double currentGamma=currentParams.gamma(c,j);
			double proposedGamma=0.0;
			// work in terms of prob of sticking with current value and switching value
			// Compute probability of sticking
			double logProbStick=0;
			double logProbSwitch=0;
			double probSwitch=0;
			if(currentParams.omega(j)==0){
				// Nothing to do - not allowed to change
				continue;
			}else{
				logProbStick+=(currentGamma*log(currentParams.rho(j))+
							(1-currentGamma)*log(1-currentParams.rho(j)));

				// Now compute probability of switching
				proposedGamma=1-currentGamma;

				logProbSwitch+=(proposedGamma*log(currentParams.rho(j))+
						(1-proposedGamma)*log(1-currentParams.rho(j)));

				double maxLogProb;
				if(logProbSwitch<logProbStick){
					maxLogProb=logProbStick;
				}else{
					maxLogProb=logProbSwitch;
				}


				probSwitch=exp(logProbSwitch-maxLogProb)/(exp(logProbStick-maxLogProb)+exp(logProbSwitch-maxLogProb));
			}
			if(unifRand(rndGenerator)<probSwitch){
				// Switching
				currentParams.gamma(c,j,proposedGamma,covariateType);

			}


		}
	}

}

// Gibbs for theta
void gibbsForThetaInActive(mcmcChain<diPBaCParams>& chain,
								unsigned int& nTry,unsigned int& nAccept,
								const mcmcModel<diPBaCParams,
												diPBaCOptions,
												diPBaCData>& model,
								diPBaCPropParams& propParams,
								baseGeneratorType& rndGenerator){

	mcmcState<diPBaCParams>& currentState = chain.currentState();
	diPBaCParams& currentParams = currentState.parameters();
	diPBaCHyperParams hyperParams = currentParams.hyperParams();
	const diPBaCData& dataset = model.dataset();
	unsigned int nCategoriesY=dataset.nCategoriesY();
	const string outcomeType = model.dataset().outcomeType();



	// Find the number of clusters
	unsigned int maxZ = currentParams.workMaxZi();
	unsigned int maxNClusters = currentParams.maxNClusters();

	nTry++;
	nAccept++;

	double location = hyperParams.muTheta();
	double scale = hyperParams.sigmaTheta();
	unsigned int dof = hyperParams.dofTheta();
	randomStudentsT studentsTRand(dof);
	for (unsigned int k=0;k<nCategoriesY;k++){
		for(unsigned int c=maxZ+1;c<maxNClusters;c++){
			double theta=location+scale*studentsTRand(rndGenerator);
			currentParams.theta(c,k,theta);
		}
	}

}


/*********** BLOCK 4 p(Theta^N|.) **********************************/
// N=Non-cluster, and Theta contains: beta, rho, omega, lambda, tau_epsilon

// Adaptive Metropolis-Hastings for beta
void metropolisHastingsForBeta(mcmcChain<diPBaCParams>& chain,
								unsigned int& nTry,unsigned int& nAccept,
								const mcmcModel<diPBaCParams,
												diPBaCOptions,
												diPBaCData>& model,
								diPBaCPropParams& propParams,
								baseGeneratorType& rndGenerator){

	mcmcState<diPBaCParams>& currentState = chain.currentState();
	diPBaCParams& currentParams = currentState.parameters();
	const string outcomeType = model.dataset().outcomeType();

	// Find the number of clusters
	unsigned int nFixedEffects = currentParams.nFixedEffects(outcomeType);

	// Find the number of categories of Y
	unsigned int nCategoriesY = currentParams.nCategoriesY();

	// Define a uniform random number generator
	randomUniform unifRand(0,1);
	// Define a normal random number generator
	randomNormal normRand(0,1);

	double betaTargetRate = propParams.betaAcceptTarget();
	unsigned int betaUpdateFreq = propParams.betaUpdateFreq();

	double currentCondLogPost = logCondPostThetaBeta(currentParams,model);

	for(unsigned int j=0;j<nFixedEffects;j++){
		for (unsigned int k=0;k<nCategoriesY;k++){
			nTry++;
			propParams.betaAddTry(j);
			double& stdDev = propParams.betaStdDev(j);
			double betaOrig = currentParams.beta(j,k);
			double betaProp = betaOrig+stdDev*normRand(rndGenerator);
			currentParams.beta(j,k,betaProp);
			double propCondLogPost = logCondPostThetaBeta(currentParams,model);
			double logAcceptRatio = propCondLogPost - currentCondLogPost;
			if(unifRand(rndGenerator)<exp(logAcceptRatio)){
				nAccept++;
				propParams.betaAddAccept(j);
				currentCondLogPost = propCondLogPost;
				// Update the std dev of the proposal
				if(propParams.nTryBeta(j)%betaUpdateFreq==0){
					stdDev += 10*(propParams.betaLocalAcceptRate(j)-betaTargetRate)/
							pow((double)(propParams.nTryBeta(j)/betaUpdateFreq)+2.0,0.75);
					propParams.betaAnyUpdates(true);
					if(stdDev>propParams.betaStdDevUpper(j)||stdDev<propParams.betaStdDevLower(j)){
						propParams.betaStdDevReset(j);
					}
					propParams.betaLocalReset(j);
				}
			}else{
				currentParams.beta(j,k,betaOrig);
				// Update the std dev of the proposal
				if(propParams.nTryBeta(j)%betaUpdateFreq==0){
					stdDev += 10*(propParams.betaLocalAcceptRate(j)-betaTargetRate)/
							pow((double)(propParams.nTryBeta(j)/betaUpdateFreq)+2.0,0.75);
					propParams.betaAnyUpdates(true);
					if(stdDev<propParams.betaStdDevLower(j)||stdDev>propParams.betaStdDevUpper(j)){
						propParams.betaStdDevReset(j);
					}
					propParams.betaLocalReset(j);
				}
			}
		}
	}
}


// Adaptive Metropolis-Hastings for lambda
void metropolisHastingsForLambda(mcmcChain<diPBaCParams>& chain,
								unsigned int& nTry,unsigned int& nAccept,
								const mcmcModel<diPBaCParams,
												diPBaCOptions,
												diPBaCData>& model,
								diPBaCPropParams& propParams,
								baseGeneratorType& rndGenerator){

	mcmcState<diPBaCParams>& currentState = chain.currentState();
	diPBaCParams& currentParams = currentState.parameters();
	const string outcomeType = model.dataset().outcomeType();

	// Find the number of subjects
	unsigned int nSubjects = currentParams.nSubjects();

	// Define a uniform random number generator
	randomUniform unifRand(0,1);
	// Define a normal random number generator
	randomNormal normRand(0,1);


	double lambdaTargetRate = propParams.lambdaAcceptTarget();
	unsigned int lambdaUpdateFreq = propParams.lambdaUpdateFreq();

	double (*logCondPostLambdai)(const diPBaCParams&,
											const mcmcModel<diPBaCParams,
															diPBaCOptions,
															diPBaCData>&,
											const unsigned int&) = NULL;

	if(outcomeType.compare("Bernoulli")==0){
		logCondPostLambdai = &logCondPostLambdaiBernoulli;
	}else if(outcomeType.compare("Binomial")==0){
		logCondPostLambdai = &logCondPostLambdaiBinomial;
	}else if(outcomeType.compare("Poisson")==0){
		logCondPostLambdai = &logCondPostLambdaiPoisson;
	}

	for(unsigned int i=0;i<nSubjects;i++){
		// Only update each lambda i with probability 0.1
		if(unifRand(rndGenerator)>0.1){
			continue;
		}

		nTry++;
		propParams.lambdaAddTry();

		double currentCondLogPost = logCondPostLambdai(currentParams,model,i);
		double& stdDev = propParams.lambdaStdDev();
		double lambdaOrig = currentParams.lambda(i);
		double lambdaProp = lambdaOrig+stdDev*normRand(rndGenerator);
		currentParams.lambda(i,lambdaProp);
		double propCondLogPost = logCondPostLambdai(currentParams,model,i);
		double logAcceptRatio = propCondLogPost - currentCondLogPost;
		if(unifRand(rndGenerator)<exp(logAcceptRatio)){
			nAccept++;
			propParams.lambdaAddAccept();
			// Update the std dev of the proposal
			if(propParams.nTryLambda()%lambdaUpdateFreq==0){
				stdDev += 10*(propParams.lambdaLocalAcceptRate()-lambdaTargetRate)/
								pow((double)(propParams.nTryLambda()/lambdaUpdateFreq)+2.0,0.75);
				propParams.lambdaAnyUpdates(true);
				if(stdDev>propParams.lambdaStdDevUpper()||stdDev<propParams.lambdaStdDevLower()){
					propParams.lambdaStdDevReset();
				}
				propParams.lambdaLocalReset();
			}
		}else{
			currentParams.lambda(i,lambdaOrig);
			// Update the std dev of the proposal
			if(propParams.nTryLambda()%lambdaUpdateFreq==0){
				stdDev += 10*(propParams.lambdaLocalAcceptRate()-lambdaTargetRate)/
								pow((double)(propParams.nTryLambda()/lambdaUpdateFreq)+2.0,0.75);
				propParams.lambdaAnyUpdates(true);
				if(stdDev<propParams.lambdaStdDevLower()||stdDev>propParams.lambdaStdDevUpper()){
					propParams.lambdaStdDevReset();
				}
				propParams.lambdaLocalReset();
			}
		}
	}

}

// Gibbs update for the precision of extra variation epsilon
void gibbsForTauEpsilon(mcmcChain<diPBaCParams>& chain,
						unsigned int& nTry,unsigned int& nAccept,
						const mcmcModel<diPBaCParams,
										diPBaCOptions,
										diPBaCData>& model,
						diPBaCPropParams& propParams,
						baseGeneratorType& rndGenerator){

	mcmcState<diPBaCParams>& currentState = chain.currentState();
	diPBaCParams& currentParams = currentState.parameters();
	diPBaCHyperParams hyperParams = currentParams.hyperParams();
	const diPBaCData& dataset = model.dataset();
	const string& outcomeType = model.dataset().outcomeType();

	unsigned int nSubjects=dataset.nSubjects();
	unsigned int nFixedEffects=dataset.nFixedEffects();

	nTry++;
	nAccept++;

	double a=hyperParams.shapeTauEpsilon(),b=hyperParams.rateTauEpsilon();

	double sumEpsilon = 0.0;
	vector<double> meanVec(nSubjects,0.0);
	if(outcomeType.compare("Poisson")==0){
		meanVec=dataset.logOffset();
	}
	for(unsigned int i=0;i<nSubjects;i++){
		int zi=currentParams.z(i);
		double meanVal=meanVec[i]+currentParams.theta(zi,0);
		for(unsigned int j=0;j<nFixedEffects;j++){
			meanVal+=currentParams.beta(j,0)*dataset.W(i,j);
		}
		double eps = currentParams.lambda(i)-meanVal;
		sumEpsilon+=pow(eps,2.0);
	}
	a+=(double)nSubjects/2.0;
	b+=sumEpsilon/2.0;

	// Boost uses shape and scale parameterisation
	randomGamma gammaRand(a,1.0/b);
	double tau = gammaRand(rndGenerator);
	currentParams.tauEpsilon(tau);

}

// Metropolis-Hastings for joint uptdate of rho and omega
void metropolisHastingsForRhoOmega(mcmcChain<diPBaCParams>& chain,
									unsigned int& nTry,unsigned int& nAccept,
									const mcmcModel<diPBaCParams,
													diPBaCOptions,
													diPBaCData>& model,
									diPBaCPropParams& propParams,
									baseGeneratorType& rndGenerator){

	mcmcState<diPBaCParams>& currentState = chain.currentState();
	diPBaCParams& currentParams = currentState.parameters();
	diPBaCHyperParams hyperParams = currentParams.hyperParams();

	// Find the number of subjects
	unsigned int nCovariates = currentParams.nCovariates();
	string covariateType = model.options().covariateType();
	string varSelectType = model.options().varSelectType();

	// Define a uniform random number generator
	randomUniform unifRand(0,1);
	// Define a normal random number generator
	randomNormal normRand(0,1);

	double currentLogPost = 0;
	if(varSelectType.compare("Continuous")==0){
		currentLogPost = logCondPostRhoOmegaj(currentParams,model,0);
	}


	double proposedLogPost = currentLogPost;
	vector<unsigned int> currentOmega = currentParams.omega();
	unsigned int proposedOmega;
	vector<double> currentRho = currentParams.rho();
	double proposedRho;

	for(unsigned int j=0;j<nCovariates;j++){

		if(varSelectType.compare("Continuous")!=0){
			currentLogPost = logCondPostRhoOmegaj(currentParams,model,j);
		}

		nTry++;

		// Propose from the priors
		double& stdDev = propParams.rhoStdDev(j);
		if(unifRand(rndGenerator)<0.5){
			// Proposing an omega 0
			if(currentOmega[j]==0){
				// Nothing to do as move to the same place
				nAccept++;
				continue;
			}
			proposedOmega=0;
			proposedRho=0.0;

			currentParams.omega(j,proposedOmega);
			currentParams.rho(j,proposedRho,covariateType,varSelectType);
			proposedLogPost = logCondPostRhoOmegaj(currentParams,model,j);
			double logAcceptRatio = proposedLogPost - currentLogPost;
			logAcceptRatio += logPdfBeta(currentRho[j],hyperParams.aRho(),hyperParams.bRho());
			if(unifRand(rndGenerator)<exp(logAcceptRatio)){
				// Move accepted
				if(varSelectType.compare("Continuous")==0){
					currentLogPost=proposedLogPost;
				}
				nAccept++;
			}else{
				// Move rejected, reset parameters
				currentParams.omega(j,currentOmega[j]);
				currentParams.rho(j,currentRho[j],covariateType,varSelectType);
			}

		}else{
			if(currentOmega[j]==1){
				proposedRho  = truncNormalRand(rndGenerator,currentRho[j],stdDev,"B",0,1);
				currentParams.rho(j,proposedRho,covariateType,varSelectType);
				proposedLogPost = logCondPostRhoOmegaj(currentParams,model,j);
				double logAcceptRatio = proposedLogPost - currentLogPost;
				logAcceptRatio += logPdfTruncatedNormal(currentRho[j],proposedRho,stdDev,"B",0,1);
				logAcceptRatio -= logPdfTruncatedNormal(proposedRho,currentRho[j],stdDev,"B",0,1);
				propParams.rhoAddTry(j);
				if(unifRand(rndGenerator)<exp(logAcceptRatio)){
					// Move accepted
					if(varSelectType.compare("Continuous")==0){
						currentLogPost=proposedLogPost;
					}

					nAccept++;
					propParams.rhoAddAccept(j);
					// Also update the proposal standard deviation
					if(propParams.nTryRho(j)%propParams.rhoUpdateFreq()==0){
						stdDev += 0.1*(propParams.rhoLocalAcceptRate(j)-propParams.rhoAcceptTarget())/
										pow((double)(propParams.nTryRho(j)/propParams.rhoUpdateFreq())+2.0,0.75);
						propParams.rhoAnyUpdates(true);
						if(stdDev>propParams.rhoStdDevUpper(j)||stdDev<propParams.rhoStdDevLower(j)){
							propParams.rhoStdDevReset(j);
						}
						propParams.rhoLocalReset(j);
					}
				}else{
					// Move rejected, reset parameters
					currentParams.omega(j,currentOmega[j]);
					currentParams.rho(j,currentRho[j],covariateType,varSelectType);
					// Also update the proposal standard deviation
					if(propParams.nTryRho(j)%propParams.rhoUpdateFreq()==0){
						stdDev += 0.1*(propParams.rhoLocalAcceptRate(j)-propParams.rhoAcceptTarget())/
										pow((double)(propParams.nTryRho(j)/propParams.rhoUpdateFreq())+2.0,0.75);
						propParams.rhoAnyUpdates(true);
					    if(stdDev>propParams.rhoStdDevUpper(j)||stdDev<propParams.rhoStdDevLower(j)){
							propParams.rhoStdDevReset(j);
						}
						propParams.rhoLocalReset(j);
					}
				}
			}else{
				proposedRho = betaRand(rndGenerator,hyperParams.aRho(),hyperParams.bRho());
				proposedOmega=1;
				currentParams.omega(j,proposedOmega);
				currentParams.rho(j,proposedRho,covariateType,varSelectType);
				proposedLogPost = logCondPostRhoOmegaj(currentParams,model,j);
				double logAcceptRatio = proposedLogPost - currentLogPost;
				logAcceptRatio -= logPdfBeta(proposedRho,hyperParams.aRho(),hyperParams.bRho());
				if(unifRand(rndGenerator)<exp(logAcceptRatio)){
					// Move accepted
					if(varSelectType.compare("Continuous")==0){
						currentLogPost=proposedLogPost;
					}
					nAccept++;
				}else{
					// Move rejected, reset parameters
					currentParams.omega(j,currentOmega[j]);
					currentParams.rho(j,currentRho[j],covariateType,varSelectType);
				}
			}
		}


	}

}

// Gibbs for update of sigmaSqY (Normal response case)
void gibbsForSigmaSqY(mcmcChain<diPBaCParams>& chain,
						unsigned int& nTry,unsigned int& nAccept,
						const mcmcModel<diPBaCParams,
										diPBaCOptions,
										diPBaCData>& model,
						diPBaCPropParams& propParams,
						baseGeneratorType& rndGenerator){

	mcmcState<diPBaCParams>& currentState = chain.currentState();
	diPBaCParams& currentParams = currentState.parameters();
	diPBaCHyperParams hyperParams = currentParams.hyperParams();

	const diPBaCData& dataset = model.dataset();

	unsigned int nSubjects = currentParams.nSubjects();
	unsigned int nFixedEffects = dataset.nFixedEffects();

	nTry++;
	nAccept++;

	double sumVal=0.0;
	for(unsigned int i=0;i<nSubjects;i++){
		int Zi = currentParams.z(i);

		double mu = currentParams.theta(Zi,0);
		for(unsigned int j=0;j<nFixedEffects;j++){
			mu+=currentParams.beta(j,0)*dataset.W(i,j);
		}

		sumVal+=pow(dataset.continuousY(i)-mu,2.0);
	}

	double posteriorShape = hyperParams.shapeSigmaSqY()+(double)nSubjects/2.0;
	double posteriorScale = hyperParams.scaleSigmaSqY()+0.5*sumVal;

	randomGamma gammaRand(posteriorShape,1.0/posteriorScale);
	double sigmaSqY=1.0/gammaRand(rndGenerator);
	currentParams.sigmaSqY(sigmaSqY);


}
/*********** BLOCK 5 p(Z|.) **********************************/

// Gibbs update for the allocation variables
void gibbsForZ(mcmcChain<diPBaCParams>& chain,
				unsigned int& nTry,unsigned int& nAccept,
				const mcmcModel<diPBaCParams,
								diPBaCOptions,
								diPBaCData>& model,
				diPBaCPropParams& propParams,
				baseGeneratorType& rndGenerator){

	mcmcState<diPBaCParams>& currentState = chain.currentState();
	diPBaCParams& currentParams = currentState.parameters();
	diPBaCHyperParams hyperParams = currentParams.hyperParams();
	const diPBaCData& dataset = model.dataset();
	const string& outcomeType = model.dataset().outcomeType();
	const string& covariateType = model.dataset().covariateType();
	const string& samplerType = model.options().samplerType();
	bool computeEntropy = model.options().computeEntropy();
	unsigned int nSubjects=dataset.nSubjects();
	unsigned int nPredictSubjects=dataset.nPredictSubjects();
	unsigned int maxNClusters=currentParams.maxNClusters();
	unsigned int nFixedEffects=dataset.nFixedEffects();
	unsigned int nCategoriesY=dataset.nCategoriesY();
	unsigned int nCovariates=dataset.nCovariates();
	vector<unsigned int>nCategories=dataset.nCategories();
	const vector<vector<bool> >& missingX=dataset.missingX();
	bool includeResponse = model.options().includeResponse();
	bool responseExtraVar = model.options().responseExtraVar();

	nTry++;
	nAccept++;

	// Define a uniform random number generator
	randomUniform unifRand(0,1);

	vector<unsigned int> nMembers(maxNClusters,0);

	vector<double> rnd(nSubjects+nPredictSubjects,0.0);
	vector<double> u(nSubjects+nPredictSubjects,0.0);
	for(unsigned int i=0;i<nSubjects+nPredictSubjects;i++){
		rnd[i] = unifRand(rndGenerator);
		u[i] = currentParams.u(i);
	}

	vector<double> testBound(maxNClusters,0.0);
	vector<double> clusterWeight(maxNClusters,0.0);
	for(unsigned int c=0;c<maxNClusters;c++){
		if(samplerType.compare("SliceDependent")==0){
			testBound[c] = exp(currentParams.logPsi(c));
			clusterWeight[c] = 0.0;
		}else if(samplerType.compare("SliceIndependent")==0){
			testBound[c] = hyperParams.workXiSlice(c);
			clusterWeight[c] = currentParams.logPsi(c)-(double)c*log(hyperParams.rSlice())-log(1-hyperParams.rSlice());
		}else if(samplerType.compare("Truncated")==0){
			testBound[c] = 1.0;
			clusterWeight[c] = currentParams.logPsi(c);
		}
	}

	// Compute the allocation probabilities in terms of the unique vectors
	vector<vector<double> > logPXiGivenZi;
	logPXiGivenZi.resize(nSubjects+nPredictSubjects);


	if(covariateType.compare("Discrete")==0){
		for(unsigned int i=0;i<nSubjects;i++){
			logPXiGivenZi[i].resize(maxNClusters,0);
			for(unsigned int c=0;c<maxNClusters;c++){
				if(u[i]<testBound[c]){
					if(currentParams.z(i)==(int)c){
						logPXiGivenZi[i][c]=currentParams.workLogPXiGivenZi(i);
					}else{
						logPXiGivenZi[i][c]=0;
						for(unsigned int j=0;j<nCovariates;j++){
							int Xij = currentParams.workDiscreteX(i,j);
							logPXiGivenZi[i][c]+=currentParams.workLogPhiStar(c,j,(unsigned int)Xij);
						}
					}
				}
			}

		}
		// For the predictive subjects we do not count missing data
		for(unsigned int i=nSubjects;i<nSubjects+nPredictSubjects;i++){
			logPXiGivenZi[i].resize(maxNClusters,0);
			for(unsigned int c=0;c<maxNClusters;c++){
				if(u[i]<testBound[c]){
					for(unsigned int j=0;j<nCovariates;j++){
						if(!missingX[i][j]){
							int Xij = currentParams.workDiscreteX(i,j);
							logPXiGivenZi[i][c]+=currentParams.workLogPhiStar(c,j,(unsigned int)Xij);
						}
					}
				}
			}
		}

	}else if(covariateType.compare("Normal")==0){
		for(unsigned int i=0;i<nSubjects;i++){
			logPXiGivenZi[i].resize(maxNClusters,0.0);
			VectorXd xi=VectorXd::Zero(nCovariates);
			for(unsigned int c=0;c<maxNClusters;c++){
				if(u[i]<testBound[c]){
					if(currentParams.z(i)==(int)c){
						logPXiGivenZi[i][c]=currentParams.workLogPXiGivenZi(i);
					}else{
						for(unsigned int j=0;j<nCovariates;j++){
							xi(j)=currentParams.workContinuousX(i,j);
						}
						logPXiGivenZi[i][c]=logPdfMultivarNormal(nCovariates,xi,currentParams.workMuStar(c),currentParams.workSqrtTau(c),currentParams.workLogDetTau(c));
					}
				}

			}
		}
		// For the predictive subjects we do not count missing data
		LLT<MatrixXd> llt;
		for(unsigned int i=nSubjects;i<nSubjects+nPredictSubjects;i++){
			logPXiGivenZi[i].resize(maxNClusters,0.0);

			unsigned int nNotMissing=dataset.nCovariatesNotMissing(i);

			for(unsigned int c=0;c<maxNClusters;c++){
				if(u[i]<testBound[c]){
					VectorXd workMuStar=currentParams.workMuStar(c);

					VectorXd xi=VectorXd::Zero(nNotMissing);
					VectorXd muStar=VectorXd::Zero(nNotMissing);
					MatrixXd sqrtTau=MatrixXd::Zero(nNotMissing,nNotMissing);
					double logDetTau =0.0;
					if(nNotMissing==nCovariates){
						muStar=workMuStar;
						sqrtTau=currentParams.workSqrtTau(c);
						logDetTau=currentParams.workLogDetTau(c);
						for(unsigned int j=0;j<nCovariates;j++){
							xi(j)=currentParams.workContinuousX(i,j);
						}
					}else{
						MatrixXd workSigma=currentParams.Sigma(c);
						MatrixXd Sigma=MatrixXd::Zero(nNotMissing,nNotMissing);
						MatrixXd Tau=MatrixXd::Zero(nNotMissing,nNotMissing);
						unsigned int j=0;
						for(unsigned int j0=0;j0<nCovariates;j0++){
							if(!missingX[i][j0]){
								xi(j)=currentParams.workContinuousX(i,j0);
								muStar(j)=workMuStar(j0);
								unsigned int r=0;
								for(unsigned int j1=0;j1<nCovariates;j1++){
									if(!missingX[i][j1]){
										Sigma(j,r)=workSigma(j0,j1);
										r++;
									}
								}
								j++;
							}
						}
						Tau = Sigma.inverse();
						sqrtTau = (llt.compute(Tau)).matrixU();
						logDetTau = log(Tau.determinant());

					}
					logPXiGivenZi[i][c]=logPdfMultivarNormal(nNotMissing,xi,muStar,sqrtTau,logDetTau);
				}
			}
		}

	}

	double (*logPYiGivenZiWi)(const diPBaCParams&, const diPBaCData&,
											const unsigned int&,const int&,
											const unsigned int&)=NULL;

	if(includeResponse){
		if(!responseExtraVar){
			if(outcomeType.compare("Bernoulli")==0){
				logPYiGivenZiWi = &logPYiGivenZiWiBernoulli;
			}else if(outcomeType.compare("Binomial")==0){
				logPYiGivenZiWi = &logPYiGivenZiWiBinomial;
			}else if(outcomeType.compare("Poisson")==0){
				logPYiGivenZiWi = &logPYiGivenZiWiPoisson;
			}else if(outcomeType.compare("Normal")==0){
				logPYiGivenZiWi = &logPYiGivenZiWiNormal;
			}else if(outcomeType.compare("Categorical")==0){
				logPYiGivenZiWi = &logPYiGivenZiWiCategorical;
			}
		}
	}

	vector<double> meanVec(nSubjects,0.0);
	if(includeResponse){
		if(outcomeType.compare("Poisson")==0){
			meanVec = dataset.logOffset();
		}
	}

	unsigned int maxZ=0;
	for(unsigned int i=0;i<nSubjects+nPredictSubjects;i++){

		vector<double> logPyXz(maxNClusters,0.0);
		// We calculate p(Z=c|y,X) \propto p(y,X,Z=c)
		// p(y,X,z=c) = p(y|Z=c)p(X|z=c)p(z=c)
		double maxLogPyXz = -(numeric_limits<double>::max());

		if(includeResponse&&i<nSubjects){
			// Response only included in allocation probs for fitting subjects, not predicting subjects
			if(responseExtraVar){
				// In this case the Y only go into this conditional
				// through lambda
				for(unsigned int c=0;c<maxNClusters;c++){
					if(u[i]<testBound[c]){
						double meanVal=meanVec[i]+currentParams.theta(c,0);
						for(unsigned int j=0;j<nFixedEffects;j++){
							meanVal+=currentParams.beta(j,0)*dataset.W(i,j);
						}

						logPyXz[c]+=logPdfNormal(currentParams.lambda(i),meanVal,
													1/sqrt(currentParams.tauEpsilon()));
					}
				}
			}else{
				// In this case the Y go in directly
				for(unsigned int c=0;c<maxNClusters;c++){
					if(u[i]<testBound[c]){
						logPyXz[c]+=logPYiGivenZiWi(currentParams,dataset,nFixedEffects,c,i);
					}
				}
			}
		}

		for(unsigned int c=0;c<maxNClusters;c++){
			if(u[i]<testBound[c]){
				// Make sure prediction subjects can only be allocated to one
				// of the non-empty clusters
				if(i<nSubjects||nMembers[c]>0){
					logPyXz[c]+=clusterWeight[c];
					logPyXz[c]+=logPXiGivenZi[i][c];
				}else{
					logPyXz[c]=-(numeric_limits<double>::max());
				}
			}else{
				logPyXz[c]=-(numeric_limits<double>::max());
			}
			if(logPyXz[c]>maxLogPyXz){
				maxLogPyXz=logPyXz[c];
			}
		}
		vector<double> pzGivenXy(maxNClusters);
		double sumVal=0;
		for(unsigned int c=0;c<maxNClusters;c++){
			double exponent = logPyXz[c] - maxLogPyXz;
			// Check for negative infinity (can only be negative)
			if(std::isinf(exponent)||std::isnan(exponent)){
				exponent=-(numeric_limits<double>::max());
			}
			pzGivenXy[c]=exp(exponent);
			sumVal+=pzGivenXy[c];
		}

		vector<double> expectedTheta(nCategoriesY);
		double entropyVal=0.0;
		vector<double> cumPzGivenXy(maxNClusters);
		for(unsigned int c=0;c<maxNClusters;c++){
			pzGivenXy[c]/=sumVal;
			if(computeEntropy){
				if(pzGivenXy[c]>0){
					entropyVal-=pzGivenXy[c]*log(pzGivenXy[c]);
				}
			}

			if(c==0){
				cumPzGivenXy[c]=pzGivenXy[c];
			}else{
				cumPzGivenXy[c]=cumPzGivenXy[c-1]+pzGivenXy[c];
			}
			if(includeResponse&&i>=nSubjects){
				if(outcomeType.compare("Categorical")==0){
					for (unsigned int k=0;k<nCategoriesY;k++){
						expectedTheta[k]+=currentParams.theta(c,k)*pzGivenXy[c];
					}
				} else {
					expectedTheta[0]+=currentParams.theta(c,0)*pzGivenXy[c];
				}
			}
		}
		unsigned int zi;
		if(maxNClusters==1){
			zi=0;
		}else{
			zi = 0;
			for(unsigned int c=0;c<maxNClusters;c++){
				if(rnd[i]<cumPzGivenXy[c]){
					zi=c;
					break;
				}
			}
		}
		if(zi>maxZ){
			maxZ=zi;
		}


		currentParams.z(i,zi,covariateType);
		if(computeEntropy){
			currentParams.workEntropy(i,entropyVal);
		}
		if(i<nSubjects){
			nMembers[zi]++;
		}else{
			if(includeResponse){
				for (unsigned int k=0;k<nCategoriesY;k++){
					currentParams.workPredictExpectedTheta(i-nSubjects,k,expectedTheta[k]);
				}
			}
		}
	}

	currentParams.workNXInCluster(nMembers);
	currentParams.workMaxZi(maxZ);

}


void updateMissingDiPBaCData(baseGeneratorType& rndGenerator,
								diPBaCParams& params,
								const diPBaCOptions& options,
								diPBaCData& dataset){

	unsigned int nSubjects = dataset.nSubjects();
	unsigned int nCovariates = dataset.nCovariates();
	vector<unsigned int> nCategories = params.nCategories();
	string covariateType = options.covariateType();

	// Define a uniform random number generator
	randomUniform unifRand(0,1);

	// For the missing variables we check which cluster the subject is allocated
	// to and then sample for the appropriate
	if(covariateType.compare("Discrete")==0){
		// We don't update the predictive subjects as their X values which
		// were missing are not used anywhere
		for(unsigned int i=0;i<nSubjects;i++){
			int zi = params.z(i);
			for(unsigned int j=0;j<nCovariates;j++){
				if(dataset.missingX(i,j)){
					vector<double> cumPhiStar(nCategories[j]);
					// Sample uniform
					double u=unifRand(rndGenerator);
					int k=0;
					cumPhiStar[k]=exp(params.workLogPhiStar(zi,j,k));
					while(u>cumPhiStar[k]){
						k++;
						// Make phi cumulative
						cumPhiStar[k]=cumPhiStar[k-1]+exp(params.workLogPhiStar(zi,j,k));
					}

					int prevX = params.workDiscreteX(i,j);
					dataset.discreteX(i,j,k);
					params.workDiscreteX(i,j,k);
					// Now we need to recompute the workLogPXiGivenZi values based
					// on the new X
					if(prevX!=k){
						double logVal = params.workLogPXiGivenZi(i);
						double oldVal, newVal;
						oldVal = params.workLogPhiStar(zi,j,prevX);
						newVal = params.workLogPhiStar(zi,j,k);
						logVal+=(newVal-oldVal);
						params.workLogPXiGivenZi(i,logVal);
					}
				}
			}
		}
	}else if(covariateType.compare("Normal")==0){
		for(unsigned int i=0;i<nSubjects;i++){
			// Check if there is anything to do
			if(dataset.nCovariatesNotMissing(i)<nCovariates){
				int zi = params.z(i);
				VectorXd newXi=multivarNormalRand(rndGenerator,params.workMuStar(zi),params.Sigma(zi));
				for(unsigned int j=0;j<nCovariates;j++){
					if(dataset.missingX(i,j)){
						dataset.continuousX(i,j,newXi(j));
						params.workContinuousX(i,j,newXi(j));
					}else{
						newXi(j)=dataset.continuousX(i,j);
					}
				}
				double logVal = logPdfMultivarNormal(nCovariates,newXi,params.workMuStar(zi),params.workSqrtTau(zi),params.workLogDetTau(zi));
				params.workLogPXiGivenZi(i,logVal);
			}
		}

	}


}


#endif /* DIPBACPROPOSALS_H_ */
