/// \file PReMiuMProposals.h
/// \author David Hastie
/// \brief Header file for model specification for PReMiuMpp

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
#include<PReMiuMOptions.h>
#include<PReMiuMModel.h>
#include<PReMiuMData.h>
#include<PReMiuMArs.h>

using namespace Eigen;

using std::vector;
using std::accumulate;
using std::numeric_limits;
using boost::math::normal_distribution;
using boost::math::students_t_distribution;
using boost::math::lgamma;


class pReMiuMPropParams{

	public:
		// Default constructor
		pReMiuMPropParams() {};

		pReMiuMPropParams(const unsigned int& nSweeps,const unsigned int& nCovariates,
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

				_kappa1StdDev = 1.0;
				_kappa1StdDevLower = 0.1;
				_kappa1StdDevUpper = 99.9;
				_nTryKappa1 = 0;
				_nAcceptKappa1 = 0;
				_nLocalAcceptKappa1 = 0;
				_nResetKappa1 = 0;
				_kappa1AcceptTarget = 0.44;
				_kappa1UpdateFreq = 25;
				_kappa1AnyUpdates = true; 

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

				_nTryTauS.resize(nCovariates);
				_nAcceptTauS.resize(nCovariates);
				_nLocalAcceptTauS.resize(nCovariates);
				_nResetTauS.resize(nCovariates);
				_TauSStdDev.resize(nCovariates);
				_TauSStdDevLower.resize(nCovariates);
				_TauSStdDevUpper.resize(nCovariates);
				for (unsigned int j = 0; j<nCovariates; j++) {
					_TauSStdDev[j] = 1.0;
					_TauSStdDevLower[j] = 0.1;
					_TauSStdDevUpper[j] = 99.9;
					_nTryTauS[j] = 0;
					_nAcceptTauS[j] = 0;
					_nLocalAcceptTauS[j] = 0;
					_nResetTauS[j] = 0;
				}
				_TauSAcceptTarget = 0.44;
				_TauSUpdateFreq = 25;
				_TauSAnyUpdates = true;


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

				_uCARStdDev=1.0;
				_uCARStdDevLower=0.1;
				_uCARStdDevUpper=9.9;
				_nTryuCAR=0;
				_nAcceptuCAR=0;
				_nLocalAcceptuCAR=0;
				_nResetuCAR=0;
				_uCARAcceptTarget=0.44;
				_uCARUpdateFreq=25;
				_uCARAnyUpdates=true;


		};

		~pReMiuMPropParams(){};


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
		vector<unsigned int> nTryTauS() const {
			return _nTryTauS;
		}

		unsigned int nTryTauS(const unsigned int& j) const {
			return _nTryTauS[j];
		}

		vector<unsigned int> nAcceptTauS() const {
			return _nAcceptTauS;
		}

		double TauSAcceptRate(const unsigned int& j) const {
			if (_nTryTauS[j]>0) {
				return (double)_nAcceptTauS[j] / (double)_nTryTauS[j];
			}
			else {
				return 0.0;
			}
		}

		unsigned int TauSUpdateFreq() const {
			return _TauSUpdateFreq;
		}

		vector<unsigned int> nLocalAcceptTauS() const {
			return _nLocalAcceptTauS;
		}

		double TauSLocalAcceptRate(const unsigned int& j) const {
			return (double)_nLocalAcceptTauS[j] / (double)_TauSUpdateFreq;
		}

		double TauSAcceptTarget() const {
			return _TauSAcceptTarget;
		}

		void TauSAddTry(const unsigned int& j) {
			_nTryTauS[j]++;
		}

		void TauSAddAccept(const unsigned int& j) {
			_nAcceptTauS[j]++;
			_nLocalAcceptTauS[j]++;
		}

		void TauSLocalReset(const unsigned int& j) {
			_nLocalAcceptTauS[j] = 0;
		}

		vector<unsigned int> nResetTauS() const {
			return _nResetTauS;
		}

		void TauSStdDevReset(const unsigned int& j) {
			_TauSStdDev[j] = 1.0;
			_nResetTauS[j]++;
			_TauSStdDevLower[j] = pow(10.0, -((double)_nResetTauS[j] + 1.0));
			_TauSStdDevUpper[j] = 100.0 - pow(10.0, -((double)_nResetTauS[j] + 1.0));
		}

		vector<double> TauSStdDev() const {
			return _TauSStdDev;
		}

		vector<double>& TauSStdDev() {
			return _TauSStdDev;
		}

		double& TauSStdDev(const unsigned int& j) {
			return _TauSStdDev[j];
		}

		const double& TauSStdDev(const unsigned int& j) const {
			return _TauSStdDev[j];
		}

		vector<double> TauSStdDevLower() const {
			return _TauSStdDevLower;
		}

		double TauSStdDevLower(const unsigned int& j) const {
			return _TauSStdDevLower[j];
		}

		vector<double> TauSStdDevUpper() const {
			return _TauSStdDevUpper;
		}

		double TauSStdDevUpper(const unsigned int& j) const {
			return _TauSStdDevUpper[j];
		}

		void TauSStdDev(const unsigned int& j, const double& sd) {
			_TauSStdDev[j] = sd;
		}

		bool TauSAnyUpdates() const {
			return _TauSAnyUpdates;
		}

		void TauSAnyUpdates(const bool& newStatus) {
			_TauSAnyUpdates = newStatus;
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

		unsigned int nTryKappa1() const {
			return _nTryKappa1;
		}

		unsigned int nAcceptKappa1() const {
			return _nAcceptKappa1;
		}

		double kappa1AcceptRate() const {
			if (_nTryKappa1>0) {
				return (double)_nAcceptKappa1 / (double)_nTryKappa1;
			}
			else {
				return 0.0;
			}
		}

		unsigned int kappa1UpdateFreq() const {
			return _kappa1UpdateFreq;
		}

		unsigned int nLocalAcceptKappa1() const {
			return _nLocalAcceptKappa1;
		}

		double kappa1LocalAcceptRate() const {
			return (double)_nLocalAcceptKappa1 / (double)_kappa1UpdateFreq;
		}

		double kappa1AcceptTarget() const {
			return _kappa1AcceptTarget;
		}

		void kappa1AddTry() {
			_nTryKappa1++;
		}

		void kappa1AddAccept() {
			_nAcceptKappa1++;
			_nLocalAcceptKappa1++;
		}

		void kappa1LocalReset() {
			_nLocalAcceptKappa1 = 0;
		}

		unsigned int nResetKappa1() const {
			return _nResetKappa1;
		}

		void kappa1StdDevReset() {
			_kappa1StdDev = 1.0;
			_nResetKappa1++;
			_kappa1StdDevLower = pow(10.0, -((double)_nResetKappa1 + 1.0));
			_kappa1StdDevUpper = 100.0 - pow(10.0, -((double)_nResetKappa1 + 1.0));
		}

		const double kappa1StdDev() const {
			return _kappa1StdDev;
		}

		double& kappa1StdDev() {
			return _kappa1StdDev;
		}

		double kappa1StdDevLower() const {
			return _kappa1StdDevLower;
		}

		double kappa1StdDevUpper() const {
			return _kappa1StdDevUpper;
		}

		void kappa1StdDev(const double& sd) {
			_kappa1StdDev = sd;
		}

		bool kappa1AnyUpdates() const {
			return _kappa1AnyUpdates;
		}

		void kappa1AnyUpdates(const bool& newStatus) {
			_kappa1AnyUpdates = newStatus;
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



		unsigned int nTryuCAR() const{
			return _nTryuCAR;
		}

		unsigned int nAcceptuCAR() const{
			return _nAcceptuCAR;
		}

		double uCARAcceptRate() const{
			if(_nTryuCAR>0){
				return (double)_nAcceptuCAR/(double)_nTryuCAR;
			}else{
				return 0.0;
			}
		}

		unsigned int uCARUpdateFreq() const{
			return _uCARUpdateFreq;
		}

		unsigned int nLocalAcceptuCAR() const{
			return _nLocalAcceptuCAR;
		}

		double uCARLocalAcceptRate() const{
			return (double)_nLocalAcceptuCAR/(double)_uCARUpdateFreq;
		}


		double uCARAcceptTarget() const{
			return _uCARAcceptTarget;
		}

		void uCARAddTry(){
			_nTryuCAR++;
		}

		void uCARAddAccept(){
			_nAcceptuCAR++;
			_nLocalAcceptuCAR++;
		}

		void uCARLocalReset(){
			_nLocalAcceptuCAR=0;
		}

		unsigned int nResetuCAR() const{
			return _nResetuCAR;
		}

		void uCARStdDevReset(){
			_uCARStdDev = 1.0;
			_nResetuCAR++;
			_uCARStdDevLower = pow(10.0,-((double)_nResetuCAR+1.0));
			_uCARStdDevUpper = 100.0-pow(10.0,-((double)_nResetuCAR+1.0));
		}

		double& uCARStdDev(){
			return _uCARStdDev;
		}

		double uCARStdDev() const{
			return _uCARStdDev;
		}

		// Member function for setting the standard deviation for
		// proposal for uCAR
		void uCARStdDev(const double& sd){
			_uCARStdDev=sd;
		}

		double uCARStdDevLower() const{
			return _uCARStdDevLower;
		}

		double uCARStdDevUpper() const{
			return _uCARStdDevUpper;
		}

		bool uCARAnyUpdates() const{
			return _uCARAnyUpdates;
		}

		void uCARAnyUpdates(const bool& newStatus){
			_uCARAnyUpdates = newStatus;
		}


		// Need to define a copy iterator
		pReMiuMPropParams& operator=(const pReMiuMPropParams& propParams){
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
			_nTryTauS = propParams.nTryTauS();
			_nAcceptTauS = propParams.nAcceptTauS();
			_nLocalAcceptTauS = propParams.nLocalAcceptTauS();
			_nResetTauS = propParams.nResetTauS();
			_TauSStdDev = propParams.TauSStdDev();
			_TauSStdDevLower = propParams.TauSStdDevLower();
			_TauSStdDevUpper = propParams.TauSStdDevUpper();
			_TauSAcceptTarget = propParams.TauSAcceptTarget();
			_TauSUpdateFreq = propParams.TauSUpdateFreq();
			_TauSAnyUpdates = propParams.TauSAnyUpdates();
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
			_nTryKappa1 = propParams.nTryKappa1();
			_nAcceptKappa1 = propParams.nAcceptKappa1();
			_nLocalAcceptKappa1 = propParams.nLocalAcceptKappa1();
			_nResetKappa1 = propParams.nResetKappa1();
			_kappa1StdDev = propParams.kappa1StdDev();
			_kappa1StdDevLower = propParams.kappa1StdDevLower();
			_kappa1StdDevUpper = propParams.kappa1StdDevUpper();
			_kappa1AcceptTarget = propParams.kappa1AcceptTarget();
			_kappa1UpdateFreq = propParams.kappa1UpdateFreq();
			_kappa1AnyUpdates = propParams.kappa1AnyUpdates(); 
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
			_nTryuCAR=propParams.nTryuCAR();
			_nAcceptuCAR=propParams.nAcceptuCAR();
			_nLocalAcceptuCAR=propParams.nLocalAcceptuCAR();
			_nResetuCAR=propParams.nResetuCAR();
			_uCARStdDev=propParams.uCARStdDev();
			_uCARStdDevLower=propParams.uCARStdDevLower();
			_uCARStdDevUpper=propParams.uCARStdDevUpper();
			_uCARAcceptTarget=propParams.uCARAcceptTarget();
			_uCARUpdateFreq=propParams.uCARUpdateFreq();
			_uCARAnyUpdates=propParams.uCARAnyUpdates();

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
		vector<unsigned int> _nTryTauS;
		vector<unsigned int> _nAcceptTauS;
		vector<unsigned int> _nLocalAcceptTauS;
		vector<unsigned int> _nResetTauS;
		vector<double> _TauSStdDev;
		vector<double> _TauSStdDevLower;
		vector<double> _TauSStdDevUpper;
		double _TauSAcceptTarget;
		unsigned int _TauSUpdateFreq;
		bool _TauSAnyUpdates;
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
		unsigned int _nTryKappa1;
		unsigned int _nAcceptKappa1;
		unsigned int _nLocalAcceptKappa1;
		unsigned int _nResetKappa1;
		double _kappa1StdDev;
		double _kappa1StdDevLower;
		double _kappa1StdDevUpper;
		double _kappa1AcceptTarget;
		unsigned int _kappa1UpdateFreq;
		bool _kappa1AnyUpdates; 
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

		unsigned int _nTryuCAR;
		unsigned int _nAcceptuCAR;
		unsigned int _nLocalAcceptuCAR;
		unsigned int _nResetuCAR;
		double _uCARStdDev;
		double _uCARStdDevLower;
		double _uCARStdDevUpper;
		double _uCARAcceptTarget;
		unsigned int _uCARUpdateFreq;
		bool _uCARAnyUpdates;


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
void gibbsForVActive(mcmcChain<pReMiuMParams>& chain,
				unsigned int& nTry,unsigned int& nAccept,
				const mcmcModel<pReMiuMParams,pReMiuMOptions,pReMiuMData>& model,
				pReMiuMPropParams& propParams,
				baseGeneratorType& rndGenerator){

	mcmcState<pReMiuMParams>& currentState = chain.currentState();
	pReMiuMParams& currentParams = currentState.parameters();

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
	double dPitmanYor = currentParams.dPitmanYor();
	for(unsigned int c=0;c<=maxZ;c++){
		double vVal = betaRand(rndGenerator,1.0+currentParams.workNXInCluster(c)-dPitmanYor,alpha+sumCPlus1ToMaxMembers[c]+dPitmanYor*(c+1));
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
void updateForPhiActive(mcmcChain<pReMiuMParams>& chain,
				unsigned int& nTry,unsigned int& nAccept,
				const mcmcModel<pReMiuMParams,pReMiuMOptions,pReMiuMData>& model,
				pReMiuMPropParams& propParams,
				baseGeneratorType& rndGenerator){

	mcmcState<pReMiuMParams>& currentState = chain.currentState();
	pReMiuMParams& currentParams = currentState.parameters();
	pReMiuMHyperParams hyperParams = currentParams.hyperParams();


	const pReMiuMData& dataset = model.dataset();
	string varSelectType = model.options().varSelectType();
	// Find the number of clusters
	unsigned int maxZ = currentParams.workMaxZi();
	// Find the number of covariates
	unsigned int nCovariates = 0;
	if(model.options().covariateType().compare("Mixed")==0){
		nCovariates = currentParams.nDiscreteCovs();
	} else {
		nCovariates = currentParams.nCovariates();
	}
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
void gibbsForMuActive(mcmcChain<pReMiuMParams>& chain,
				unsigned int& nTry,unsigned int& nAccept,
				const mcmcModel<pReMiuMParams,pReMiuMOptions,pReMiuMData>& model,
				pReMiuMPropParams& propParams,
				baseGeneratorType& rndGenerator){

	mcmcState<pReMiuMParams>& currentState = chain.currentState();
	pReMiuMParams& currentParams = currentState.parameters();
	pReMiuMHyperParams hyperParams = currentParams.hyperParams();

	const pReMiuMData& dataset = model.dataset();
	bool useIndependentNormal = model.options().useIndependentNormal();
	bool useHyperpriorR1 = model.options().useHyperpriorR1();
	bool useSeparationPrior = model.options().useSeparationPrior();

	// Find the number of clusters
	unsigned int maxZ = currentParams.workMaxZi();
	// Find the number of covariates
	unsigned int nCovariates = 0;
	if(model.options().covariateType().compare("Mixed")==0){
		nCovariates = currentParams.nContinuousCovs();
	} else {
		nCovariates = currentParams.nCovariates();
	}
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
			gammaMat[c](j,j)=currentParams.gamma(c,currentParams.nDiscreteCovs()+j);
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
		VectorXd meanVec(nCovariates);

		if (useHyperpriorR1) {
			covMat = (currentParams.Tau00() + nXInC*gammaMat[c] * currentParams.Tau(c)*gammaMat[c]).inverse();	
			meanVec = currentParams.Tau00()*currentParams.mu00() +
				nXInC*gammaMat[c] * currentParams.Tau(c)*(meanX[c] - oneMinusGammaMat[c] * currentParams.nullMu());
			meanVec = covMat*meanVec;
		}
		else if (useSeparationPrior) {
			covMat = (hyperParams.Tau00() + nXInC*gammaMat[c] * currentParams.Tau(c)*gammaMat[c]).inverse();
			meanVec = hyperParams.Tau00()*currentParams.mu00() +
				nXInC*gammaMat[c] * currentParams.Tau(c)*(meanX[c] - oneMinusGammaMat[c] * currentParams.nullMu());
			meanVec = covMat*meanVec;
		} 
		else if (!useIndependentNormal){
			covMat = (hyperParams.Tau0() + nXInC*gammaMat[c] * currentParams.Tau(c)*gammaMat[c]).inverse();
			meanVec = hyperParams.Tau0()*hyperParams.mu0() +
				nXInC*gammaMat[c] * currentParams.Tau(c)*(meanX[c] - oneMinusGammaMat[c] * currentParams.nullMu());
			meanVec = covMat*meanVec;
		}

		VectorXd mu(nCovariates);
		// We sample from this posterior
		mu = multivarNormalRand(rndGenerator, meanVec, covMat);

		// We store our sample
		currentParams.mu(c,mu, useIndependentNormal);
	}
}

// Gibbs update for mu in Normal covariate case and use of normal inverse prior
void gibbsForMuActiveNIWP(mcmcChain<pReMiuMParams>& chain,
				unsigned int& nTry,unsigned int& nAccept,
				const mcmcModel<pReMiuMParams,pReMiuMOptions,pReMiuMData>& model,
				pReMiuMPropParams& propParams,
				baseGeneratorType& rndGenerator){

	mcmcState<pReMiuMParams>& currentState = chain.currentState();
	pReMiuMParams& currentParams = currentState.parameters();
	pReMiuMHyperParams hyperParams = currentParams.hyperParams();

	const pReMiuMData& dataset = model.dataset();
	bool useIndependentNormal = model.options().useIndependentNormal();

	// Find the number of clusters
	unsigned int maxZ = currentParams.workMaxZi();
	// Find the number of covariates
	unsigned int nCovariates = 0;
	if(model.options().covariateType().compare("Mixed")==0){
		nCovariates = currentParams.nContinuousCovs();
	} else {
		nCovariates = currentParams.nCovariates();
	}
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
			gammaMat[c](j,j)=currentParams.gamma(c,currentParams.nDiscreteCovs()+j);
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
		covMat=(gammaMat[c]*currentParams.Sigma(c)*gammaMat[c])/(hyperParams.nu0()+nXInC);
		//There are 0 on the diagonal elements of the covariance matrix when the covariate j is not selected
		//we replace them by 0.1, the value does not care, since the sampling values will be replaced by \bar{x}_j
		for (unsigned int j=0; j<nCovariates; j++){
           		if (covMat(j,j)==0) covMat(j,j)=0.1;
		}
		VectorXd meanVec(nCovariates);
        	meanVec = hyperParams.nu0()*hyperParams.mu0()+
					nXInC*(gammaMat[c]*meanX[c]-oneMinusGammaMat[c]*currentParams.nullMu());
		meanVec /= (hyperParams.nu0()+nXInC);
		VectorXd mu(nCovariates);
		// We sample from this posterior
		mu = multivarNormalRand(rndGenerator,meanVec,covMat);

		// We store our sample
		currentParams.mu(c, mu, useIndependentNormal);

	}

}



// Gibbs update for mu in independent Normal covariate case 
void gibbsForMuActiveIndep(mcmcChain<pReMiuMParams>& chain,
	unsigned int& nTry, unsigned int& nAccept,
	const mcmcModel<pReMiuMParams, pReMiuMOptions, pReMiuMData>& model,
	pReMiuMPropParams& propParams,
	baseGeneratorType& rndGenerator) {

	mcmcState<pReMiuMParams>& currentState = chain.currentState();
	pReMiuMParams& currentParams = currentState.parameters();
	pReMiuMHyperParams hyperParams = currentParams.hyperParams();

	const pReMiuMData& dataset = model.dataset();
	bool useIndependentNormal = model.options().useIndependentNormal(); 

	// Find the number of clusters
	unsigned int maxZ = currentParams.workMaxZi();
	// Find the number of covariates
	unsigned int nCovariates = 0;
	if (model.options().covariateType().compare("Mixed") == 0) {
		nCovariates = currentParams.nContinuousCovs();
	}
	else {
		nCovariates = currentParams.nCovariates();
	}
	// Find the number of subjects
	unsigned int nSubjects = dataset.nSubjects();


	nTry++;
	nAccept++;

	// In the following it is useful to have the rows of X as
	// Eigen dynamic vectors
	vector<VectorXd> xi(nSubjects);
	for (unsigned int i = 0; i<nSubjects; i++) {
		xi[i].setZero(nCovariates);
		for (unsigned int j = 0; j<nCovariates; j++) {
			xi[i](j) = dataset.continuousX(i, j);
		}
	}

	// We begin by computing the mean X for individuals in each cluster
	vector<VectorXd> meanX(maxZ + 1);
	for (unsigned int c = 0; c <= maxZ; c++) {
		meanX[c].setZero(nCovariates);
	}

	for (unsigned int i = 0; i<nSubjects; i++) {
		meanX[currentParams.z(i)] = meanX[currentParams.z(i)] + xi[i];
	}
	for (unsigned int c = 0; c <= maxZ; c++) {
		int nXInC = currentParams.workNXInCluster(c);
		if (nXInC>0) {
			meanX[c] = meanX[c] / (double)nXInC;
		}
		else {
			meanX[c].setZero(nCovariates);
		}
	}

	// Having computed this we can calcuate the posterior mean
	// and posterior covariance for each mu_c,j

	//initialize gamma_cj and 1-gamma_cj used for variable selection 
	double gamma_cj = 0.0;
	double oneMinusGamma_cj = 0.0;
	VectorXd mu0 = hyperParams.mu0();
	VectorXd Tau0 = hyperParams.Tau0_Indep();
	VectorXd nullMu = currentParams.nullMu();
	
	for (unsigned int c = 0; c <= maxZ; c++) {

		int nXInC = currentParams.workNXInCluster(c);
		VectorXd mu(nCovariates);

		// Loop over the covariates
		for (unsigned int j = 0; j<nCovariates; j++) {

			gamma_cj=currentParams.gamma(c, currentParams.nDiscreteCovs() + j);
			oneMinusGamma_cj = 1 - gamma_cj;

			double denom = 0.0;
			denom = nXInC*(1.0 / Tau0(j))*gamma_cj*gamma_cj + (1.0 / currentParams.Tau_Indep(c, j));
			double meanNum = 0.0;
			meanNum = nXInC*(1.0 / Tau0(j))*meanX[c](j)*gamma_cj + (1.0 / currentParams.Tau_Indep(c, j))*mu0(j)
				- nXInC*(1.0 / Tau0(j))*gamma_cj*oneMinusGamma_cj*nullMu(j);
			double variance = (1.0 / currentParams.Tau_Indep(c, j))*(1.0 / Tau0(j)) / denom;
			double mean = meanNum / denom;
				    
			// We sample from this posterior
			mu (j)= NormalRand(rndGenerator, mean, variance);

		}
		// We store our sample
		currentParams.mu(c, mu, useIndependentNormal);
		
	}
}



// Gibbs update for Tau in the Normal covariate case
void gibbsForTauActive(mcmcChain<pReMiuMParams>& chain,
				unsigned int& nTry,unsigned int& nAccept,
				const mcmcModel<pReMiuMParams,pReMiuMOptions,pReMiuMData>& model,
				pReMiuMPropParams& propParams,
				baseGeneratorType& rndGenerator){

	mcmcState<pReMiuMParams>& currentState = chain.currentState();
	pReMiuMParams& currentParams = currentState.parameters();
	pReMiuMHyperParams hyperParams = currentParams.hyperParams();
	bool useHyperpriorR1 = model.options().useHyperpriorR1();

	const pReMiuMData& dataset = model.dataset();

	// Find the number of clusters
	unsigned int maxZ = currentParams.workMaxZi();
	// Find the number of covariates
	unsigned int nCovariates = 0;
	if(model.options().covariateType().compare("Mixed")==0){
		nCovariates = currentParams.nContinuousCovs();
	} else {
		nCovariates = currentParams.nCovariates();
	}
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

	if (useHyperpriorR1){
		for(unsigned int c=0;c<=maxZ;c++){
			Rc[c]=(currentParams.R1().inverse()+Rc[c]).inverse();
			MatrixXd Tau = wishartRand(rndGenerator,Rc[c],currentParams.workNXInCluster(c)+currentParams.kappa11());
	
			currentParams.Tau(c,Tau);
		}
	} else {
		for(unsigned int c=0;c<=maxZ;c++){
			Rc[c]=(hyperParams.R0().inverse()+Rc[c]).inverse();
			MatrixXd Tau = wishartRand(rndGenerator,Rc[c],currentParams.workNXInCluster(c)+hyperParams.kappa0());
	
			currentParams.Tau(c,Tau);
		}
	}


}

// Gibbs update for TauR in the Normal covariate case when separation strategy is applied
void gibbsForTauRActive(mcmcChain<pReMiuMParams>& chain,
	unsigned int& nTry, unsigned int& nAccept,
	const mcmcModel<pReMiuMParams, pReMiuMOptions, pReMiuMData>& model,
	pReMiuMPropParams& propParams,
	baseGeneratorType& rndGenerator) {

	mcmcState<pReMiuMParams>& currentState = chain.currentState();
	pReMiuMParams& currentParams = currentState.parameters();
	pReMiuMHyperParams hyperParams = currentParams.hyperParams();

	const pReMiuMData& dataset = model.dataset();

	// Find the number of clusters
	unsigned int maxZ = currentParams.workMaxZi();
	// Find the number of covariates
	unsigned int nCovariates = 0;
	if (model.options().covariateType().compare("Mixed") == 0) {
		nCovariates = currentParams.nContinuousCovs();
	}
	else {
		nCovariates = currentParams.nCovariates();
	}
	// Find the number of subjects
	unsigned int nSubjects = dataset.nSubjects();

	nTry++;
	nAccept++;

	// In the following it is useful to have the rows of X as
	// Eigen dynamic vectors
	vector<VectorXd> xi(nSubjects);
	for (unsigned int i = 0; i<nSubjects; i++) {
		xi[i].setZero(nCovariates);
		for (unsigned int j = 0; j<nCovariates; j++) {
			xi[i](j) = dataset.continuousX(i, j);
		}
	}

	vector<MatrixXd> Rc(maxZ + 1);
	for (unsigned int c = 0; c <= maxZ; c++) {
		Rc[c].setZero(nCovariates, nCovariates);
	}

	for (unsigned int i = 0; i<nSubjects; i++) {
		unsigned int zi = currentParams.z(i);
		Rc[zi] = Rc[zi] + (xi[i] - currentParams.workMuStar(zi))*((xi[i] - currentParams.workMuStar(zi)).transpose());
	}

	for (unsigned int c = 0; c <= maxZ; c++) {
		Rc[c] = (hyperParams.R0().inverse() + currentParams.TauS(c)*Rc[c]*currentParams.TauS(c)).inverse();
		MatrixXd TauR = wishartRand(rndGenerator, Rc[c], currentParams.workNXInCluster(c) + currentParams.kappa11());

		currentParams.TauR(c, TauR);
	}

}

// Adaptive Metropolis Hastings move for taus
void metropolisHastingsForTauS(mcmcChain<pReMiuMParams>& chain,
	unsigned int& nTry, unsigned int& nAccept,
	const mcmcModel<pReMiuMParams, pReMiuMOptions, pReMiuMData>& model,
	pReMiuMPropParams& propParams,
	baseGeneratorType& rndGenerator) {

	mcmcState<pReMiuMParams>& currentState = chain.currentState();
	pReMiuMParams& currentParams = currentState.parameters();
	pReMiuMHyperParams hyperParams = currentParams.hyperParams();

	// Find the number of clusters
	unsigned int maxZ = currentParams.workMaxZi();
	// Find the number of covariates
	unsigned int nCovariates = 0;
	if (model.options().covariateType().compare("Mixed") == 0) {
		nCovariates = currentParams.nContinuousCovs();
	}
	else {
		nCovariates = currentParams.nCovariates();
	}

	randomUniform unifRand(0, 1);

	double TauSTargetRate = propParams.TauSAcceptTarget();
	unsigned int TauSUpdateFreq = propParams.TauSUpdateFreq();

	for (unsigned int j = 0; j < nCovariates; j++) {
		for (unsigned int c = 0; c <= maxZ; c++) {
			nTry++;
			propParams.TauSAddTry(j);

			double& stdDev = propParams.TauSStdDev(j);
			double TauSOrig = currentParams.TauS(c,j);
			double currentCondLogPost = logCondPostTauS(currentParams, model, c, j);

			double TauSProp;
			TauSProp = truncNormalRand(rndGenerator, TauSOrig, stdDev, "L", 0, 0);
			currentParams.TauS(c, j, TauSProp);
			double propCondLogPost = logCondPostTauS(currentParams, model, c, j);
			double logAcceptRatio = propCondLogPost - currentCondLogPost;
			
			// Add the proposal contribution
			logAcceptRatio += logPdfTruncatedNormal(TauSOrig, TauSProp, stdDev, "L", 0, 0);
			logAcceptRatio -= logPdfTruncatedNormal(TauSProp, TauSOrig, stdDev, "L", 0, 0);

			if (unifRand(rndGenerator)<exp(logAcceptRatio)) {
				nAccept++;
				propParams.TauSAddAccept(j);
				// Update the std dev of the proposal
				if (propParams.nTryTauS(j) % TauSUpdateFreq == 0) {
					stdDev += 10 * (propParams.TauSLocalAcceptRate(j) - TauSTargetRate) /
						pow((double)(propParams.nTryTauS(j) / TauSUpdateFreq) + 2.0, 0.75);
					propParams.TauSAnyUpdates(true);
					if (stdDev>propParams.TauSStdDevUpper(j) || stdDev<propParams.TauSStdDevLower(j)) {
						propParams.TauSStdDevReset(j);
					}
					propParams.TauSLocalReset(j);
				}
			}
			else {
				currentParams.TauS(c, j, TauSOrig);
				// Update the std dev of the proposal
				if (propParams.nTryTauS(j) % TauSUpdateFreq == 0) {
					stdDev += 10 * (propParams.TauSLocalAcceptRate(j) - TauSTargetRate) /
						pow((double)(propParams.nTryTauS(j) / TauSUpdateFreq) + 2.0, 0.75);
					propParams.TauSAnyUpdates(true);
					if (stdDev<propParams.TauSStdDevLower(j) || stdDev>propParams.TauSStdDevUpper(j)) {
						propParams.TauSStdDevReset(j);
					}
					propParams.TauSLocalReset(j);
				}
			}
		 }	
	}
}



// Gibbs update for Tau in the independent Normal case
void gibbsForTauActiveIndep(mcmcChain<pReMiuMParams>& chain,
	unsigned int& nTry, unsigned int& nAccept,
	const mcmcModel<pReMiuMParams, pReMiuMOptions, pReMiuMData>& model,
	pReMiuMPropParams& propParams,
	baseGeneratorType& rndGenerator) {

	mcmcState<pReMiuMParams>& currentState = chain.currentState();
	pReMiuMParams& currentParams = currentState.parameters();
	pReMiuMHyperParams hyperParams = currentParams.hyperParams();

	const pReMiuMData& dataset = model.dataset();

	// Find the number of clusters
	unsigned int maxZ = currentParams.workMaxZi();
	// Find the number of covariates
	unsigned int nCovariates = 0;
	if (model.options().covariateType().compare("Mixed") == 0) {
		nCovariates = currentParams.nContinuousCovs();
	}
	else {
		nCovariates = currentParams.nCovariates();
	}
	// Find the number of subjects
	unsigned int nSubjects = dataset.nSubjects();

	nTry++;
	nAccept++;

	// In the following it is useful to have the rows of X as
	// Eigen dynamic vectors
	vector<VectorXd> xi(nSubjects);
	for (unsigned int i = 0; i<nSubjects; i++) {
		xi[i].setZero(nCovariates);
		for (unsigned int j = 0; j<nCovariates; j++) {
			xi[i](j) = dataset.continuousX(i, j);
		}
	}

	vector<VectorXd> XiMinusMuStarSq(nSubjects);
	for (unsigned int i = 0; i < nSubjects; i++) {
		XiMinusMuStarSq[i].setZero(nCovariates);
		VectorXd muStar = currentParams.workMuStar(currentParams.z(i));
		for (unsigned int j = 0; j < nCovariates; j++) {
			XiMinusMuStarSq[i](j) = (xi[i](j) - muStar(j))*(xi[i](j) - muStar(j));
		}
	}

	// We begin by computing the sum of X minus mu_star in each cluster
	vector<VectorXd> sumXiMinusMuStarSq(maxZ + 1);
	for (unsigned int c = 0; c <= maxZ; c++) {
		sumXiMinusMuStarSq[c].setZero(nCovariates);
	}

	for (unsigned int i = 0; i<nSubjects; i++) {
		sumXiMinusMuStarSq[currentParams.z(i)] = sumXiMinusMuStarSq[currentParams.z(i)] + XiMinusMuStarSq[i];
	}


	for (unsigned int c = 0; c <= maxZ; c++) {
		VectorXd tau(nCovariates);
		int nXInC = currentParams.workNXInCluster(c);
		for (unsigned int j = 0; j < nCovariates; j++) {
			double kappaNew = (double) nXInC / 2 + hyperParams.kappa1();
			double rNew = (sumXiMinusMuStarSq[c](j) + 2 * currentParams.R1_Indep(j))/2;

			randomGamma gammaRand(kappaNew, 1.0 / rNew);
			tau(j) = gammaRand(rndGenerator);
		}
		currentParams.Tau_Indep(c, tau);
	}

}

// Gibbs update for update of gamma (only used in the binary variable selection case)
void gibbsForGammaActive(mcmcChain<pReMiuMParams>& chain,
					unsigned int& nTry,unsigned int& nAccept,
					const mcmcModel<pReMiuMParams,
									pReMiuMOptions,
									pReMiuMData>& model,
					pReMiuMPropParams& propParams,
					baseGeneratorType& rndGenerator){


	mcmcState<pReMiuMParams>& currentState = chain.currentState();
	pReMiuMParams& currentParams = currentState.parameters();
	pReMiuMHyperParams hyperParams = currentParams.hyperParams();
	

	// Find the number of subjects
	unsigned int nCovariates = currentParams.nCovariates();
	unsigned int nSubjects = currentParams.nSubjects();
	unsigned int maxZ = currentParams.workMaxZi();
	string covariateType = model.options().covariateType();
	string varSelectType = model.options().varSelectType();
	bool useIndependentNormal = model.options().useIndependentNormal();

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
				currentParams.gamma(c,j,currentGamma[j],covariateType,useIndependentNormal);

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
				currentParams.gamma(c,j,currentGamma[j],covariateType, useIndependentNormal);
			}
			// Otherwise switching but nothing to do here as we had already done the
			// switch in the calculations

		}
	}

}


// Adaptive Metropolis-Hastings for theta
void metropolisHastingsForThetaActive(mcmcChain<pReMiuMParams>& chain,
								unsigned int& nTry,unsigned int& nAccept,
								const mcmcModel<pReMiuMParams,
												pReMiuMOptions,
												pReMiuMData>& model,
								pReMiuMPropParams& propParams,
								baseGeneratorType& rndGenerator){

	mcmcState<pReMiuMParams>& currentState = chain.currentState();
	pReMiuMParams& currentParams = currentState.parameters();
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
void metropolisHastingsForLabels123(mcmcChain<pReMiuMParams>& chain,
				unsigned int& nTry,unsigned int& nAccept,
				const mcmcModel<pReMiuMParams,pReMiuMOptions,pReMiuMData>& model,
				pReMiuMPropParams& propParams,
				baseGeneratorType& rndGenerator){

	mcmcState<pReMiuMParams>& currentState = chain.currentState();
	pReMiuMParams& currentParams = currentState.parameters();

	unsigned int maxZ = currentParams.workMaxZi();
	if(maxZ==0){
		// If there is only one cluster with individuals in, don't do anything
		return;
	}
	string varSelectType = model.options().varSelectType();
	string covariateType = model.options().covariateType();
	bool useIndependentNormal = model.options().useIndependentNormal();
	bool useSeparationPrior = model.options().useSeparationPrior();

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
	nTry++;
	unsigned int i1=(unsigned int)nNotEmpty*unifRand(rndGenerator);
	unsigned int c1=nonEmptyIndices[i1];
	nonEmptyIndices.erase(nonEmptyIndices.begin()+i1);
	unsigned int i2=(unsigned int)(nNotEmpty-1)*unifRand(rndGenerator);
	unsigned int c2=nonEmptyIndices[i2];

	// Check whether we accept the move
	double logAcceptRatio = ((double)currentParams.workNXInCluster(c2)-
		(double)currentParams.workNXInCluster(c1))
		*(currentParams.logPsi(c1)-currentParams.logPsi(c2));

	if(unifRand(rndGenerator)<exp(logAcceptRatio)){
//		nAccept++;
		// Switch the labels
		currentParams.switchLabels(c1, c2, covariateType, varSelectType, useIndependentNormal, useSeparationPrior);
	}

	// Move 2 - swap labels of 2 randomly selected neighbouring clusters,
	//			also swapping the v at the same time
//	nTry++;
	c1=(unsigned int)maxZ*unifRand(rndGenerator);

	logAcceptRatio=(double)currentParams.workNXInCluster(c1)*
		log(1-currentParams.v(c1+1))
		- (double)currentParams.workNXInCluster(c1+1)*log(1-currentParams.v(c1));

	if(unifRand(rndGenerator)<exp(logAcceptRatio)){
		nAccept++;

		// Switch the labels
		currentParams.switchLabels(c1, c1 + 1, covariateType, varSelectType, useIndependentNormal, useSeparationPrior);

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
				maxZ=c1;
			}
		}


	}

	// Move 3

//	nTry++;
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
//		nAccept++;
		currentParams.switchLabels(c1, c1 + 1, covariateType, varSelectType, useIndependentNormal, useSeparationPrior);
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

void metropolisHastingsForLabels12(mcmcChain<pReMiuMParams>& chain,
				unsigned int& nTry,unsigned int& nAccept,
				const mcmcModel<pReMiuMParams,pReMiuMOptions,pReMiuMData>& model,
				pReMiuMPropParams& propParams,
				baseGeneratorType& rndGenerator){

	mcmcState<pReMiuMParams>& currentState = chain.currentState();
	pReMiuMParams& currentParams = currentState.parameters();

	unsigned int maxZ = currentParams.workMaxZi();
	if(maxZ==0){
		// If there is only one cluster with individuals in, don't do anything
		return;
	}
	string varSelectType = model.options().varSelectType();
	string covariateType = model.options().covariateType();
	bool useIndependentNormal = model.options().useIndependentNormal();
	bool useSeparationPrior = model.options().useSeparationPrior();

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
	nTry++;
	unsigned int i1=(unsigned int)nNotEmpty*unifRand(rndGenerator);
	unsigned int c1=nonEmptyIndices[i1];
	nonEmptyIndices.erase(nonEmptyIndices.begin()+i1);
	unsigned int i2=(unsigned int)(nNotEmpty-1)*unifRand(rndGenerator);
	unsigned int c2=nonEmptyIndices[i2];

	// Check whether we accept the move
	double logAcceptRatio = ((double)currentParams.workNXInCluster(c2)-(double)currentParams.workNXInCluster(c1))
								*(currentParams.logPsi(c1)-currentParams.logPsi(c2));

	if(unifRand(rndGenerator)<exp(logAcceptRatio)){
		nAccept++;
		// Switch the labels
		currentParams.switchLabels(c1, c2, covariateType, varSelectType, useIndependentNormal, useSeparationPrior);
	}

	// Move 2 - swap labels of 2 randomly selected neighbouring clusters,
	//			also swapping the v at the same time
//	nTry++;
	c1=(unsigned int)maxZ*unifRand(rndGenerator);

	logAcceptRatio=(double)currentParams.workNXInCluster(c1)*log(1-currentParams.v(c1+1))
							- (double)currentParams.workNXInCluster(c1+1)*log(1-currentParams.v(c1));

	if(unifRand(rndGenerator)<exp(logAcceptRatio)){
//		nAccept++;
		// Switch the labels
		currentParams.switchLabels(c1, c1 + 1, covariateType, varSelectType, useIndependentNormal, useSeparationPrior);

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

}

void metropolisHastingsForLabels3(mcmcChain<pReMiuMParams>& chain,
				unsigned int& nTry,unsigned int& nAccept,
				const mcmcModel<pReMiuMParams,pReMiuMOptions,pReMiuMData>& model,
				pReMiuMPropParams& propParams,
				baseGeneratorType& rndGenerator){

	mcmcState<pReMiuMParams>& currentState = chain.currentState();
	pReMiuMParams& currentParams = currentState.parameters();

	unsigned int maxZ = currentParams.workMaxZi();
	if(maxZ==0){
		// If there is only one cluster with individuals in, don't do anything
		return;
	}
	string varSelectType = model.options().varSelectType();
	string covariateType = model.options().covariateType();
	bool useIndependentNormal = model.options().useIndependentNormal();
	bool useSeparationPrior= model.options().useSeparationPrior();

	randomUniform unifRand(0,1);

	// Move 3

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
	nTry++;
	unsigned int i1=(unsigned int)nNotEmpty*unifRand(rndGenerator);
	unsigned int c1=nonEmptyIndices[i1];
	nonEmptyIndices.erase(nonEmptyIndices.begin()+i1);

	// Check whether we accept the move
	double logAcceptRatio=0;

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
					log(exp(currentParams.logPsi(c1+1))*const1+exp(currentParams.logPsi(c1))*const2);
	logAcceptRatio+=(double)(currentParams.workNXInCluster(c1+1))*log(const1);
	logAcceptRatio+=(double)(currentParams.workNXInCluster(c1))*log(const2);

	if(unifRand(rndGenerator)<exp(logAcceptRatio)){
		nAccept++;
		currentParams.switchLabels(c1, c1 + 1, covariateType, varSelectType, useIndependentNormal, useSeparationPrior);
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

void gibbsForU(mcmcChain<pReMiuMParams>& chain,
				unsigned int& nTry,unsigned int& nAccept,
				const mcmcModel<pReMiuMParams,pReMiuMOptions,pReMiuMData>& model,
				pReMiuMPropParams& propParams,
				baseGeneratorType& rndGenerator){

	mcmcState<pReMiuMParams>& currentState = chain.currentState();
	pReMiuMParams& currentParams = currentState.parameters();
	pReMiuMHyperParams hyperParams = currentParams.hyperParams();
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
void metropolisHastingsForAlpha(mcmcChain<pReMiuMParams>& chain,
								unsigned int& nTry,unsigned int& nAccept,
								const mcmcModel<pReMiuMParams,
												pReMiuMOptions,
												pReMiuMData>& model,
								pReMiuMPropParams& propParams,
								baseGeneratorType& rndGenerator){

	mcmcState<pReMiuMParams>& currentState = chain.currentState();
	pReMiuMParams& currentParams = currentState.parameters();
	pReMiuMHyperParams hyperParams = currentParams.hyperParams();

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
void gibbsForVInActive(mcmcChain<pReMiuMParams>& chain,
				unsigned int& nTry,unsigned int& nAccept,
				const mcmcModel<pReMiuMParams,pReMiuMOptions,pReMiuMData>& model,
				pReMiuMPropParams& propParams,
				baseGeneratorType& rndGenerator){

	mcmcState<pReMiuMParams>& currentState = chain.currentState();
	pReMiuMParams& currentParams = currentState.parameters();
	pReMiuMHyperParams hyperParams = currentParams.hyperParams();
	string samplerType = model.options().samplerType();
	string covariateType = model.options().covariateType();
	bool useIndependentNormal = model.options().useIndependentNormal();
	bool useSeparationPrior= model.options().useSeparationPrior();

	nTry++;
	nAccept++;

	unsigned int maxZ = currentParams.workMaxZi();
	unsigned int maxNClusters = currentParams.maxNClusters();

	double minUi = currentParams.workMinUi();

	vector<double> vNew=currentParams.v();
	vector<double> logPsiNew=currentParams.logPsi();

	double alpha = currentParams.alpha();
	double dPitmanYor = currentParams.dPitmanYor();

	if(samplerType.compare("Truncated")==0){
		// Just sample from the prior

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
		currentParams.maxNClusters(maxNClusters, covariateType, useIndependentNormal, useSeparationPrior);

	}

	currentParams.v(vNew);
	currentParams.logPsi(logPsiNew);

}

/*********** BLOCK 3 p(Theta^I|.) **********************************/
// I=Inactive. Sample the inactive cluster variables from the prior
// Theta contains phi, mu, Tau, gamma, theta. Only need to sample
// up to maxNClusters = max_i{Ci}. Several different routines here for
// each of the variables

// Adaptive Metropolis Hastings move for kappa1
void metropolisHastingsForKappa1(mcmcChain<pReMiuMParams>& chain,
	unsigned int& nTry, unsigned int& nAccept,
	const mcmcModel<pReMiuMParams,
	pReMiuMOptions,
	pReMiuMData>& model,
	pReMiuMPropParams& propParams,
	baseGeneratorType& rndGenerator) {

	mcmcState<pReMiuMParams>& currentState = chain.currentState();
	pReMiuMParams& currentParams = currentState.parameters();
	pReMiuMHyperParams hyperParams = currentParams.hyperParams();

	// Find the number of clusters
	unsigned int maxZ = currentParams.workMaxZi();
	// Find the number of covariates
	unsigned int nCovariates = 0;
	if (model.options().covariateType().compare("Mixed") == 0) {
		nCovariates = currentParams.nContinuousCovs();
	}
	else {
		nCovariates = currentParams.nCovariates();
	}

	// Define a uniform random number generator
	randomUniform unifRand(0, 1);

	double& stdDev = propParams.kappa1StdDev();
	double kappa1Current = currentParams.kappa11();

	double kappa1Prop;
	kappa1Prop = truncNormalRand(rndGenerator, kappa1Current, stdDev, "L", nCovariates, 0);

	double logAcceptRatio = 0.0;
	double workLogDetR1 = currentParams.workLogDetR1();
	for (unsigned int c = 0; c <= maxZ; c++) {
		double workLogDetTau = currentParams.workLogDetTau(c);
		logAcceptRatio += -logMultivarGammaFn(kappa1Prop / 2.0, nCovariates) +
			(-(kappa1Prop*nCovariates) / 2.0)*log(2.0) + (-kappa1Prop / 2.0)*workLogDetR1 +
			((kappa1Prop - nCovariates - 1.0) / 2.0)*workLogDetTau;
		logAcceptRatio -= -logMultivarGammaFn(kappa1Current / 2.0, nCovariates) +
			(-(kappa1Current*nCovariates) / 2.0)*log(2.0) + (-kappa1Current / 2.0)*workLogDetR1 +
			((kappa1Current - nCovariates - 1.0) / 2.0)*workLogDetTau;
	}

	logAcceptRatio += logPdfInverseGamma(kappa1Prop - nCovariates , hyperParams.shapeKappa1(), hyperParams.scaleKappa1());
	logAcceptRatio -= logPdfInverseGamma(kappa1Current - nCovariates , hyperParams.shapeKappa1(), hyperParams.scaleKappa1());

	// Add the proposal contribution
	logAcceptRatio += logPdfTruncatedNormal(kappa1Current, kappa1Prop, stdDev, "L", nCovariates , 0);
	logAcceptRatio -= logPdfTruncatedNormal(kappa1Prop, kappa1Current, stdDev, "L", nCovariates , 0);

	propParams.kappa1AddTry();
	nTry++;
	if (unifRand(rndGenerator)<exp(logAcceptRatio)) {
		nAccept++;
		propParams.kappa1AddAccept();
		// If the move was accepted update the state
		currentParams.kappa11(kappa1Prop);

		// Also update the proposal standard deviation
		if (propParams.nTryKappa1() % propParams.kappa1UpdateFreq() == 0) {
			stdDev += 10 * (propParams.kappa1LocalAcceptRate() - propParams.kappa1AcceptTarget()) /
				pow((double)(propParams.nTryKappa1() / propParams.kappa1UpdateFreq()) + 2.0, 0.75);
			propParams.kappa1AnyUpdates(true);
			if (stdDev>propParams.kappa1StdDevUpper() || stdDev<propParams.kappa1StdDevLower()) {
				propParams.kappa1StdDevReset();
			}
			propParams.kappa1LocalReset();
		}
	}
	else {
		// Otherwise update the proposal standard deviation
		if (propParams.nTryKappa1() % propParams.kappa1UpdateFreq() == 0) {
			stdDev += 10 * (propParams.kappa1LocalAcceptRate() - propParams.kappa1AcceptTarget()) /
				pow((double)(propParams.nTryKappa1() / propParams.kappa1UpdateFreq()) + 2.0, 0.75);
			propParams.kappa1AnyUpdates(true);
			if (stdDev>propParams.kappa1StdDevUpper() || stdDev<propParams.kappa1StdDevLower()) {
				propParams.kappa1StdDevReset();
			}
			propParams.kappa1LocalReset();
		}

	}

} 

// Adaptive Metropolis Hastings move for kappa1 when the separation prior is used 
void metropolisHastingsForKappa1SP(mcmcChain<pReMiuMParams>& chain,
	unsigned int& nTry, unsigned int& nAccept,
	const mcmcModel<pReMiuMParams,
	pReMiuMOptions,
	pReMiuMData>& model,
	pReMiuMPropParams& propParams,
	baseGeneratorType& rndGenerator) {

	mcmcState<pReMiuMParams>& currentState = chain.currentState();
	pReMiuMParams& currentParams = currentState.parameters();
	pReMiuMHyperParams hyperParams = currentParams.hyperParams();

	// Find the number of clusters
	unsigned int maxZ = currentParams.workMaxZi();
	// Find the number of covariates
	unsigned int nCovariates = 0;
	if (model.options().covariateType().compare("Mixed") == 0) {
		nCovariates = currentParams.nContinuousCovs();
	}
	else {
		nCovariates = currentParams.nCovariates();
	}

	// Define a uniform random number generator
	randomUniform unifRand(0, 1);

	double& stdDev = propParams.kappa1StdDev();
	double kappa1Current = currentParams.kappa11();

	double kappa1Prop;
	kappa1Prop = truncNormalRand(rndGenerator, kappa1Current, stdDev, "L", nCovariates, 0);

	double logAcceptRatio = 0.0;
	double workLogDetR0 = hyperParams.workLogDetR0();
	for (unsigned int c = 0; c <= maxZ; c++) {
		double workLogDetTauR = currentParams.workLogDetTauR(c);
		logAcceptRatio += -logMultivarGammaFn(kappa1Prop / 2.0, nCovariates) +
			(-(kappa1Prop*nCovariates) / 2.0)*log(2.0) + (-kappa1Prop / 2.0)*workLogDetR0 +
			((kappa1Prop - nCovariates - 1.0) / 2.0)*workLogDetTauR;
		logAcceptRatio -= -logMultivarGammaFn(kappa1Current / 2.0, nCovariates) +
			(-(kappa1Current*nCovariates) / 2.0)*log(2.0) + (-kappa1Current / 2.0)*workLogDetR0 +
			((kappa1Current - nCovariates - 1.0) / 2.0)*workLogDetTauR;
	}

	logAcceptRatio += logPdfInverseGamma(kappa1Prop - nCovariates, hyperParams.shapeKappa1(), hyperParams.scaleKappa1());
	logAcceptRatio -= logPdfInverseGamma(kappa1Current - nCovariates, hyperParams.shapeKappa1(), hyperParams.scaleKappa1());

	// Add the proposal contribution
	logAcceptRatio += logPdfTruncatedNormal(kappa1Current, kappa1Prop, stdDev, "L", nCovariates, 0);
	logAcceptRatio -= logPdfTruncatedNormal(kappa1Prop, kappa1Current, stdDev, "L", nCovariates, 0);

	propParams.kappa1AddTry();
	nTry++;
	if (unifRand(rndGenerator)<exp(logAcceptRatio)) {
		nAccept++;
		propParams.kappa1AddAccept();
		// If the move was accepted update the state
		currentParams.kappa11(kappa1Prop);

		// Also update the proposal standard deviation
		if (propParams.nTryKappa1() % propParams.kappa1UpdateFreq() == 0) {
			stdDev += 10 * (propParams.kappa1LocalAcceptRate() - propParams.kappa1AcceptTarget()) /
				pow((double)(propParams.nTryKappa1() / propParams.kappa1UpdateFreq()) + 2.0, 0.75);
			propParams.kappa1AnyUpdates(true);
			if (stdDev>propParams.kappa1StdDevUpper() || stdDev<propParams.kappa1StdDevLower()) {
				propParams.kappa1StdDevReset();
			}
			propParams.kappa1LocalReset();
		}
	}
	else {
		// Otherwise update the proposal standard deviation
		if (propParams.nTryKappa1() % propParams.kappa1UpdateFreq() == 0) {
			stdDev += 10 * (propParams.kappa1LocalAcceptRate() - propParams.kappa1AcceptTarget()) /
				pow((double)(propParams.nTryKappa1() / propParams.kappa1UpdateFreq()) + 2.0, 0.75);
			propParams.kappa1AnyUpdates(true);
			if (stdDev>propParams.kappa1StdDevUpper() || stdDev<propParams.kappa1StdDevLower()) {
				propParams.kappa1StdDevReset();
			}
			propParams.kappa1LocalReset();
		}

	}

}



// Gibbs move for updating R1
void gibbsForR1(mcmcChain<pReMiuMParams>& chain,
				unsigned int& nTry,unsigned int& nAccept,
				const mcmcModel<pReMiuMParams,pReMiuMOptions,pReMiuMData>& model,
				pReMiuMPropParams& propParams,
				baseGeneratorType& rndGenerator){

	mcmcState<pReMiuMParams>& currentState = chain.currentState();
	pReMiuMParams& currentParams = currentState.parameters();
	pReMiuMHyperParams hyperParams = currentParams.hyperParams();

	// Find the number of clusters
	unsigned int maxZ = currentParams.workMaxZi();
	// Find the number of covariates
	unsigned int nCovariates = 0;
	if(model.options().covariateType().compare("Mixed")==0){
		nCovariates = currentParams.nContinuousCovs();
	} else {
		nCovariates = currentParams.nCovariates();
	}

	nTry++;
	nAccept++;


	MatrixXd SumTau = MatrixXd::Zero(nCovariates,nCovariates);
	unsigned int workNactive=0;
	for (unsigned int c=0; c<maxZ+1;c++){
		SumTau += currentParams.Tau(c);
		workNactive += 1;
	}
	SumTau += hyperParams.R0();
	MatrixXd R0Star = SumTau.inverse();

	MatrixXd inverseR1 =  wishartRand(rndGenerator,R0Star,workNactive*currentParams.kappa11()+hyperParams.kappa0());
	MatrixXd R1 = inverseR1.inverse();
	currentParams.R1(R1);
	
}

// Gibbs update for mu00 in Normal covariate case
void gibbsForMu00(mcmcChain<pReMiuMParams>& chain,
	unsigned int& nTry, unsigned int& nAccept,
	const mcmcModel<pReMiuMParams, pReMiuMOptions, pReMiuMData>& model,
	pReMiuMPropParams& propParams,
	baseGeneratorType& rndGenerator) {

	mcmcState<pReMiuMParams>& currentState = chain.currentState();
	pReMiuMParams& currentParams = currentState.parameters();
	pReMiuMHyperParams hyperParams = currentParams.hyperParams();
	bool useSeparationPrior = model.options().useSeparationPrior();

	// Find the number of clusters
	unsigned int maxZ = currentParams.workMaxZi();
	// Find the number of covariates
	unsigned int nCovariates = 0;
	if (model.options().covariateType().compare("Mixed") == 0) {
		nCovariates = currentParams.nContinuousCovs();
	}
	else {
		nCovariates = currentParams.nCovariates();
	}

	nTry++;
	nAccept++;
	
	VectorXd SumMu;
	SumMu.setZero(nCovariates);
	unsigned int workNactive = 0;
	for (unsigned int c = 0; c <= maxZ; c++) {
		SumMu = SumMu + currentParams.mu(c);
		workNactive += 1;
	}

	MatrixXd sigmaMu00(nCovariates, nCovariates);
	VectorXd meanMu00(nCovariates);

	if (useSeparationPrior) {
		sigmaMu00 = (workNactive*hyperParams.Tau00() + hyperParams.Tau0()).inverse();
		meanMu00 = sigmaMu00*(hyperParams.Tau00()*SumMu + hyperParams.Tau0()*hyperParams.mu0());
	}
	else {
		sigmaMu00 = (workNactive*currentParams.Tau00() + hyperParams.Tau0()).inverse();
		meanMu00 = sigmaMu00*(currentParams.Tau00()*SumMu + hyperParams.Tau0()*hyperParams.mu0());
	}
	

	VectorXd mu00(nCovariates);
	mu00 = multivarNormalRand(rndGenerator, meanMu00, sigmaMu00);

	currentParams.mu00(mu00);

}

// Gibbs update for Tau00 in Normal covariate case
void gibbsForTau00(mcmcChain<pReMiuMParams>& chain,
	unsigned int& nTry, unsigned int& nAccept,
	const mcmcModel<pReMiuMParams, pReMiuMOptions, pReMiuMData>& model,
	pReMiuMPropParams& propParams,
	baseGeneratorType& rndGenerator) {

	mcmcState<pReMiuMParams>& currentState = chain.currentState();
	pReMiuMParams& currentParams = currentState.parameters();
	pReMiuMHyperParams hyperParams = currentParams.hyperParams();

	// Find the number of clusters
	unsigned int maxZ = currentParams.workMaxZi();
	// Find the number of covariates
	unsigned int nCovariates = 0;
	if (model.options().covariateType().compare("Mixed") == 0) {
		nCovariates = currentParams.nContinuousCovs();
	}
	else {
		nCovariates = currentParams.nCovariates();
	}

	nTry++;
	nAccept++;

	MatrixXd SumMuDiff;
	SumMuDiff.setZero(nCovariates, nCovariates);
	unsigned int workNactive = 0;
	for (unsigned int c = 0; c <= maxZ; c++) {
		SumMuDiff = SumMuDiff + (currentParams.mu(c) - currentParams.mu00())*((currentParams.mu(c) - currentParams.mu00()).transpose());
		workNactive += 1;
	}

	SumMuDiff += hyperParams.R00().inverse();

	MatrixXd RUpadated(nCovariates, nCovariates);
	RUpadated= SumMuDiff.inverse();

	MatrixXd Tau00(nCovariates, nCovariates);
	Tau00 = wishartRand(rndGenerator, RUpadated, workNactive + hyperParams.kappa00());

	currentParams.Tau00(Tau00);

}

// Gibbs move for updating BetaTauS when the separation prior is used 
void gibbsForBetaTauS(mcmcChain<pReMiuMParams>& chain,
	unsigned int& nTry, unsigned int& nAccept,
	const mcmcModel<pReMiuMParams, pReMiuMOptions, pReMiuMData>& model,
	pReMiuMPropParams& propParams,
	baseGeneratorType& rndGenerator) {

	mcmcState<pReMiuMParams>& currentState = chain.currentState();
	pReMiuMParams& currentParams = currentState.parameters();
	pReMiuMHyperParams hyperParams = currentParams.hyperParams();

	// Find the number of clusters
	unsigned int maxZ = currentParams.workMaxZi();
	// Find the number of covariates
	unsigned int nCovariates = 0;
	if (model.options().covariateType().compare("Mixed") == 0) {
		nCovariates = currentParams.nContinuousCovs();
	}
	else {
		nCovariates = currentParams.nCovariates();
	}

	nTry++;
	nAccept++;

	VectorXd sumTau = VectorXd::Zero(nCovariates);
	unsigned int workNactive = 0;
	for (unsigned int c = 0; c < maxZ + 1; c++) {
		for (unsigned int j = 0; j < nCovariates; j++) {
			sumTau(j) = sumTau(j) + currentParams.TauS(c,j);
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
	currentParams.beta_taus(betas_new);

}


// Gibbs move for updating R1 when the independent normal conditional likelihood is used 
void gibbsForR1Indep(mcmcChain<pReMiuMParams>& chain,
	unsigned int& nTry, unsigned int& nAccept,
	const mcmcModel<pReMiuMParams, pReMiuMOptions, pReMiuMData>& model,
	pReMiuMPropParams& propParams,
	baseGeneratorType& rndGenerator) {

	mcmcState<pReMiuMParams>& currentState = chain.currentState();
	pReMiuMParams& currentParams = currentState.parameters();
	pReMiuMHyperParams hyperParams = currentParams.hyperParams();

	// Find the number of clusters
	unsigned int maxZ = currentParams.workMaxZi();
	// Find the number of covariates
	unsigned int nCovariates = 0;
	if (model.options().covariateType().compare("Mixed") == 0) {
		nCovariates = currentParams.nContinuousCovs();
	}
	else {
		nCovariates = currentParams.nCovariates();
	}

	nTry++;
	nAccept++;

	VectorXd sumTau= VectorXd::Zero(nCovariates);
	unsigned int workNactive = 0;
	for (unsigned int c = 0; c < maxZ + 1; c++) {
		sumTau = sumTau + currentParams.Tau_Indep(c);
		workNactive += 1;
	}

	VectorXd r1(nCovariates);
	VectorXd R0_Indep = hyperParams.R0_Indep();
	for (unsigned int j = 0; j < nCovariates; j++) {
		double kappaNew = workNactive * (double)hyperParams.kappa1() + (double)hyperParams.kappa0();
		double rNew = sumTau(j) + R0_Indep(j);

		randomGamma gammaRand(kappaNew, 1.0 / rNew);
		r1(j) = gammaRand(rndGenerator);

	}

	currentParams.R1_Indep(r1);

}

// Gibbs move for updating phi
void gibbsForPhiInActive(mcmcChain<pReMiuMParams>& chain,
				unsigned int& nTry,unsigned int& nAccept,
				const mcmcModel<pReMiuMParams,pReMiuMOptions,pReMiuMData>& model,
				pReMiuMPropParams& propParams,
				baseGeneratorType& rndGenerator){

	mcmcState<pReMiuMParams>& currentState = chain.currentState();
	pReMiuMParams& currentParams = currentState.parameters();
	pReMiuMHyperParams hyperParams = currentParams.hyperParams();

	string varSelectType = model.options().varSelectType();
	// Find the number of clusters
	unsigned int maxZ = currentParams.workMaxZi();
	unsigned int maxNClusters = currentParams.maxNClusters();

	// Find the number of covariates
	unsigned int nCovariates = 0;
	if(model.options().covariateType().compare("Mixed")==0){
		nCovariates = currentParams.nDiscreteCovs();
	} else {
		nCovariates = currentParams.nCovariates();
	}

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
void gibbsForMuInActive(mcmcChain<pReMiuMParams>& chain,
				unsigned int& nTry,unsigned int& nAccept,
				const mcmcModel<pReMiuMParams,pReMiuMOptions,pReMiuMData>& model,
				pReMiuMPropParams& propParams,
				baseGeneratorType& rndGenerator){

	mcmcState<pReMiuMParams>& currentState = chain.currentState();
	pReMiuMParams& currentParams = currentState.parameters();
	pReMiuMHyperParams hyperParams = currentParams.hyperParams();

	bool useIndependentNormal = model.options().useIndependentNormal();
	bool useHyperpriorR1 = model.options().useHyperpriorR1();
	bool useSeparationPrior = model.options().useSeparationPrior();


	// Find the number of clusters
	unsigned int maxZ = currentParams.workMaxZi();
	unsigned int maxNClusters = currentParams.maxNClusters();
	// Find the number of covariates
	unsigned int nCovariates = 0;
	if(model.options().covariateType().compare("Mixed")==0){
		nCovariates = currentParams.nContinuousCovs();
	} else {
		nCovariates = currentParams.nCovariates();
	}

	nTry++;
	nAccept++;

	MatrixXd covMat(nCovariates,nCovariates);
	VectorXd meanVec(nCovariates);

	if (useHyperpriorR1) {
		covMat = currentParams.Tau00().inverse();
		meanVec = currentParams.mu00();
	}
	else if (useSeparationPrior) {
		covMat = hyperParams.Tau00().inverse();
		meanVec = currentParams.mu00();
	} 
	else {
		covMat = hyperParams.Tau0().inverse();
		meanVec = hyperParams.mu0();
	}

	for(unsigned int c=maxZ+1;c<maxNClusters;c++){
		VectorXd mu(nCovariates);

		// We sample from this posterior
		mu = multivarNormalRand(rndGenerator, meanVec, covMat);

		// We store our sample
		currentParams.mu(c, mu, useIndependentNormal);
	}

}



// Gibbs update for mu in Independent Normal covariate case
void gibbsForMuInActiveIndep(mcmcChain<pReMiuMParams>& chain,
	unsigned int& nTry, unsigned int& nAccept,
	const mcmcModel<pReMiuMParams, pReMiuMOptions, pReMiuMData>& model,
	pReMiuMPropParams& propParams,
	baseGeneratorType& rndGenerator) {

	mcmcState<pReMiuMParams>& currentState = chain.currentState();
	pReMiuMParams& currentParams = currentState.parameters();
	pReMiuMHyperParams hyperParams = currentParams.hyperParams();

	bool useIndependentNormal = model.options().useIndependentNormal();

	// Find the number of clusters
	unsigned int maxZ = currentParams.workMaxZi();
	unsigned int maxNClusters = currentParams.maxNClusters();
	// Find the number of covariates
	unsigned int nCovariates = 0;
	if (model.options().covariateType().compare("Mixed") == 0) {
		nCovariates = currentParams.nContinuousCovs();
	}
	else {
		nCovariates = currentParams.nCovariates();
	}

	nTry++;
	nAccept++;

	VectorXd mu0(nCovariates);
	mu0 = hyperParams.mu0();
	VectorXd Tau0(nCovariates);
	Tau0 = hyperParams.Tau0_Indep();

	for (unsigned int c = maxZ + 1; c<maxNClusters; c++) {
		VectorXd mu(nCovariates);
		for (unsigned int j = 0; j < nCovariates; j++) {
			double mean = mu0(j);
			double variance = 1.0 / Tau0(j);
			mu(j) = NormalRand(rndGenerator, mean, variance);	
		}
		currentParams.mu(c, mu, useIndependentNormal);
	}

}


// Gibbs update for mu in Normal covariate case when the normal inverse Wishart prior is used
void gibbsForMuInActiveNIWP(mcmcChain<pReMiuMParams>& chain,
				unsigned int& nTry,unsigned int& nAccept,
				const mcmcModel<pReMiuMParams,pReMiuMOptions,pReMiuMData>& model,
				pReMiuMPropParams& propParams,
				baseGeneratorType& rndGenerator){

	mcmcState<pReMiuMParams>& currentState = chain.currentState();
	pReMiuMParams& currentParams = currentState.parameters();
	pReMiuMHyperParams hyperParams = currentParams.hyperParams();

	bool useIndependentNormal = model.options().useIndependentNormal();

	// Find the number of clusters
	unsigned int maxZ = currentParams.workMaxZi();
	unsigned int maxNClusters = currentParams.maxNClusters();
	// Find the number of covariates
	unsigned int nCovariates = 0;
	if(model.options().covariateType().compare("Mixed")==0){
		nCovariates = currentParams.nContinuousCovs();
	} else {
		nCovariates = currentParams.nCovariates();
	}

	nTry++;
	nAccept++;


	VectorXd meanVec(nCovariates);
	meanVec = hyperParams.mu0();

	for(unsigned int c=maxZ+1;c<maxNClusters;c++){
        MatrixXd covMat(nCovariates,nCovariates);
        covMat = currentParams.Sigma(c)/hyperParams.nu0();

		VectorXd mu(nCovariates);
		// We sample from this posterior
		mu = multivarNormalRand(rndGenerator,meanVec,covMat);
		// We store our sample
		currentParams.mu(c, mu, useIndependentNormal);
	}

}

// Gibbs update for Tau in the Normal covariate case
void gibbsForTauInActive(mcmcChain<pReMiuMParams>& chain,
				unsigned int& nTry,unsigned int& nAccept,
				const mcmcModel<pReMiuMParams,pReMiuMOptions,pReMiuMData>& model,
				pReMiuMPropParams& propParams,
				baseGeneratorType& rndGenerator){

	mcmcState<pReMiuMParams>& currentState = chain.currentState();
	pReMiuMParams& currentParams = currentState.parameters();
	pReMiuMHyperParams hyperParams = currentParams.hyperParams();
	bool useHyperpriorR1 = model.options().useHyperpriorR1();

	// Find the number of clusters
	unsigned int maxZ = currentParams.workMaxZi();
	unsigned int maxNClusters = currentParams.maxNClusters();

	nTry++;
	nAccept++;

	if (useHyperpriorR1){
		for(unsigned int c=maxZ+1;c<maxNClusters;c++){
			MatrixXd Tau = wishartRand(rndGenerator,currentParams.R1(), currentParams.kappa11());
			currentParams.Tau(c,Tau);
		}
	} else {
		for(unsigned int c=maxZ+1;c<maxNClusters;c++){
			MatrixXd Tau = wishartRand(rndGenerator,hyperParams.R0(),hyperParams.kappa0());
			currentParams.Tau(c,Tau);
		}
	}


}


// Gibbs update for TauR in the Normal covariate case when separation prior is used 
void gibbsForTauRInActive(mcmcChain<pReMiuMParams>& chain,
	unsigned int& nTry, unsigned int& nAccept,
	const mcmcModel<pReMiuMParams, pReMiuMOptions, pReMiuMData>& model,
	pReMiuMPropParams& propParams,
	baseGeneratorType& rndGenerator) {

	mcmcState<pReMiuMParams>& currentState = chain.currentState();
	pReMiuMParams& currentParams = currentState.parameters();
	pReMiuMHyperParams hyperParams = currentParams.hyperParams();

	// Find the number of clusters
	unsigned int maxZ = currentParams.workMaxZi();
	unsigned int maxNClusters = currentParams.maxNClusters();

	nTry++;
	nAccept++;

	for (unsigned int c = maxZ + 1; c<maxNClusters; c++) {
		MatrixXd TauR = wishartRand(rndGenerator, hyperParams.R0(), currentParams.kappa11());
		currentParams.TauR(c, TauR);
	}

}


// Gibbs update for TauS in the Normal covariate case when separation prior is used 
void gibbsForTauSInActive(mcmcChain<pReMiuMParams>& chain,
	unsigned int& nTry, unsigned int& nAccept,
	const mcmcModel<pReMiuMParams, pReMiuMOptions, pReMiuMData>& model,
	pReMiuMPropParams& propParams,
	baseGeneratorType& rndGenerator) {

	mcmcState<pReMiuMParams>& currentState = chain.currentState();
	pReMiuMParams& currentParams = currentState.parameters();
	pReMiuMHyperParams hyperParams = currentParams.hyperParams();

	// Find the number of clusters
	unsigned int maxZ = currentParams.workMaxZi();
	unsigned int maxNClusters = currentParams.maxNClusters();

	unsigned int nCovariates = 0;
	if (model.options().covariateType().compare("Mixed") == 0) {
		nCovariates = currentParams.nContinuousCovs();
	}
	else {
		nCovariates = currentParams.nCovariates();
	}

	nTry++;
	nAccept++;

	for (unsigned int c = maxZ + 1; c<maxNClusters; c++) {
		for (unsigned int j = 0; j < nCovariates; j++) {
			double alpha_taus = hyperParams.alpha_taus();
			double beta_tausj = currentParams.beta_taus(j);
			randomGamma gammaRand(alpha_taus, 1.0 / beta_tausj);
			double tausj = gammaRand(rndGenerator);
			currentParams.TauS(c, j, tausj);
		}
	}

}


//Gibbs update for Tau in the Independent Normal covariate case
void gibbsForTauInActiveIndep(mcmcChain<pReMiuMParams>& chain,
	unsigned int& nTry, unsigned int& nAccept,
	const mcmcModel<pReMiuMParams, pReMiuMOptions, pReMiuMData>& model,
	pReMiuMPropParams& propParams,
	baseGeneratorType& rndGenerator) {

	mcmcState<pReMiuMParams>& currentState = chain.currentState();
	pReMiuMParams& currentParams = currentState.parameters();
	pReMiuMHyperParams hyperParams = currentParams.hyperParams();

	// Find the number of clusters
	unsigned int maxZ = currentParams.workMaxZi();
	unsigned int maxNClusters = currentParams.maxNClusters();

	unsigned int nCovariates = 0;
	if (model.options().covariateType().compare("Mixed") == 0) {
		nCovariates = currentParams.nContinuousCovs();
	}
	else {
		nCovariates = currentParams.nCovariates();
	}

	nTry++;
	nAccept++;

	VectorXd Tau(nCovariates);
	for (unsigned int c = maxZ + 1; c<maxNClusters; c++) {
		for (unsigned int j = 0; j < nCovariates; j++) {
			double kappa = hyperParams.kappa1();
			double r = currentParams.R1_Indep(j);
			randomGamma gammaRand(kappa, 1.0 / r);
			Tau(j) = gammaRand(rndGenerator);
		}
		currentParams.Tau_Indep(c, Tau);
	}

}

// Gibbs update for update of gamma (only used in the binary variable selection case)
void gibbsForGammaInActive(mcmcChain<pReMiuMParams>& chain,
					unsigned int& nTry,unsigned int& nAccept,
					const mcmcModel<pReMiuMParams,
									pReMiuMOptions,
									pReMiuMData>& model,
					pReMiuMPropParams& propParams,
					baseGeneratorType& rndGenerator){


	mcmcState<pReMiuMParams>& currentState = chain.currentState();
	pReMiuMParams& currentParams = currentState.parameters();
	pReMiuMHyperParams hyperParams = currentParams.hyperParams();

	// Find the number of subjects
	unsigned int nCovariates = currentParams.nCovariates();
	unsigned int maxZ = currentParams.workMaxZi();
	unsigned int maxNClusters = currentParams.maxNClusters();
	string covariateType = model.options().covariateType();
	string varSelectType = model.options().varSelectType();
	bool useIndependentNormal = model.options().useIndependentNormal();

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
				currentParams.gamma(c,j,proposedGamma,covariateType, useIndependentNormal);

			}


		}
	}

}

// Gibbs for theta
void gibbsForThetaInActive(mcmcChain<pReMiuMParams>& chain,
								unsigned int& nTry,unsigned int& nAccept,
								const mcmcModel<pReMiuMParams,
												pReMiuMOptions,
												pReMiuMData>& model,
								pReMiuMPropParams& propParams,
								baseGeneratorType& rndGenerator){

	mcmcState<pReMiuMParams>& currentState = chain.currentState();
	pReMiuMParams& currentParams = currentState.parameters();
	pReMiuMHyperParams hyperParams = currentParams.hyperParams();
	const pReMiuMData& dataset = model.dataset();
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

// Gibbs for nu inactive (for survival case)
void gibbsForNuInActive(mcmcChain<pReMiuMParams>& chain,
								unsigned int& nTry,unsigned int& nAccept,
								const mcmcModel<pReMiuMParams,
												pReMiuMOptions,
												pReMiuMData>& model,
								pReMiuMPropParams& propParams,
								baseGeneratorType& rndGenerator){

	mcmcState<pReMiuMParams>& currentState = chain.currentState();
	pReMiuMParams& currentParams = currentState.parameters();
	pReMiuMHyperParams hyperParams = currentParams.hyperParams();
	const string outcomeType = model.dataset().outcomeType();

	// Find the number of clusters
	unsigned int maxZ = currentParams.workMaxZi();
	unsigned int maxNClusters = currentParams.maxNClusters();

	nTry++;
	nAccept++;

	randomGamma gammaRand(hyperParams.shapeNu(),hyperParams.scaleNu());
	for(unsigned int c=maxZ+1;c<maxNClusters;c++){
		double nu=gammaRand(rndGenerator);
		currentParams.nu(c,nu);
	}
}

/*********** BLOCK 4 p(Theta^N|.) **********************************/
// N=Non-cluster, and Theta contains: beta, rho, omega, lambda, tau_epsilon, uCAR and TauCAR

// Adaptive Metropolis-Hastings for beta
void metropolisHastingsForBeta(mcmcChain<pReMiuMParams>& chain,
								unsigned int& nTry,unsigned int& nAccept,
								const mcmcModel<pReMiuMParams,
												pReMiuMOptions,
												pReMiuMData>& model,
								pReMiuMPropParams& propParams,
								baseGeneratorType& rndGenerator){

	mcmcState<pReMiuMParams>& currentState = chain.currentState();
	pReMiuMParams& currentParams = currentState.parameters();
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
void metropolisHastingsForLambda(mcmcChain<pReMiuMParams>& chain,
								unsigned int& nTry,unsigned int& nAccept,
								const mcmcModel<pReMiuMParams,
												pReMiuMOptions,
												pReMiuMData>& model,
								pReMiuMPropParams& propParams,
								baseGeneratorType& rndGenerator){

	mcmcState<pReMiuMParams>& currentState = chain.currentState();
	pReMiuMParams& currentParams = currentState.parameters();
	const string outcomeType = model.dataset().outcomeType();

	// Find the number of subjects
	unsigned int nSubjects = currentParams.nSubjects();

	// Define a uniform random number generator
	randomUniform unifRand(0,1);
	// Define a normal random number generator
	randomNormal normRand(0,1);


	double lambdaTargetRate = propParams.lambdaAcceptTarget();
	unsigned int lambdaUpdateFreq = propParams.lambdaUpdateFreq();

	double (*logCondPostLambdai)(const pReMiuMParams&,
											const mcmcModel<pReMiuMParams,
															pReMiuMOptions,
															pReMiuMData>&,
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
void gibbsForTauEpsilon(mcmcChain<pReMiuMParams>& chain,
						unsigned int& nTry,unsigned int& nAccept,
						const mcmcModel<pReMiuMParams,
										pReMiuMOptions,
										pReMiuMData>& model,
						pReMiuMPropParams& propParams,
						baseGeneratorType& rndGenerator){

	mcmcState<pReMiuMParams>& currentState = chain.currentState();
	pReMiuMParams& currentParams = currentState.parameters();
	pReMiuMHyperParams hyperParams = currentParams.hyperParams();
	const pReMiuMData& dataset = model.dataset();
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
void metropolisHastingsForRhoOmega(mcmcChain<pReMiuMParams>& chain,
									unsigned int& nTry,unsigned int& nAccept,
									const mcmcModel<pReMiuMParams,
													pReMiuMOptions,
													pReMiuMData>& model,
									pReMiuMPropParams& propParams,
									baseGeneratorType& rndGenerator){

	mcmcState<pReMiuMParams>& currentState = chain.currentState();
	pReMiuMParams& currentParams = currentState.parameters();
	pReMiuMHyperParams hyperParams = currentParams.hyperParams();

	// Find the number of subjects
	unsigned int nCovariates = currentParams.nCovariates();
	string covariateType = model.options().covariateType();
	string varSelectType = model.options().varSelectType();
	bool useIndependentNormal = model.options().useIndependentNormal();

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
		if(unifRand(rndGenerator)>hyperParams.atomRho()){
			// Proposing an omega 0
			if(currentOmega[j]==0){
				// Nothing to do as move to the same place
				nAccept++;
				continue;
			}
			proposedOmega=0;
			proposedRho=0.0;

			currentParams.omega(j,proposedOmega);
			currentParams.rho(j, proposedRho, covariateType, varSelectType, useIndependentNormal);
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
				currentParams.rho(j, currentRho[j], covariateType, varSelectType, useIndependentNormal);
			}

		}else{
			if(currentOmega[j]==1){
				proposedRho  = truncNormalRand(rndGenerator,currentRho[j],stdDev,"B",0,1);
				currentParams.rho(j, proposedRho, covariateType, varSelectType, useIndependentNormal);
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
					currentParams.rho(j, currentRho[j], covariateType, varSelectType, useIndependentNormal);
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
				currentParams.rho(j,proposedRho,covariateType,varSelectType, useIndependentNormal);
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
					currentParams.rho(j, currentRho[j], covariateType, varSelectType, useIndependentNormal);
				}
			}
		}


	}

}

// Gibbs for update of sigmaSqY (Normal response case)
void gibbsForSigmaSqY(mcmcChain<pReMiuMParams>& chain,
						unsigned int& nTry,unsigned int& nAccept,
						const mcmcModel<pReMiuMParams,
										pReMiuMOptions,
										pReMiuMData>& model,
						pReMiuMPropParams& propParams,
						baseGeneratorType& rndGenerator){
	mcmcState<pReMiuMParams>& currentState = chain.currentState();
	pReMiuMParams& currentParams = currentState.parameters();
	pReMiuMHyperParams hyperParams = currentParams.hyperParams();

	const pReMiuMData& dataset = model.dataset();

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

// Gibbs for update of sigmaSqYQuantile (Quantile response case)
void gibbsForSigmaSqYQuantile(mcmcChain<pReMiuMParams>& chain,
				unsigned int& nTry,unsigned int& nAccept,
				const mcmcModel<pReMiuMParams, pReMiuMOptions, pReMiuMData>& model, pReMiuMPropParams& propParams, baseGeneratorType& rndGenerator){
	mcmcState<pReMiuMParams>& currentState = chain.currentState();
	pReMiuMParams& currentParams = currentState.parameters();
	pReMiuMHyperParams hyperParams = currentParams.hyperParams();

	const pReMiuMData& dataset = model.dataset();

	unsigned int nSubjects = currentParams.nSubjects();
	unsigned int nFixedEffects = dataset.nFixedEffects();
	double pQuantile = hyperParams.pQuantile();

	nTry++;
	nAccept++;

	double sumVal=0.0;
	for(unsigned int i=0;i<nSubjects;i++){
		int Zi = currentParams.z(i);
		double mu = currentParams.theta(Zi,0);
		for(unsigned int j=0;j<nFixedEffects;j++){
			mu+=currentParams.beta(j,0)*dataset.W(i,j);
		}
		// here used different parametrisation of ALD than in paper (see distribution.h for further details on this parametrisation)
		sumVal+=(abs(dataset.continuousY(i)-mu)+(2*pQuantile-1)*(dataset.continuousY(i)-mu))/2;
	}
      
	double posteriorShape = hyperParams.shapeSigmaSqY()+(double)nSubjects;
	double posteriorScale = hyperParams.scaleSigmaSqY()+sumVal;
        
	randomGamma gammaRand(posteriorShape,1.0/posteriorScale);
	double sigmaSqY=1.0/gammaRand(rndGenerator);
	currentParams.sigmaSqY(sigmaSqY);

}


// Gibbs for update of nu (survival response case) using adaptive rejection sampling
void gibbsForNu(mcmcChain<pReMiuMParams>& chain,
		unsigned int& nTry,unsigned int& nAccept,
		const mcmcModel<pReMiuMParams,
		pReMiuMOptions,
		pReMiuMData>& model,
		pReMiuMPropParams& propParams,
		baseGeneratorType& rndGenerator){

	mcmcState<pReMiuMParams>& currentState = chain.currentState();
	pReMiuMParams& currentParams = currentState.parameters();
	pReMiuMHyperParams hyperParams = currentParams.hyperParams();
	const bool weibullFixedShape=model.options().weibullFixedShape();
	// Find the number of clusters
	unsigned int maxZ = currentParams.workMaxZi();

	nTry++;
	nAccept++;
	if (weibullFixedShape){
		double nu = ARSsampleNu(currentParams, model, 0,logNuPostSurvival,rndGenerator);
		currentParams.nu(0,nu);
	} else {
		for (unsigned int c=0;c<=maxZ;c++){
			double nu = ARSsampleNu(currentParams, model, c,logNuPostSurvival,rndGenerator);
			currentParams.nu(c,nu);
		}
	}
}

// Gibbs update for the precision of spatial random term
void gibbsForTauCAR(mcmcChain<pReMiuMParams>& chain,
						unsigned int& nTry,unsigned int& nAccept,
						const mcmcModel<pReMiuMParams,
										pReMiuMOptions,
										pReMiuMData>& model,
						pReMiuMPropParams& propParams,
						baseGeneratorType& rndGenerator){

	mcmcState<pReMiuMParams>& currentState = chain.currentState();
	pReMiuMParams& currentParams = currentState.parameters();
	pReMiuMHyperParams hyperParams = currentParams.hyperParams();
	const pReMiuMData& dataset = model.dataset();

	//Rprintf("TauCAR before update is %f \n.", currentParams.TauCAR());

	unsigned int nSubjects=dataset.nSubjects();

	nTry++;
	nAccept++;

	double a=hyperParams.shapeTauCAR(),b=hyperParams.rateTauCAR();

	double sumCAR1 = 0.0;
	double sumCAR2 = 0.0;
	for (unsigned int i=0; i<nSubjects; i++){
		double uCARi = currentParams.uCAR(i);
		int nNeighi = dataset.nNeighbours(i);
		sumCAR1+= uCARi*uCARi*nNeighi;
		for (int j = 0; j<nNeighi; j++){
			unsigned int nj = dataset.neighbours(i,j);
			double ucarj = currentParams.uCAR(nj-1);
			sumCAR2+=uCARi*ucarj;
	        }
	}
	double sumCAR=sumCAR1-sumCAR2;

	a+=(double)(nSubjects-1)/2.0;
	b+=sumCAR/2.0;


	// Boost uses shape and scale parameterisation
	randomGamma gammaRand(a,1.0/b);
	double tau = gammaRand(rndGenerator);
	currentParams.TauCAR(tau);
	//Rprintf("TauCAR after update is %f \n .", currentParams.TauCAR());
}

// Gibbs update for spatial random term using adaptive rejection sampling for Poisson outcome
void adaptiveRejectionSamplerForUCARPoisson(mcmcChain<pReMiuMParams>& chain,
						unsigned int& nTry,unsigned int& nAccept,
						const mcmcModel<pReMiuMParams,
										pReMiuMOptions,
										pReMiuMData>& model,
						pReMiuMPropParams& propParams,
						baseGeneratorType& rndGenerator){

	mcmcState<pReMiuMParams>& currentState = chain.currentState();
	pReMiuMParams& currentParams = currentState.parameters();
	const pReMiuMData& dataset = model.dataset();
	unsigned int nSubjects=dataset.nSubjects();


	nTry++;
	nAccept++;

	vector<double> tempU;
	tempU.resize(nSubjects);
	for (unsigned int iSub=0; iSub<nSubjects; iSub++){
		double ui=ARSsampleCAR(currentParams, model, iSub,logUiPostPoissonSpatial,rndGenerator);
		tempU[iSub]=ui;
	}
	double meanU=0.0;
	for (unsigned int i=0; i<nSubjects; i++){meanU+=tempU[i];}
	meanU/=nSubjects;
	for (unsigned int i=0; i<nSubjects; i++){tempU[i]-=meanU;}
	currentParams.uCAR(tempU);

}

// Random Walk Metropolis for Poisson outcome with spatial random effect
void metropolisForUCARPoisson(mcmcChain<pReMiuMParams>& chain,
						unsigned int& nTry,unsigned int& nAccept,
						const mcmcModel<pReMiuMParams,
										pReMiuMOptions,
										pReMiuMData>& model,
						pReMiuMPropParams& propParams,
						baseGeneratorType& rndGenerator){

	mcmcState<pReMiuMParams>& currentState = chain.currentState();
	pReMiuMParams& currentParams = currentState.parameters();
	const pReMiuMData& dataset = model.dataset();
	unsigned int nSubjects=dataset.nSubjects();
	unsigned int nFixedEffects=dataset.nFixedEffects();

	// Define a uniform random number generator
	randomUniform unifRand(0,1);
	// Define a normal random number generator
	randomNormal normRand(0,1);

	double uCARTargetRate = propParams.uCARAcceptTarget();
	unsigned int uCARUpdateFreq = propParams.uCARUpdateFreq();

	vector<double> tempU;
	tempU.resize(nSubjects);
	for (unsigned int iSub=0; iSub<nSubjects; iSub++){
		nTry++;
		propParams.uCARAddTry();
		double& stdDev = propParams.uCARStdDev();
			
		// likelihood	
		double uCAROrig = currentParams.uCAR(iSub);
		double uCARProp = uCAROrig +stdDev*normRand(rndGenerator);
		int zi = currentParams.z(iSub);
		double currentCondLogPost = logPYiGivenZiWiPoissonSpatial(currentParams,dataset,nFixedEffects,zi,iSub);
		currentParams.uCAR(iSub,uCARProp);

		double propCondLogPost = logPYiGivenZiWiPoissonSpatial(currentParams,dataset,nFixedEffects,zi,iSub);
		// prior variance
		unsigned int nNeighi = dataset.nNeighbours(iSub);
		double priorVar = currentParams.TauCAR()/nNeighi;

		// prior mean
		double priorMean=0.0;
		for (unsigned int j = 0; j<nNeighi; j++){
        		unsigned int nj = dataset.neighbours(iSub,j);
        		double ucarj = currentParams.uCAR(nj-1);
        		priorMean+=ucarj;
		}
		priorMean/=nNeighi;

		double propPrior = 0.5*pow((uCARProp-priorMean),2)/sqrt(priorVar);
		double currentPrior = 0.5*pow((uCAROrig-priorMean),2)/sqrt(priorVar);

		double logAcceptRatio = propCondLogPost - currentCondLogPost - propPrior + currentPrior; // Duncan

		if(unifRand(rndGenerator)<exp(logAcceptRatio)){
			nAccept++;
			tempU[iSub]=uCARProp;
			propParams.uCARAddAccept();
			currentCondLogPost = propCondLogPost;
			// Update the std dev of the proposal
			if(propParams.nTryuCAR()%uCARUpdateFreq==0){
				stdDev += 10*(propParams.uCARLocalAcceptRate()-uCARTargetRate)/
					pow((double)(propParams.nTryuCAR()/uCARUpdateFreq)+2.0,0.75);
				propParams.uCARAnyUpdates(true);
				if(stdDev>propParams.uCARStdDevUpper()||stdDev<propParams.uCARStdDevLower()){
					propParams.uCARStdDevReset();
				}
				propParams.thetaLocalReset();
			}
			// make sure that the spatial random effects sum up to 0
			double meanU=0.0;
			for (unsigned int kk=0; kk<nSubjects; kk++){meanU+=tempU[kk];}
			meanU/=nSubjects;
			for (unsigned int kk=0; kk<nSubjects; kk++){tempU[kk]-=meanU;}
			currentParams.uCAR(tempU);
		}else{
			tempU[iSub]=uCAROrig;
			currentParams.uCAR(iSub,uCAROrig);
			// Update the std dev of the proposal
			if(propParams.nTryuCAR()%uCARUpdateFreq==0){
				stdDev += 10*(propParams.uCARLocalAcceptRate()-uCARTargetRate)/
					pow((double)(propParams.nTryuCAR()/uCARUpdateFreq)+2.0,0.75);
				propParams.uCARAnyUpdates(true);
				if(stdDev<propParams.uCARStdDevLower()||stdDev>propParams.uCARStdDevUpper()){
					propParams.uCARStdDevReset();
				}
				propParams.uCARLocalReset();
			}
		}
	}
	double meanU=0.0;
	for (unsigned int i=0; i<nSubjects; i++){meanU+=tempU[i];}
	meanU/=nSubjects;
	for (unsigned int i=0; i<nSubjects; i++){tempU[i]-=meanU;}
	currentParams.uCAR(tempU);		
}


// Gibbs update for spatial random for Normal case
void gibbsForUCARNormal(mcmcChain<pReMiuMParams>& chain,
						unsigned int& nTry,unsigned int& nAccept,
						const mcmcModel<pReMiuMParams,
										pReMiuMOptions,
										pReMiuMData>& model,
						pReMiuMPropParams& propParams,
						baseGeneratorType& rndGenerator){

	mcmcState<pReMiuMParams>& currentState = chain.currentState();
	pReMiuMParams& currentParams = currentState.parameters();
	const pReMiuMData& dataset = model.dataset();
	unsigned int nSubjects=dataset.nSubjects();
	unsigned int nFixedEffects=dataset.nFixedEffects();

	vector<double> tempU;
	tempU.resize(nSubjects);

	nTry++;
	nAccept++;
	for (unsigned int iSub=0; iSub<nSubjects; iSub++){
		int nNeighi = dataset.nNeighbours(iSub);
		double sigmaSqUCAR = 1/(1/currentParams.sigmaSqY()+currentParams.TauCAR()*nNeighi);
		int Zi = currentParams.z(iSub);
		double betaW = 0.0;
		for(unsigned int j=0;j<nFixedEffects;j++){
			betaW+=currentParams.beta(j,0)*dataset.W(iSub,j);
		}
			double meanUi=0.0;
		for (int j = 0; j<nNeighi; j++){
        		unsigned int nj = dataset.neighbours(iSub,j);
        		double ucarj = currentParams.uCAR(nj-1);
	        		meanUi+=ucarj;
		}
		meanUi/=nNeighi;	
		double mUCAR = 1/currentParams.sigmaSqY()*(dataset.continuousY(iSub)-currentParams.theta(Zi,0)-betaW)+currentParams.TauCAR()*nNeighi*meanUi;
		mUCAR = mUCAR * sigmaSqUCAR;
		randomNormal normRand(0,1);
		tempU[iSub]=sqrt(sigmaSqUCAR)*normRand(rndGenerator)+mUCAR;
	}
	double meanU=0.0;
	for (unsigned int i=0; i<nSubjects; i++){meanU+=tempU[i];}
	meanU/=nSubjects;
	for (unsigned int i=0; i<nSubjects; i++){tempU[i]-=meanU;}
	currentParams.uCAR(tempU);
}



/*********** BLOCK 5 p(Z|.) **********************************/

// Gibbs update for the allocation variables
void gibbsForZ(mcmcChain<pReMiuMParams>& chain,
				unsigned int& nTry,unsigned int& nAccept,
				const mcmcModel<pReMiuMParams,
								pReMiuMOptions,
								pReMiuMData>& model,
				pReMiuMPropParams& propParams,
				baseGeneratorType& rndGenerator){

	mcmcState<pReMiuMParams>& currentState = chain.currentState();
	pReMiuMParams& currentParams = currentState.parameters();
	pReMiuMHyperParams hyperParams = currentParams.hyperParams();
	const pReMiuMData& dataset = model.dataset();
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
	unsigned int nDiscreteCovs=dataset.nDiscreteCovs();
	unsigned int nContinuousCovs=dataset.nContinuousCovs();
	vector<unsigned int>nCategories=dataset.nCategories();
	const vector<vector<bool> >& missingX=dataset.missingX();
	bool includeResponse = model.options().includeResponse();
	bool responseExtraVar = model.options().responseExtraVar();
	const bool includeCAR=model.options().includeCAR();
	const string& predictType = model.options().predictType();
	bool useIndependentNormal = model.options().useIndependentNormal();


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
						if (useIndependentNormal) {
							VectorXd muStar = currentParams.workMuStar(c);
							VectorXd sigma_cj(nCovariates);
							for (unsigned int j = 0; j < nCovariates; j++) {
								sigma_cj(j)= sqrt(1.0 / currentParams.Tau_Indep(c, j));
							}
							for (unsigned int j = 0; j<nCovariates; j++) {
								logPXiGivenZi[i][c] += logPdfNormal(xi(j), muStar(j), sigma_cj(j));
							}
						}
						else {
							logPXiGivenZi[i][c] = logPdfMultivarNormal(nCovariates, xi, currentParams.workMuStar(c), currentParams.workSqrtTau(c), currentParams.workLogDetTau(c));
						}
						
					}
				}

			}
		}
		// For the predictive subjects we do not count missing data
		LLT<MatrixXd> llt;
		for(unsigned int i=nSubjects;i<nSubjects+nPredictSubjects;i++){
			logPXiGivenZi[i].resize(maxNClusters,0.0);

			unsigned int nNotMissing=dataset.nContinuousCovariatesNotMissing(i);

			for(unsigned int c=0;c<maxNClusters;c++){
				if(u[i]<testBound[c]){
					VectorXd workMuStar=currentParams.workMuStar(c);

					VectorXd xi=VectorXd::Zero(nNotMissing);
					VectorXd muStar=VectorXd::Zero(nNotMissing);
					VectorXd sigma_cj = VectorXd::Zero(nNotMissing);
					MatrixXd sqrtTau = MatrixXd::Zero(nNotMissing, nNotMissing);
					double logDetTau = 0.0;
					
					if(nNotMissing==nCovariates){
						muStar=workMuStar;
						if (useIndependentNormal) {
							for (unsigned int j = 0; j < nCovariates; j++) {
								sigma_cj(j) = sqrt(1.0 / currentParams.Tau_Indep(c,j));
							}
						}
						else {
							sqrtTau = currentParams.workSqrtTau(c);
							logDetTau = currentParams.workLogDetTau(c);
						}
						for(unsigned int j=0;j<nCovariates;j++){
							xi(j)=currentParams.workContinuousX(i,j);
						}
					}else{

						if (useIndependentNormal) {
							VectorXd workSigma = currentParams.Sigma_Indep(c);
							unsigned int j = 0;
							for (unsigned int j0 = 0; j0<nCovariates; j0++) {
								if (!missingX[i][nDiscreteCovs + j0]) {
									xi(j) = currentParams.workContinuousX(i, j0);
									muStar(j) = workMuStar(j0);
									sigma_cj(j) = sqrt(workSigma(j0));
									j++;
								}
							}
						}
						else {
							MatrixXd workSigma = currentParams.Sigma(c);
							MatrixXd Sigma = MatrixXd::Zero(nNotMissing, nNotMissing);
							MatrixXd Tau = MatrixXd::Zero(nNotMissing, nNotMissing);
							unsigned int j = 0;
							for (unsigned int j0 = 0; j0<nCovariates; j0++) {
								if (!missingX[i][nDiscreteCovs + j0]) {
									xi(j) = currentParams.workContinuousX(i, j0);
									muStar(j) = workMuStar(j0);
									unsigned int r = 0;
									for (unsigned int j1 = 0; j1<nCovariates; j1++) {
										if (!missingX[i][nDiscreteCovs + j1]) {
											Sigma(j, r) = workSigma(j0, j1);
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

					}
						
					if (useIndependentNormal) {
						for (unsigned int j = 0; j<nNotMissing; j++) {
							logPXiGivenZi[i][c] += logPdfNormal(xi(j), muStar(j), sigma_cj(j));
						}
					}
					else {
						logPXiGivenZi[i][c] = logPdfMultivarNormal(nNotMissing, xi, muStar, sqrtTau, logDetTau);
					}

				}
			}
		}

	}else if(covariateType.compare("Mixed")==0){
		for(unsigned int i=0;i<nSubjects;i++){
			VectorXd xi=VectorXd::Zero(nContinuousCovs);
			logPXiGivenZi[i].resize(maxNClusters,0);
			for(unsigned int c=0;c<maxNClusters;c++){
				if(u[i]<testBound[c]){
					if(currentParams.z(i)==(int)c){
						logPXiGivenZi[i][c]=currentParams.workLogPXiGivenZi(i);
					}else{
						logPXiGivenZi[i][c]=0;
						for(unsigned int j=0;j<nDiscreteCovs;j++){
							int Xij = currentParams.workDiscreteX(i,j);
							logPXiGivenZi[i][c]+=currentParams.workLogPhiStar(c,j,(unsigned int)Xij);
						}
						for(unsigned int j=0;j<nContinuousCovs;j++){
							xi(j)=currentParams.workContinuousX(i,j);
						}

						if (useIndependentNormal) {
							VectorXd muStar = currentParams.workMuStar(c);
							VectorXd sigma_cj(nContinuousCovs);
							for (unsigned int j = 0; j < nContinuousCovs; j++) {
								sigma_cj(j)= sqrt(1.0 / currentParams.Tau_Indep(c,j));
							}
							for (unsigned int j = 0; j<nContinuousCovs; j++) {
								logPXiGivenZi[i][c] += logPdfNormal(xi(j), muStar(j), sigma_cj(j));
							}
						}
						else {
							logPXiGivenZi[i][c] += logPdfMultivarNormal(nContinuousCovs, xi, currentParams.workMuStar(c), currentParams.workSqrtTau(c), currentParams.workLogDetTau(c));						
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
					for(unsigned int j=0;j<nDiscreteCovs;j++){
						if(!missingX[i][j]){
							int Xij = currentParams.workDiscreteX(i,j);
							logPXiGivenZi[i][c]+=currentParams.workLogPhiStar(c,j,(unsigned int)Xij);
						}
					}
				}
			}
		}

		// For the predictive subjects we do not count missing data
		LLT<MatrixXd> llt;
		for(unsigned int i=nSubjects;i<nSubjects+nPredictSubjects;i++){

			unsigned int nNotMissing=dataset.nContinuousCovariatesNotMissing(i);

			for(unsigned int c=0;c<maxNClusters;c++){
				if(u[i]<testBound[c]){
					VectorXd workMuStar=currentParams.workMuStar(c);
					VectorXd xi=VectorXd::Zero(nNotMissing);
					VectorXd muStar=VectorXd::Zero(nNotMissing);
					VectorXd sigma_cj = VectorXd::Zero(nNotMissing);
					MatrixXd sqrtTau = MatrixXd::Zero(nNotMissing, nNotMissing);
					double logDetTau = 0.0;
		
					if(nNotMissing==nContinuousCovs){
						muStar=workMuStar;
						if (useIndependentNormal) {
							for (unsigned int j = 0; j < nContinuousCovs; j++) {
								sigma_cj(j) = sqrt(1.0 / currentParams.Tau_Indep(c, j));
							}
						}
						else {
							sqrtTau = currentParams.workSqrtTau(c);
							logDetTau = currentParams.workLogDetTau(c);
						}
						for(unsigned int j=0;j<nContinuousCovs;j++){
							xi(j)=currentParams.workContinuousX(i,j);
						}
					}else{

						if (useIndependentNormal) {
							VectorXd workSigma = currentParams.Sigma_Indep(c);
							unsigned int j = 0;
							for (unsigned int j0 = 0; j0<nContinuousCovs; j0++) {
								if (!missingX[i][nDiscreteCovs + j0]) {
									xi(j) = currentParams.workContinuousX(i, j0);
									muStar(j) = workMuStar(j0);
									sigma_cj(j) = sqrt(workSigma(j0));
									j++;
								}
							}
						}
						else {
							MatrixXd workSigma = currentParams.Sigma(c);
							MatrixXd Sigma = MatrixXd::Zero(nNotMissing, nNotMissing);
							MatrixXd Tau = MatrixXd::Zero(nNotMissing, nNotMissing);
							unsigned int j = 0;
							for (unsigned int j0 = 0; j0<nContinuousCovs; j0++) {
								if (!missingX[i][nDiscreteCovs + j0]) {
									xi(j) = currentParams.workContinuousX(i, j0);
									muStar(j) = workMuStar(j0);
									unsigned int r = 0;
									for (unsigned int j1 = 0; j1<nContinuousCovs; j1++) {
										if (!missingX[i][nDiscreteCovs + j1]) {
											Sigma(j, r) = workSigma(j0, j1);
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

					}

					if (useIndependentNormal) {
						for (unsigned int j = 0; j<nNotMissing; j++) {
							logPXiGivenZi[i][c] += logPdfNormal(xi(j), muStar(j), sigma_cj(j));
						}
					}
					else {
						logPXiGivenZi[i][c] += logPdfMultivarNormal(nNotMissing, xi, muStar, sqrtTau, logDetTau);
					}

					
				}
			}
		}

	}
	double (*logPYiGivenZiWi)(const pReMiuMParams&, const pReMiuMData&,
											const unsigned int&,const int&,
											const unsigned int&)=NULL;
	if(includeResponse){
		if(!responseExtraVar){
			if(outcomeType.compare("Bernoulli")==0){
				logPYiGivenZiWi = &logPYiGivenZiWiBernoulli;
			}else if(outcomeType.compare("Binomial")==0){
				logPYiGivenZiWi = &logPYiGivenZiWiBinomial;
			}else if(outcomeType.compare("Poisson")==0){
				if (includeCAR){
					logPYiGivenZiWi = &logPYiGivenZiWiPoissonSpatial;
				}else{
					logPYiGivenZiWi = &logPYiGivenZiWiPoisson;
				}
			}else if(outcomeType.compare("Normal")==0){
				logPYiGivenZiWi = &logPYiGivenZiWiNormal;
			}else if(outcomeType.compare("Quantile")==0){
				logPYiGivenZiWi = &logPYiGivenZiWiQuantile;
			}else if(outcomeType.compare("Categorical")==0){
				logPYiGivenZiWi = &logPYiGivenZiWiCategorical;
			}else if(outcomeType.compare("Survival")==0){
				logPYiGivenZiWi = &logPYiGivenZiWiSurvival;
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
			if (predictType.compare("RaoBlackwell")==0){
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
		}
		if(includeResponse&&i>=nSubjects){
			if (predictType.compare("random")==0){
				// choose which component of the mixture we are sampling from
				double u=unifRand(rndGenerator);
				unsigned int c=0;
				while(cumPzGivenXy[c]<=u){
					c++;	
				}
				if(outcomeType.compare("Quantile")==0){
					// draw from the ALD distribution (Yu et al, 2005)
					// X1, X2 distributed Exp(1) then X1/p-X2/(1-p) has ALD(0,1;p)
					// if X distr ALD(0,1;p) then mu+sigma X has distribution ALD(mu,sigma;p)  
					double u1=unifRand(rndGenerator);
					double u2=unifRand(rndGenerator);
					double EXP1 = -log(u1);
					double EXP2 = -log(u2);
					double r= EXP1/hyperParams.pQuantile()-EXP2/(1-hyperParams.pQuantile());
					expectedTheta[0]=currentParams.sigmaSqY()*r+currentParams.theta(c,0);
					// line below should be a mistake because sigmaSqY is not a square for the ALD					
					//expectedTheta[0]=sqrt(currentParams.sigmaSqY())*r+currentParams.theta(c,0);
				}else{ 
					// draw from the normal distribution of that sample (only for yModel=Normal)
					// Create a normal random generator
					randomNormal normRand(0,1);
					expectedTheta[0]=sqrt(currentParams.sigmaSqY())*normRand(rndGenerator)+currentParams.theta(c,0);
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


		currentParams.z(i, zi, covariateType, useIndependentNormal);
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


void updateMissingPReMiuMData(baseGeneratorType& rndGenerator,
								pReMiuMParams& params,
								const pReMiuMOptions& options,
								pReMiuMData& dataset){

	unsigned int nSubjects = dataset.nSubjects();
	unsigned int nCovariates = dataset.nCovariates();
	unsigned int nDiscreteCovs = dataset.nDiscreteCovs();
	unsigned int nContinuousCovs = dataset.nContinuousCovs();
	vector<unsigned int> nCategories = params.nCategories();
	string covariateType = options.covariateType();

	// Define a uniform random number generator
	randomUniform unifRand(0,1);

	// For the missing variables we check which cluster the subject is allocated
	// to and then sample for the appropriate
	if(covariateType.compare("Discrete")==0){
		// We don't update the predictive subjects as their X values which
		// were missing are not used anywhere 
		// change: we now compute them because we are interested in looking at their posterior predictive distributions
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
			if(dataset.nContinuousCovariatesNotMissing(i)<nCovariates){
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

	}else if(covariateType.compare("Mixed")==0){
		// Discrete part of mixed type covariates
		// We don't update the predictive subjects as their X values which
		// were missing are not used anywhere
		for(unsigned int i=0;i<nSubjects;i++){
			int zi = params.z(i);
			for(unsigned int j=0;j<nDiscreteCovs;j++){
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

		// Normal part of mixed type covariates
		for(unsigned int i=0;i<nSubjects;i++){
			// Check if there is anything to do
			if(dataset.nContinuousCovariatesNotMissing(i)<nContinuousCovs){
				int zi = params.z(i);
				VectorXd newXi=multivarNormalRand(rndGenerator,params.workMuStar(zi),params.Sigma(zi));
				for(unsigned int j=0;j<nContinuousCovs;j++){
					if(dataset.missingX(i,nDiscreteCovs+j)){
						dataset.continuousX(i,j,newXi(j));
						params.workContinuousX(i,j,newXi(j));
					}else{
						newXi(j)=dataset.continuousX(i,j);
					}
				}
				double logVal = logPdfMultivarNormal(nContinuousCovs,newXi,params.workMuStar(zi),params.workSqrtTau(zi),params.workLogDetTau(zi));
				params.workLogPXiGivenZi(i,logVal);
			}
		}


	}


}


#endif /* DIPBACPROPOSALS_H_ */
