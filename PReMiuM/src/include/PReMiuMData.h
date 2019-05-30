/// \file PReMiuMData.h
/// \author David Hastie
/// \brief Header file for data specification for PReMiuMpp

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

#ifndef DIPBACDATA_H_
#define DIPBACDATA_H_

// Standard includes
#include<cmath>
#include<vector>
#include<iostream>
#include<fstream>
#include<string>
#include<stdexcept>

using std::vector;
using std::ifstream;
using std::string;

/// \class pReMiuMData PReMiuMData.h "PReMiuMModel.h"
/// \brief A class for PReMiuMpp Data
class pReMiuMData{

	public:
		/// \brief Default constructor
		pReMiuMData(): _nSubjects(0), _nCovariates(0), _nFixedEffects(0), _nCategoriesY(0), _nPredictSubjects(0) {};

		/// \brief Default destructor
		~pReMiuMData(){};

		/// \brief Return the number of subjects
		unsigned int nSubjects() const{
			return _nSubjects;
		}

		/// \brief Return the number of subjects
		unsigned int& nSubjects(){
			return _nSubjects;
		}

		/// \brief Set the number of subjects
		void nSubjects(const unsigned int& nSubj){
			_nSubjects = nSubj;
		}

		unsigned int size() const{
			return _nSubjects;
		}

		/// \brief Return the number of covariates
		unsigned int nCovariates() const{
			return _nCovariates;
		}

		/// \brief Return the number of covariates
		unsigned int& nCovariates(){
			return _nCovariates;
		}

		/// \brief Set the number of covariates
		void nCovariates(const unsigned int& nCov){
			_nCovariates = nCov;
		}

		/// \brief Return the number of discrete covariates
		unsigned int nDiscreteCovs() const{
			return _nDiscreteCovs;
		}

		/// \brief Return the number of discrete covariates
		unsigned int& nDiscreteCovs(){
			return _nDiscreteCovs;
		}

		/// \brief Set the number of discrete covariates
		void nDiscreteCovs(const unsigned int& nDiscrCovs){
			_nDiscreteCovs = nDiscrCovs;
		}
		/// \brief Return the number of continuous covariates
		unsigned int nContinuousCovs() const{
			return _nContinuousCovs;
		}

		/// \brief Return the number of continuous covariates
		unsigned int& nContinuousCovs(){
			return _nContinuousCovs;
		}


		/// \brief Set the number of continuous covariates
		void nContinuousCovs(const unsigned int& nContCovs){
			_nContinuousCovs = nContCovs;
		}


		/// \brief Return the number of fixed effectss
		unsigned int nFixedEffects() const{
			return _nFixedEffects;
		}

		/// \brief Return the number of fixed effectss
		unsigned int& nFixedEffects(){
			return _nFixedEffects;
		}

		/// \brief Set the number of fixed effectss
		void nFixedEffects(const unsigned int& nConf){
			_nFixedEffects = nConf;
		}

		/// \brief Return the number of categories
		unsigned int nCategoriesY() const{
			return _nCategoriesY;
		}

		/// \brief Return the number of categories
		unsigned int& nCategoriesY(){
			return _nCategoriesY;
		}

		/// \brief Set the number of categories
		void nCategoriesY(const unsigned int& nCat){
			_nCategoriesY = nCat;
		}

		/// \brief Return the number of subjects
		unsigned int nPredictSubjects() const{
			return _nPredictSubjects;
		}

		/// \brief Return the number of subjects
		unsigned int& nPredictSubjects(){
			return _nPredictSubjects;
		}

		/// \brief Set the number of subjects
		void nPredictSubjects(const unsigned int& nPredSubj){
			_nPredictSubjects = nPredSubj;
		}



		/// \brief Return the vector of the number of categories
		vector<unsigned int> nCategories() const{
			return _nCategories;
		}

		/// \brief Return the vector of the number of categories
		vector<unsigned int>& nCategories(){
			return _nCategories;
		}


		/// \brief Set the vector of the number of categories
		void nCategories(const vector<unsigned int>& nCats){
			_nCategories.clear();
			_nCategories.resize(nCats.size());
			_nCategories.insert(_nCategories.begin(),nCats.begin(),nCats.end());
		}

		/// \brief Return the number of categories for covariate j
		unsigned int nCategories(const unsigned int& j) const{
			if(j>_nCovariates){
				throw std::range_error("nCategories subscript j out of range");
			}
			return _nCategories[j];
		}

		/// \brief Return the vector of the covariate names
		vector<string> covariateNames() const{
			return _covariateNames;
		}

		/// \brief Return the vector of the covariate names
		vector<string>& covariateNames(){
			return _covariateNames;
		}


		/// \brief Set the vector of the covariate names
		void covariateNames(const vector<string>& covNames){
			_covariateNames.clear();
			_covariateNames.resize(covNames.size());
			_covariateNames.insert(_covariateNames.begin(),covNames.begin(),covNames.end());
		}

		/// \brief Return name for covariate j
		string covariateNames(const unsigned int& j) const{
			return _covariateNames[j];
		}

		/// \brief Return the vector of the fixed effects names
		vector<string> fixedEffectNames() const{
			return _fixedEffectNames;
		}

		/// \brief Return the vector of the fixed effects names
		vector<string>& fixedEffectNames(){
			return _fixedEffectNames;
		}

		/// \brief Set the vector of the fixed effects names
		void fixedEffectNames(const vector<string>& fixEffNames){
			_fixedEffectNames.clear();
			_fixedEffectNames.resize(fixEffNames.size());
			_fixedEffectNames.insert(_fixedEffectNames.begin(),fixEffNames.begin(),fixEffNames.end());
		}

		/// \brief Return name for fixed effects j
		string fixedEffectNames(const unsigned int& j) const{
			return _fixedEffectNames[j];
		}


		/// \brief Return the outcome model type
		const string& outcomeType() const{
			return _outcomeType;
		}

		/// \brief Set the outcome model type
		void outcomeType(const string& outType){
			_outcomeType=outType;
		}

		/// \brief Return the outcome model type
		const string& covariateType() const{
			return _covariateType;
		}

		/// \brief Set the outcome model type
		void covariateType(const string& covType){
			_covariateType=covType;
		}

		/// \brief Return the output vector
		const vector<unsigned int>& discreteY() const{
			return _discreteY;
		}

		/// \brief Return the output vector
		vector<unsigned int>& discreteY() {
			return _discreteY;
		}

		/// \brief Set the output vector
		void discreteY(const vector<unsigned int>& yVec){
			_discreteY.clear();
			_discreteY.resize(yVec.size());
			_discreteY.insert(_discreteY.begin(),yVec.begin(),yVec.end());
		}

		/// \brief Return the output value for the ith subject
		unsigned int discreteY(const unsigned int& i) const{
			if(i>_nSubjects){
				throw std::range_error("y subscript i out of range");
			}
			return _discreteY[i];
		}

		/// \brief Return the output vector
		const vector<double>& continuousY() const{
			return _continuousY;
		}

		/// \brief Return the output vector
		vector<double>& continuousY() {
			return _continuousY;
		}

		/// \brief Set the output vector
		void continuousY(const vector<double>& yVec){
			_continuousY.clear();
			_continuousY.resize(yVec.size());
			_continuousY.insert(_continuousY.begin(),yVec.begin(),yVec.end());
		}

		/// \brief Return the output value for the ith subject
		double continuousY(const unsigned int& i) const{
			if(i>_nSubjects){
				throw std::range_error("y subscript i out of range");
			}
			return _continuousY[i];
		}

		/// \brief Return the covariate matrix
		const vector<vector<int> >& discreteX() const{
			return _discreteX;
		}

		/// \brief Return the covariate matrix
		vector<vector<int> >& discreteX(){
			return _discreteX;
		}

		/// \brief Return the jth covariate for subject i
		int discreteX(const unsigned int& i,const unsigned int& j) const{
			return _discreteX[i][j];
		}

		/// \brief Set the jth covariate for subject i
		void discreteX(const unsigned int& i,const unsigned int& j,const int& x){
			_discreteX[i][j]=x;
		}

		/// \brief Return the covariate matrix
		const vector<vector<double> >& continuousX() const{
			return _continuousX;
		}

		/// \brief Return the covariate matrix
		vector<vector<double> >& continuousX(){
			return _continuousX;
		}

		/// \brief Return the jth covariate for subject i
		double continuousX(const unsigned int& i,const unsigned int& j) const{
			return _continuousX[i][j];
		}

		/// \brief Set the jth covariate for subject i
		void continuousX(const unsigned int& i,const unsigned int& j,const double& x){
			_continuousX[i][j]=x;
		}

		/// \brief Return the missing covariate matrix
		const vector<vector<bool> >& missingX() const{
			return _missingX;
		}

		/// \brief Return the missing covariate matrix
		vector<vector<bool> >& missingX(){
			return _missingX;
		}

		/// \brief Return the jth covariate for subject i
		bool missingX(const unsigned int& i,const unsigned int& j) const{
			return _missingX[i][j];
		}

		/// \brief Return the number of covariates not missing for subject i
		unsigned int nContinuousCovariatesNotMissing(const unsigned int& i) const{
			return _nContinuousCovariatesNotMissing[i];
		}

		/// \brief Return the number of covariates not missing for each subject
		vector<unsigned int>& nContinuousCovariatesNotMissing(){
			return _nContinuousCovariatesNotMissing;
		}

		/// \brief Return the fixed effects matrix
		const vector<vector<double> >& W() const{
			return _W;
		}

		/// \brief Return the fixed effects matrix
		vector<vector<double> >& W(){
			return _W;
		}

		/// \brief Return the jth covariate for subject i
		double W(const unsigned int& i,const unsigned int& j) const{
			return _W[i][j];
		}

		/// \brief Return the logOffset vector
		const vector<double>& logOffset() const{
			return _logOffset;
		}

		/// \brief Return the logOffset vector
		vector<double>& logOffset(){
			return _logOffset;
		}

		/// \brief Return the logOffset for subject i
		double logOffset(const unsigned int& i) const{
			return _logOffset[i];
		}

		/// \brief Return the n vector
		const vector<unsigned int>& nTrials() const{
			return _nTrials;
		}

		/// \brief Return the n vector
		vector<unsigned int>& nTrials(){
			return _nTrials;
		}

		/// \brief Return n for subject i
		unsigned int nTrials(const unsigned int& i) const{
			return _nTrials[i];
		}


		/// \brief Return the n vector for Survival data
		const vector<unsigned int>& censoring() const{
			return _censoring;
		}

		/// \brief Return the n vector for Survival data 
		vector<unsigned int>& censoring(){
			return _censoring;
		}

		/// \brief Return n for subject i for Survival data
		unsigned int censoring(const unsigned int& i) const{
			return _censoring[i];
		}
		/// \brief Return the vector of lists of neighbours
		const vector<vector<unsigned int> >& neighbours() const{
			return _neighbours;
		}

		/// \brief Return the vector of lists of neighbours
		vector<vector<unsigned int> >& neighbours() {
			return _neighbours;
		}

		/// \brief Return the list of neighbours for subject i
		vector<unsigned int> neighbours(const unsigned int& i) const{
			return _neighbours[i];
		}

		/// \brief Return the j-th neighbour for subject i
		unsigned int neighbours(const unsigned int& i, const unsigned& j) const{
			return _neighbours[i][j];
		}

		/// \brief Set the vector of lists of neighbours
		void neighbours(const vector<vector<unsigned int> >&  neighvec){
			_neighbours = neighvec;
		}

		/// \brief Set the list of neighbours for subject i
		void neighbours(const vector<unsigned int>&  neighvec, const unsigned int& i){
			_neighbours[i] = neighvec;
		}

		/// \brief Return the number of neighbours for all subjects
		const vector<unsigned int>& nNeighbours() const{
			return _nNeighbours;
		}

		/// \brief Return the number of neighbours for all subjects
		vector<unsigned int>& nNeighbours(){
			return _nNeighbours;
		}

		/// \brief Return the number of neighbours for subject i
		unsigned int nNeighbours(const unsigned int& i) const{
			return _nNeighbours[i];
		}

		/// \brief Set the number of neighbours for all subjects
		void nNeighbours(const vector<unsigned int>& nNeigh){
			_nNeighbours=nNeigh;
		}

		/// \brief Set the number of neighbours for subject i
		void nNeighbours(const unsigned int& nNeigh, const unsigned int& i ){
			_nNeighbours[i]=nNeigh;
		}

		/// \brief Set includeCAR
		void includeCAR(const bool& incl ){
			_includeCAR=incl;
		}

		/// \brief return includeCAR
		bool& includeCAR(){
			return _includeCAR;
		}

		/// \brief return includeCAR
		bool includeCAR() const{
			return _includeCAR;
		}

	private:
		/// \brief The number of subjects
		unsigned int _nSubjects;

		/// \brief The number of covariates
		unsigned int _nCovariates;

		/// \brief The number of discrete covariates
		unsigned int _nDiscreteCovs;

		/// \brief The number of continuous covariates
		unsigned int _nContinuousCovs;

		/// \brief The number of fixed effects covariates
		unsigned int _nFixedEffects;

		/// \brief The number of categories for outcome discreteY when outcome is categorical
		unsigned int _nCategoriesY;

		/// \brief The number of subjects we are making predictions for
		unsigned int _nPredictSubjects;

		/// \brief A vector of the number of categories for each covariate
		vector<unsigned int> _nCategories;

		/// \brief A string describing the model for y
		string _outcomeType;

		/// \brief A string describing the model for X
		string _covariateType;


		/// \brief A vector of the output variables
		vector<unsigned int> _discreteY;

		/// \brief A vector of the output variables
		vector<double> _continuousY;


		/// \brief A matrix (vector of vectors) of the covariate data
		/// \note this is a signed int because missing values are typically stored
		/// as negative values
		vector<vector<int> > _discreteX;

		/// \brief A matrix (vector of vectors) of the covariate data
		/// \note this is a signed int because missing values are typically stored
		/// as negative values
		vector<vector<double> > _continuousX;

		/// \brief A vector of covariate names
		vector<string> _covariateNames;

		/// \brief A matrix (vector of vectors) of where there are missing
		/// covariate values
		vector<vector<bool> > _missingX;

		/// \brief A matrix of the number of non missing covariates for each subject
		vector<unsigned int> _nContinuousCovariatesNotMissing;

		/// \brief A matrix of the fixed effects covariates
		/// \note This may need to changed to be signed or double
		vector<vector<double> > _W;

		/// \brief A vector of fixed effects names
		vector<string> _fixedEffectNames;

		/// \brief A vector of logOffsets (only used in the Poisson model)
		vector<double> _logOffset;

		/// \brief A vector of n for each individual (only used in the Binomial model)
		vector<unsigned int> _nTrials;

		/// \brief A vector of n for each individual (only used in the survival model)
		vector<unsigned int> _censoring;

		/// \brief A vector of vector of neighbours for each subject
		vector<vector<unsigned int> > _neighbours;

		/// \brief A containing the number of neighbours for each subject
		vector<unsigned int> _nNeighbours;

		/// \brief Is the CAR term is included
		bool _includeCAR;

		/// \brief Is the CAR term is included
		bool _includeuCARinit;
};


#endif //DIPBACDATA_H_
