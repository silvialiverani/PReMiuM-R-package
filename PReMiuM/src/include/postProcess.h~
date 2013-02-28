#ifndef _PReMiuM_POST_PROCESS
#define _PReMiuM_POST_PROCESS

#include <Rcpp.h>

// Note that the functions Gamma and LogGamma are mutually dependent.
double LogGamma(double);
double Gamma(double);

//RcppExport SEXP profRegr(SEXP input, SEXP output, SEXP hyper, SEXP predict, SEXP nSweeps, SEXP nBurn, SEXP nProgress, SEXP nFilter, SEXP nClusInit, 	SEXP seed, SEXP yModel, SEXP xModel, SEXP sampler, SEXP alpha, SEXP excludeY, SEXP extraYVar, SEXP varSelect, SEXP entropy);
RcppExport SEXP profRegr(SEXP inputString);

RcppExport SEXP calcDisSimMat(SEXP fileName, SEXP nSweeps, SEXP nBurn, SEXP nFilter,SEXP nSubjects,SEXP nPredictSubjects);

RcppExport SEXP pYGivenZW(SEXP betaIn,SEXP thetaIn,SEXP zAlloc,SEXP sigmaBeta,
		SEXP sigmaTheta, SEXP dofTheta, SEXP dofBeta, SEXP nSubjects,
		SEXP yMat,SEXP betaW,SEXP nFixedEffects,SEXP nNames,SEXP constants,
		SEXP maxnNames);

RcppExport SEXP GradpYGivenZW(SEXP betaIn,SEXP thetaIn,SEXP zAlloc, SEXP nSubjects,
		SEXP betaW,SEXP yMat, SEXP nFixedEffects,SEXP nNames,SEXP maxnNames);

RcppExport SEXP pZpX(SEXP nClusters, SEXP nCategories,SEXP aPhi, SEXP n,SEXP nCovariates,SEXP zAlloc, SEXP xMat, SEXP nNames,SEXP alpha);

#endif
