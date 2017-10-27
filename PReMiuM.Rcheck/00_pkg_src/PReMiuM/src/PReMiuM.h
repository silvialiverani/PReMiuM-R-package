#include <Rcpp.h>
#include <Rinternals.h> // for SEXP

// Note that the functions Gamma and LogGamma are mutually dependent.
double LogGamma(double);
double Gamma(double);

RcppExport SEXP profRegr(SEXP inputString);

RcppExport SEXP calcDisSimMat(SEXP fileName, SEXP nSweeps, SEXP nBurn, SEXP nFilter,SEXP nSubjects,SEXP nPredictSubjects, SEXP onlyLS);

RcppExport SEXP pYGivenZW(SEXP betaIn,SEXP thetaIn,SEXP zAlloc,SEXP sigmaBeta,
		SEXP sigmaTheta, SEXP dofTheta, SEXP dofBeta, SEXP nSubjects,
		SEXP yMat,SEXP betaW,SEXP nFixedEffects,SEXP nNames,SEXP constants,
		SEXP maxnNames);

RcppExport SEXP GradpYGivenZW(SEXP betaIn,SEXP thetaIn,SEXP zAlloc, SEXP nSubjects,
		SEXP betaW,SEXP yMat, SEXP nFixedEffects,SEXP nNames,SEXP maxnNames);

RcppExport SEXP pZpX(SEXP nClusters, SEXP nCategories,SEXP aPhi, SEXP n,SEXP nCovariates,SEXP zAlloc, SEXP xMat, SEXP nNames,SEXP alpha);


