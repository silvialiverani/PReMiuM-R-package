#include <vector>
#include <ctime>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <stdexcept>

#include "include/postProcess.h"

using std::string;
using std::ifstream;
using std::vector;

SEXP calcDisSimMat(SEXP fileName, SEXP nSweeps, SEXP nBurn, SEXP nFilter,SEXP nSubjects,SEXP nPredictSubjects){
    
    string fName = Rcpp::as<string>(fileName);

    int nS = Rcpp::as<int>(nSweeps);
    int nB = Rcpp::as<int>(nBurn);
    int nF = Rcpp::as<int>(nFilter);
    unsigned long int nSj = Rcpp::as<int>(nSubjects);
    int nPSj = Rcpp::as<int>(nPredictSubjects);

    // Calculate how many samples we will have
    int nLines = 1 + (nS+nB)/nF;
    // Calculate when the burn in ends
    int firstLine = 2+nB/nF;

    vector<double> disSimMat((nSj*(nSj-1))/2+nPSj*nSj,1.0);

    // Open the file with sample in
    ifstream zFile;
    zFile.open(fName.c_str());

    // The denominator for the contributions
    double denom = 1+(double)nLines-(double)firstLine;

    // A temporary vector for storing the cluster from each sweep
    vector<int> clusterData(nSj+nPSj);
    for(int k=1;k<=nLines;k++){
    	if(k<firstLine){
        	// Ignore the burn in
    		for(int i=0;i<nSj+nPSj;i++){
    			int tmp=0;
    			zFile >> tmp;
    		}
     	}else{
    		if((1+k-firstLine)==1||(1+k-firstLine)%1000==0){
				std::cout << "Stage 1:" << 1+k-firstLine << " samples out of " << 1+nLines-firstLine << std::endl;
			}
			for(int i=0;i<nSj+nPSj;i++){
    			// Fill up the cluster data for this sweep
    			zFile >> clusterData[i];
    		}


    		// Now we need to populate the dissimilarity matrix
    		int r=0;
			for(int i=0;i<nSj-1;i++){
    			for(int ii=i+1;ii<nSj;ii++){
    				if(clusterData[i]==clusterData[ii]){
        				disSimMat[r]-=1.0/denom;
        			}
    				r++;
    			}
    		}

			for(int i=0;i<nPSj;i++){
				for(int ii=0;ii<nSj;ii++){
					if(clusterData[nSj+i]==clusterData[ii]){
						disSimMat[r]-=1.0/denom;
					}
					r++;
				}
			}
    	}
    }

    // Cheaper (for large datasets) to re-read everything in again, instead
    // of storing things
    zFile.seekg(0,std::ios::beg);
    int minIndex=0;
    double currMinSum=nSj*(nSj-1)/2.0;
    for(int k=1;k<=nLines;k++){
    	if(k<firstLine){
    		// Ignore the burn in
    		for(int i=0;i<nSj+nPSj;i++){
    			int tmp=0;
    			zFile >> tmp;
        	}
        }else{
        	if((1+k-firstLine)==1||(1+k-firstLine)%1000==0){
        		std::cout << "Stage 2:" << 1+k-firstLine << " samples out of " << 1+nLines-firstLine << std::endl;
    		}
        	for(int i=0;i<nSj+nPSj;i++){
        		// Fill up the cluster data for this sweep
        		zFile >> clusterData[i];
        	}

        	double tmpSum=0.0;
    		int r=0;
			for(int i=0;i<nSj-1;i++){
    			for(int ii=i+1;ii<nSj;ii++){
    				if(clusterData[i]==clusterData[ii]){
        				tmpSum+=pow(disSimMat[r],2.0);
        			}else{
        				tmpSum+=pow(1.0-disSimMat[r],2.0);
        			}
    				r++;
    			}
    		}
			if(tmpSum<currMinSum){
				minIndex=k;
				currMinSum=tmpSum;
			}
        }
    }

    zFile.close();

    return Rcpp::List::create(Rcpp::Named("lsOptSweep")=minIndex,
    						  Rcpp::Named("disSimMat")=Rcpp::wrap<vector<double> >(disSimMat));
}

SEXP pYGivenZW(SEXP betaIn,SEXP thetaIn,SEXP zAlloc,SEXP sigmaBeta,
		SEXP sigmaTheta, SEXP dofTheta, SEXP dofBeta, SEXP nSubjects,
		SEXP yMat,SEXP betaW,SEXP nFixedEffects,SEXP nNames,SEXP constants,
		SEXP maxnNames){

	Rcpp::NumericVector beta(betaIn);
	Rcpp::NumericVector bW(betaW);
    Rcpp::NumericVector theta(thetaIn);
    Rcpp::IntegerVector z(zAlloc);
    unsigned long int nSj = Rcpp::as<int>(nSubjects);
    Rcpp::NumericVector yM(yMat);
    int nFE = Rcpp::as<int>(nFixedEffects);
    double con = Rcpp::as<double>(constants);
    Rcpp::IntegerVector nN(nNames);
    int maxnN = Rcpp::as<int>(maxnNames);
    double sigmaB = Rcpp::as<double>(sigmaBeta);
    double sigmaT = Rcpp::as<double>(sigmaTheta);
    double dofB = Rcpp::as<double>(dofBeta);
    double dofT = Rcpp::as<double>(dofTheta);

    vector<double> thetaTmp(maxnN+1);
    unsigned int k=0;
	for (unsigned int i=0; i<nN.size();i++) {
		thetaTmp.at(nN[i]) = theta[k];
		k++;
	}
	vector<double> pred(nSj);
	for (unsigned int i=0; i<nSj;i++) {
		pred[i] = thetaTmp.at(z(i));
	}
	vector<double> predictorAll(nSj);
	for (unsigned int i=0; i<nSj;i++) {
		predictorAll[i] = pred[i] + bW[i];
	}
	double output = 0.0;
	for (unsigned int i=0;i<nSj;i++){
		output = output + yM[i] * predictorAll.at(i)-log(1+exp(predictorAll.at(i)));
	}
	// contribution from theta
	double paramsThetaTmp = 1/pow(sigmaT,2);
	double outputTheta = 0.0;
	for (unsigned int i=0;i<theta.size();i++){
		outputTheta = outputTheta+ log(dofT+pow(theta[i],2)*paramsThetaTmp);
	}
	outputTheta = outputTheta*0.5*(dofT+1);
	// contribution from beta
	double paramsBetaTmp = 1/pow(sigmaB,2);
	double outputBeta = 0.0;
	for (unsigned int i=0;i<beta.size();i++){
		outputBeta = outputBeta+ log(dofB+pow(beta[i],2)*paramsBetaTmp);
	}
	outputBeta = outputBeta*0.5*(dofB+1);
	//can i drop constants?
	double out = -(output-outputTheta-outputBeta+con)/nSj;
    return Rcpp::List::create(Rcpp::Named("pYGivenZW")=out);
}


SEXP GradpYGivenZW(SEXP betaIn,SEXP thetaIn,SEXP zAlloc, SEXP nSubjects,
		SEXP betaW,SEXP yMat, SEXP nFixedEffects,SEXP nNames,SEXP maxnNames){
	Rcpp::NumericVector beta(betaIn);
	Rcpp::NumericVector bW(betaW);
    Rcpp::NumericVector theta(thetaIn);
    Rcpp::IntegerVector z(zAlloc);
    unsigned long int nSj = Rcpp::as<int>(nSubjects);
    int nFE = Rcpp::as<int>(nFixedEffects);
    Rcpp::IntegerVector nN(nNames);
    int maxnN = Rcpp::as<int>(maxnNames);
    Rcpp::NumericVector yM(yMat);

    vector<double> thetaTmp(maxnN+1);
    unsigned int k=0;
	for (unsigned int i=0; i<nN.size();i++) {
		thetaTmp.at(nN[i]) = theta[k];
		k++;
	}
	vector<double> pred(nSj);
	for (unsigned int i=0; i<nSj;i++) {
		pred[i] = thetaTmp.at(z(i));
	}
	Rcpp::NumericVector yPred(nSj);
	for (unsigned int i=0;i<nSj;i++){
		yPred[i] = yM[i]-1/(1+1/exp(pred[i]+bW[i]));
	}

    return Rcpp::wrap(yPred);
}

SEXP pZpX(SEXP nClusters, SEXP nCategories,SEXP aPhi, SEXP n,SEXP nCovariates,SEXP zAlloc, SEXP xMat, SEXP nNames,SEXP alpha){

	int nC = Rcpp::as<int>(nClusters);
	Rcpp::NumericVector nCat(nCategories);
	Rcpp::NumericVector nn(n);
	double aP = Rcpp::as<double>(aPhi);
	int nCov = Rcpp::as<int>(nCovariates);
	Rcpp::NumericVector xM(xMat);
	Rcpp::NumericVector z(zAlloc);
	Rcpp::NumericVector names(nNames);
	double a = Rcpp::as<double>(alpha);

	unsigned int nSj = z.size();

	double out = 0.0;
	for (unsigned int k=0;k<nCov;k++){
		out = out+nC*(LogGamma(nCat[k]*aP)-nCat[k]*LogGamma(aP));
	}

	double sumNCat = std::accumulate(nCat.begin(),nCat.end(),0);
	vector<double> nCatK(sumNCat*nC,0);
	for (unsigned int c=0;c<nC;c++){
		unsigned int index = 0;
		for (unsigned int k=0;k<nCov;k++){
			out = out-LogGamma(nCat[k]*aP+nn[c]);
		} 
		for (unsigned int j=0;j<nCov;j++){
			for (unsigned int kj=0;kj<nCat[j];kj++){
				for (unsigned int zi=0;zi<nSj;zi++){
					if (z[zi] == names[c]){
						if (xM[j*nSj+zi]==kj){
							nCatK[c*sumNCat+index]++;
						}
					}			
				}
				index = index+1;
			}	
		}		
	}
	for (unsigned int i=0; i< nCatK.size();i++){
		out = out + LogGamma(aP + nCatK[i]);
	}

	out = out + nC * log(a) + LogGamma(nSj+1);
	vector<double> counts(nSj,0);
	for (unsigned int i=0;i<nC;i++){
		counts[nn[i]-1] = counts[nn[i]-1]+1;
	}
	for (unsigned int i=0;i<nSj;i++){
		out = out -log(a+i)-counts[i]*log(i+1)-LogGamma(counts[i]+1);
	}

	return Rcpp::wrap(out);
}


// Note that the functions Gamma and LogGamma are mutually dependent.

double Gamma
(
    double x    // We require x > 0
)
{
	if (x <= 0.0)
	{
		std::stringstream os;
        os << "Invalid input argument " << x <<  ". Argument must be positive.";
        throw std::invalid_argument( os.str() ); 
	}

    // Split the function domain into three intervals:
    // (0, 0.001), [0.001, 12), and (12, infinity)

    ///////////////////////////////////////////////////////////////////////////
    // First interval: (0, 0.001)
	//
	// For small x, 1/Gamma(x) has power series x + gamma x^2  - ...
	// So in this range, 1/Gamma(x) = x + gamma x^2 with error on the order of x^3.
	// The relative error over this interval is less than 6e-7.

	const double gamma = 0.577215664901532860606512090; // Euler's gamma constant

    if (x < 0.001)
        return 1.0/(x*(1.0 + gamma*x));

    ///////////////////////////////////////////////////////////////////////////
    // Second interval: [0.001, 12)
    
	if (x < 12.0)
    {
        // The algorithm directly approximates gamma over (1,2) and uses
        // reduction identities to reduce other arguments to this interval.
		
		double y = x;
        int n = 0;
        bool arg_was_less_than_one = (y < 1.0);

        // Add or subtract integers as necessary to bring y into (1,2)
        // Will correct for this below
        if (arg_was_less_than_one)
        {
            y += 1.0;
        }
        else
        {
            n = static_cast<int> (floor(y)) - 1;  // will use n later
            y -= n;
        }

        // numerator coefficients for approximation over the interval (1,2)
        static const double p[] =
        {
            -1.71618513886549492533811E+0,
             2.47656508055759199108314E+1,
            -3.79804256470945635097577E+2,
             6.29331155312818442661052E+2,
             8.66966202790413211295064E+2,
            -3.14512729688483675254357E+4,
            -3.61444134186911729807069E+4,
             6.64561438202405440627855E+4
        };

        // denominator coefficients for approximation over the interval (1,2)
        static const double q[] =
        {
            -3.08402300119738975254353E+1,
             3.15350626979604161529144E+2,
            -1.01515636749021914166146E+3,
            -3.10777167157231109440444E+3,
             2.25381184209801510330112E+4,
             4.75584627752788110767815E+3,
            -1.34659959864969306392456E+5,
            -1.15132259675553483497211E+5
        };

        double num = 0.0;
        double den = 1.0;
        int i;

        double z = y - 1;
        for (i = 0; i < 8; i++)
        {
            num = (num + p[i])*z;
            den = den*z + q[i];
        }
        double result = num/den + 1.0;

        // Apply correction if argument was not initially in (1,2)
        if (arg_was_less_than_one)
        {
            // Use identity gamma(z) = gamma(z+1)/z
            // The variable "result" now holds gamma of the original y + 1
            // Thus we use y-1 to get back the orginal y.
            result /= (y-1.0);
        }
        else
        {
            // Use the identity gamma(z+n) = z*(z+1)* ... *(z+n-1)*gamma(z)
            for (i = 0; i < n; i++)
                result *= y++;
        }

		return result;
    }

    ///////////////////////////////////////////////////////////////////////////
    // Third interval: [12, infinity)

    if (x > 171.624)
    {
		// Correct answer too large to display. Force +infinity.
		double temp = DBL_MAX;
		return temp*2.0;
    }

    return exp(LogGamma(x));
}

double LogGamma
(
    double x    // x must be positive
)
{
	if (x <= 0.0)
	{
		std::stringstream os;
        os << "Invalid input argument " << x <<  ". Argument must be positive.";
        throw std::invalid_argument( os.str() ); 
	}

    if (x < 12.0)
    {
        return log(fabs(Gamma(x)));
    }

	// Abramowitz and Stegun 6.1.41
    // Asymptotic series should be good to at least 11 or 12 figures
    // For error analysis, see Whittiker and Watson
    // A Course in Modern Analysis (1927), page 252

    static const double c[8] =
    {
		 1.0/12.0,
		-1.0/360.0,
		1.0/1260.0,
		-1.0/1680.0,
		1.0/1188.0,
		-691.0/360360.0,
		1.0/156.0,
		-3617.0/122400.0
    };
    double z = 1.0/(x*x);
    double sum = c[7];
    for (int i=6; i >= 0; i--)
    {
        sum *= z;
        sum += c[i];
    }
    double series = sum/x;

    static const double halfLogTwoPi = 0.91893853320467274178032973640562;
    double logGamma = (x - 0.5)*log(x) - x + halfLogTwoPi + series;    
	return logGamma;
}

