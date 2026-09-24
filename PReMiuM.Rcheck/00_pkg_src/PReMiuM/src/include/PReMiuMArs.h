/*** PReMiuMArs.h ***/
//
//    AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//    Modificado por: Paulino Perez Rodriguez
//    Fecha: 17/02/07
//    PURPOSE: Adaptive rejection sampling
//    Modified by Aurore J. Lavigne on 28/06/2014
//    Modified by Silvia Liverani on 24/07/2014
//
/* ********************************************************************************* */
#ifndef _ARSPReMiuM_H_
#define _ARSPReMiuM_H_


//standard include
#include <iostream>

//PReMiuM includes
#include<PReMiuMData.h>
#include<MCMC/model.h>
#include<Math/ars2.h>

//using namespace itpp;
using namespace std;


// Adaptive rejection sampler for spatial CAR model
double ARSsampleCAR(pReMiuMParams params,
                 const mcmcModel<pReMiuMParams,pReMiuMOptions,pReMiuMData>& model,
                 const unsigned int& iSub,
                 void (*evalhxhprimax)(const pReMiuMParams&,const mcmcModel<pReMiuMParams,pReMiuMOptions,pReMiuMData>&, const unsigned int&,const double&, double*, double*),
                 baseGeneratorType& rndGenerator)
{
	//initialise sampler
	const int ns=200; // ns=100 originally
	vector<double> xTmp;
	double ui=0, xlb, xub;
	int lb, ub;
	const int m = 5;
	xTmp.resize(m);
	ui=params.uCAR(iSub);

// Aurore initially had put
//	xTmp[0] = ui-50;
//        xTmp[1] = ui-2;
//        xTmp[2] = ui+0;
//        xTmp[3] = ui+2;
//        xTmp[4] = ui+50;
// with m=5 above

	xTmp[0] = ui-50;
        xTmp[1] = ui-2;
        xTmp[2] = ui+0;
        xTmp[3] = ui+2;
        xTmp[4] = ui+50;
	    xlb=0;
	    xub=0;
	    lb=0; //false
	    ub=0; //false
	double* x = &xTmp[0];
	double hx[m];
	double hpx[m];
	double y1=0;
	double y2=0;
    for (int i=0; i<m; i++){
        (*evalhxhprimax)(params, model, iSub, x[i], &y1, &y2);
        hx[i]=y1;
        hpx[i]=y2;
    }
    double emax=64;
    int iwv[ns+7];
    double rwv[6*(ns+1)+9];
    int ifault=0;

    initial_(&ns, &m, &emax, x, hx, hpx, &lb, &xlb, &ub, &xub, &ifault, iwv, rwv);

    //Check initialization is done properly
    if (ifault!=0){
	Rprintf("Error in the Adaptive Rejection Sampler");
	Rprintf("Error in ARS, cannot update uCAR for subject %d \n", iSub);
        Rprintf("Error in subroutine initial, ifault equals %d \n",ifault);
        return 0;
    }else{
        double beta=0;
        sample_( iwv, rwv, &beta,  &ifault, params, model, iSub, evalhxhprimax,rndGenerator );
        if (ifault!=0){
	Rprintf("Error in the Adaptive Rejection Sampler");
	Rprintf("Error in ARS, cannot update uCAR for subject %d \n", iSub);
            Rprintf("Error in subroutine sample, ifault equals %d \n",ifault);
            return 0;
        }else{
        return beta;
        }
    }//end of no error u=in subroutine initial
}

double ARSsampleNu(pReMiuMParams params,
                 const mcmcModel<pReMiuMParams,pReMiuMOptions,pReMiuMData>& model,
                 const unsigned int& cluster,
                 void (*evalhxhprimax)(const pReMiuMParams&,const mcmcModel<pReMiuMParams,pReMiuMOptions,pReMiuMData>&, const unsigned int&,const double&, double*, double*),
                 baseGeneratorType& rndGenerator)
{
	//initialise sampler
	const int ns=30; // number of attempts (max_attempts = ns*3)
	vector<double> xTmp;
	//double ui=0;
	double xlb, xub;
	int lb, ub;
	const int m = 7;
	xTmp.resize(m);
	//ui=params.nu(cluster);	
	xlb=0;
	xub=0;
	lb=1; //true
	ub=0; //false
	xTmp[0] = 0.0001;
        xTmp[1] = 0.05;
        xTmp[2] = 0.5;
        xTmp[3] = 1;
        xTmp[4] = 2;
        xTmp[5] = 5;
        xTmp[6] = 10;
//std::cout<<xTmp[0]<<" "<<xTmp[1]<<" "<<xTmp[2]<<" "<<xTmp[3]<<" "<<xTmp[4]<<" "<<std::endl;
// can try to remove some of these points on the x axis to improve efficiency
	double* x = &xTmp[0];
	double hx[m];
	double hpx[m];
	double y1=0;
	double y2=0;
    for (int i=0; i<m; i++){
        (*evalhxhprimax)(params, model, cluster, x[i], &y1, &y2);
        hx[i]=y1;
        hpx[i]=y2;

    }
// summary:
//uncomment below, the evaluations of the function; when it is not concave there is error ifault=4
// otherwise it works at the moment

// need to have small fixed effects (normalised) so that beta does not get too big

//std::cout<<"y "<<hx[0]<<" "<<hx[1]<<" "<<hx[2]<<" "<<hx[3]<<" "<<hx[4]<<" "<<std::endl;

    double emax=64;
    int iwv[ns+7];
    double rwv[6*(ns+1)+9];
    int ifault=0;

    initial_(&ns, &m, &emax, x, hx, hpx, &lb, &xlb, &ub, &xub, &ifault, iwv, rwv);

    //Check initialization is done properly
    if (ifault!=0){
	Rprintf("Error in the Adaptive Rejection Sampler");
	Rprintf("Error in ARS, cannot update nu (survival response)");
        Rprintf("Error in subroutine initial, ifault equals %d \n",ifault);
        return 0;
    }else{
        double beta=0;
        sample_( iwv, rwv, &beta,  &ifault, params, model, cluster, evalhxhprimax,rndGenerator );
        if (ifault!=0){
	Rprintf("Error in the Adaptive Rejection Sampler");
		Rprintf("Error in ARS, cannot update nu (survival response)");
            Rprintf("Error in subroutine sample, ifault equals %d \n",ifault);
            return 0;
        }else{
        return beta;
        }
    }//end of no error u=in subroutine initial
}


#endif
