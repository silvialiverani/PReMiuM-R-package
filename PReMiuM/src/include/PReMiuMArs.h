/*** PReMiuMArs.h ***/
//
//    AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//    Modificado por: Paulino Perez Rodriguez
//    Fecha: 17/02/07
//    PURPOSE: Adaptive rejection sampling
//    Modified by Aurore J. Lavigne on 28/06/2014
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

double ARSsample(pReMiuMParams params,
                 const mcmcModel<pReMiuMParams,pReMiuMOptions,pReMiuMData>& model,
                 const unsigned int& iSub,
                 void (*evalhxhprimax)(const pReMiuMParams&,const mcmcModel<pReMiuMParams,pReMiuMOptions,pReMiuMData>&, const unsigned int&,const double&, double*, double*),
                 baseGeneratorType& rndGenerator)
{
    //initialise sampler
    const int ns=100;
    const int m=5;
    double x[m];//pointeur vers le tableau des valeurs intiales
    double hx[m];
    double hpx[m];
    double ui=params.uCAR(iSub);
        x[0] = ui-50;
        x[1] = ui-2;
        x[2] = ui+0;
        x[3] = ui+2;
        x[4] = ui+50;
    double y1=0;
    double y2=0;
    for (int i=0; i<m; i++){
        (*evalhxhprimax)(params, model, iSub, x[i], &y1, &y2);
        hx[i]=y1;
        hpx[i]=y2;
    }
    double emax=64;
    double xlb=0;
    double xub=0;
    int lb=0; //false
    int ub=0; //false
    int iwv[ns+7];
    double rwv[6*(ns+1)+9];
    int ifault=0;

    initial_(&ns, &m, &emax, x, hx, hpx, &lb, &xlb, &ub, &xub, &ifault, iwv, rwv);

    //Check initialization is done properly
    if (ifault!=0){
        Rprintf("Error in ARS, cannot update uCAR for subject %d \n", iSub);
        Rprintf("Error in subroutine initial, ifault equals %d \n",ifault);
        return 0;
    }else{
        double beta=0;
        sample_( iwv, rwv, &beta,  &ifault, params, model, iSub, evalhxhprimax,rndGenerator );
        if (ifault!=0){
            Rprintf("Error in ARS, cannot update uCAR for subject %d \n", iSub);
            Rprintf("Error in subroutine sample, ifault equals %d \n",ifault);
            return 0;
        }else{
        return beta;
        }
    }//end of no error u=in subroutine initial
}




#endif
