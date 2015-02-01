#include <iostream>
#include <stdio.h>
#include <vector>

#include "limits.h"
#include "math.h"
#include "likelihood.h"
#include "idealG.h"
#include "missed_eventsG.h"
#include "occupancies.h"
#include "idealG.h"
#include "errors.h"

extern "C" int add(int a, int b) {
    int c = a + b;
    return c;
}

extern "C" int dcpDetWs ( double * determinant , double * Q, int nopen, int k, double tres, double s, bool isOpen ) {

    /* create Q-matrix */
    DCProgs::t_rmatrix qmatrix(k ,k);
    for (int i=0; i < k; i++){
        for (int j=0; j < k;j++){
            qmatrix(i,j) = Q[(j*k)+i];
        }
    }

    if (tres < 0) {
        std::cerr << "tres cannot be negative" << std::endl;
        return 1;
    }

    try {
        DCProgs::DeterminantEq detWs (DCProgs::QMatrix(qmatrix , nopen), tres);
        if (isOpen)
            *determinant = detWs(s);
        else {
            DCProgs::DeterminantEq detWsT = detWs.transpose();
            *determinant = detWsT(s);
        }

    } catch (const DCProgs::errors::Root& e){
        std::cerr << "[WARN]: Error in dcprogs::dcpDetWs" << std::endl;
        std::cerr << e.what() << std::endl;
        return 1;
    }
    return 0;
}

extern "C" int dcpDerivDetWs ( double * deriv , double * Q, int nopen, int k, double tres, double s, bool isOpen ) {

    /* create Q-matrix */
    DCProgs::t_rmatrix qmatrix(k ,k);
    for (int i=0; i < k; i++){
        for (int j=0; j < k;j++){
            qmatrix(i,j) = Q[(j*k)+i];
        }
    }

    if (tres < 0) {
        std::cerr << "tres cannot be negative" << std::endl;
        return 1;
    }

    DCProgs::t_rmatrix detderiv;
    try {
        DCProgs::DeterminantEq detWs (DCProgs::QMatrix(qmatrix,nopen), tres);
        if (isOpen)
            detderiv = detWs.s_derivative(s);
        else {
            DCProgs::DeterminantEq detWsT = detWs.transpose();
            detderiv = detWsT.s_derivative(s);
        }

    } catch (const DCProgs::errors::Root& e){
        std::cerr << "[WARN]: Error in dcprogs::dcpDerivDetWs" << std::endl;
        std::cerr << e.what() << std::endl;
        return 1;
    } 

    int elements;
    isOpen ? (elements = nopen * nopen) : (elements = ( k - nopen ) * ( k - nopen ));

    for ( int index = 0; index < elements; index++ ) {
        deriv[index] = detderiv(index);
    }
    return 0;
}

extern "C" int dcpSurvivorXs ( double * survivor, double * Q, int nopen, int k, double tres, double s, bool isOpen){
    /* create Q-matrix */
    DCProgs::t_rmatrix qmatrix(k ,k);
    for (int i=0; i < k; i++){
        for (int j=0; j < k;j++){
            qmatrix(i,j) = Q[(j*k)+i];
        }
    }

    if (tres < 0) {
        std::cerr << "tres cannot be negative" << std::endl;
        return 1;
    }


    DCProgs::t_rmatrix pSurvive;

    try {
        DCProgs::LaplaceSurvivor survivor (DCProgs::QMatrix(qmatrix , nopen));
        if (isOpen)
            pSurvive = survivor.H(s, tres);
        else {
            DCProgs::LaplaceSurvivor survivorT = survivor.transpose();
            pSurvive = survivorT.H(s, tres);
        }


    } catch (const DCProgs::errors::Root& e){
        std::cerr << "[WARN]: Error in dcprogs::dcpSurvivorXs" << std::endl;
        std::cerr << e.what() << std::endl;
        return 1;
    }

    int elements;
    isOpen ? (elements = nopen * nopen) : (elements = ( k - nopen ) * ( k - nopen ));

    for ( int index = 0; index < elements; index++ ) {
        survivor[index] = pSurvive(index);
    }
    return 0;
}

extern "C" int dcpDerivSurvivorXs( double * survivor, double * Q, int nopen, int k, double tres, double s, bool isOpen ) {

    DCProgs::t_rmatrix qmatrix(k ,k);
    for (int i=0; i < k; i++){
        for (int j=0; j < k;j++){
            qmatrix(i,j) = Q[(j*k)+i];
        }
    }

    if (tres < 0) {
        std::cerr << "tres cannot be negative" << std::endl;
        return 1;
    }


    DCProgs::t_rmatrix pSurvive;

    try {
        DCProgs::LaplaceSurvivor survivor (DCProgs::QMatrix(qmatrix , nopen));
        if (isOpen)
            pSurvive = survivor.s_derivative(s, tres);
        else {
            DCProgs::LaplaceSurvivor survivorT = survivor.transpose();
            pSurvive = survivorT.s_derivative(s, tres);
        }

    } catch (const DCProgs::errors::Root& e){
        std::cerr << "[WARN]: Error in dcprogs::dcpSurvivorXs" << std::endl;
        std::cerr << e.what() << std::endl;
        return 1;
    }

    int elements;
    isOpen ? (elements = nopen * nopen) : (elements = ( k - nopen ) * ( k - nopen ));

    for ( int index = 0; index < elements; index++ ) {
        survivor[index] = pSurvive(index);
    }
    return 0;
}

extern "C" int dcpApproxSurvivorComponent( double * survivor, double * root, double * Q, int nopen, int k, double tres, int ith, bool isOpen, double * dcpOptions) {

    /*
        returns the i'th component of the AR matrices where the i'th iterate is the i'th open or closed state
    */

     /* create Q-matrix */
    DCProgs::t_rmatrix qmatrix(k ,k);
    for (int i=0; i < k; i++){
        for (int j=0; j < k;j++){
            qmatrix(i,j) = Q[(j*k)+i];
        }
    }

    if (tres < 0) {
        std::cerr << "tres cannot be negative" << std::endl;
        return 1;
    }

    if ( ith < 0 || (isOpen && ith >= nopen) || (!isOpen && ith >= (k - nopen))) {
        std::cerr << "ith cannot be negative or out of bounds for number of open or closed states" << std::cout;
        return 1;
    }


    DCProgs::Asymptotes::t_MatrixAndRoot component;
    DCProgs::t_rmatrix ARi;

    try { 
        // Create survivor object A_S - ignore first parameter for rootfinding
        double xtol = dcpOptions[1];
        double rtol = dcpOptions[2];
        int itermax = (int)dcpOptions[3];
        double lower_bound = dcpOptions[4];
        double upper_bound = dcpOptions[5];    
        
        DCProgs::ApproxSurvivor survivor(qmatrix , tres, xtol, rtol, itermax,  lower_bound, upper_bound);
        
        if (isOpen){
            ARi = std::get<0>(survivor.get_af_components(ith));
            *root = std::get<1>(survivor.get_af_components(ith));
        }
        else {
            ARi = std::get<0>(survivor.get_fa_components(ith));
            *root = std::get<1>(survivor.get_fa_components(ith));
        }

    } catch (const DCProgs::errors::Root& e) {
        std::cerr << "[WARN]: Error in dcprogs::dcpApproxSurvivorComponent" << std::endl;
        std::cerr << e.what() << std::endl;
        return 1;
    }

    //return the result for the AR(ith) matrix 

    int elements;
    isOpen ? (elements = nopen * nopen) : (elements = ( k - nopen ) * ( k - nopen ));

    for ( int index = 0; index < elements; index++ ) {
        survivor[index] = ARi(index);
    }
 
    return 0;

}

extern "C" int dcpFindRoots( double * roots, double * Q, int nopen, int k, double tres, bool isOpen, double * dcpOptions ) {

     /* create Q-matrix */
    DCProgs::t_rmatrix qmatrix(k ,k);
    for (int i=0; i < k; i++){
        for (int j=0; j < k; j++){
            qmatrix(i,j) = Q[(j*k)+i];
        }
    }

    if (tres < 0) {
        std::cerr << "tres cannot be negative" << std::endl;
        return 1;
    }

    std::vector<DCProgs::Root> asymptoticroots;

    try {
        DCProgs::DeterminantEq detWs(DCProgs::QMatrix(qmatrix , nopen), tres);

        double xtol = dcpOptions[1];
        double rtol = dcpOptions[2];
        int itermax = dcpOptions[3];
        double lower_bound = dcpOptions[4];
        double upper_bound = dcpOptions[5];

        if (isOpen)
            asymptoticroots = DCProgs::find_roots(detWs,xtol,rtol,itermax,lower_bound,upper_bound);
        else {
            DCProgs::DeterminantEq detWsT = detWs.transpose();
            asymptoticroots = DCProgs::find_roots(detWsT,xtol,rtol,itermax,lower_bound,upper_bound);
        }

    } catch (const DCProgs::errors::Root& e) {
        std::cerr << "[WARN]: error in dcpFindRoots" << std::endl;
        std::cerr << e.what() << std::endl;
        return 1;
    }

    int index = 0;
    for(DCProgs::Root const &root: asymptoticroots) {
        roots[index] = root.root;
        std::cout << root.root <<std::endl;
        index++;
    }
    return 0;
}

extern "C" int dcpExactSurvivorXt ( double * survivor, double * Q, int nopen, int k, double tres, double t, bool isOpen) {

    /* create Q-matrix */
    DCProgs::t_rmatrix qmatrix(k ,k);
    for (int i=0; i < k; i++){
        for (int j=0; j < k;j++){
            qmatrix(i,j) = Q[(j*k)+i];
        }
    }

    if (tres < 0) {
        std::cerr << "tres cannot be negative" << std::endl;
        return 1;
    }

    if (t < 0) {
        std::cerr << "t cannot be negative" << std::cout;
        return 1;
    }

    DCProgs::t_rmatrix pSurvive;

    try {
        DCProgs::ExactSurvivor survivor (DCProgs::QMatrix(qmatrix , nopen), tres);
        if (isOpen)
            pSurvive = survivor.af(t);
        else
            pSurvive = survivor.fa(t);


    } catch (const DCProgs::errors::Root& e){
        std::cerr << "[WARN]: Error in dcprogs::dcpSurvivorXt" << std::endl;
        std::cerr << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "[WARN]: Unknown Error in dcprogs::dcpSurvivorXt" << std::endl;
        return 1;
    }

    int elements;
    isOpen ? (elements = nopen * nopen) : (elements = ( k - nopen ) * ( k - nopen ));
    for ( int index = 0; index < elements; index++ ) {
        survivor[index] = pSurvive(index);
    }
    return 0;
}

extern "C" int dcpExactSurvivorRecursiveMatrix ( double * recursiveMatrix , double * Q, int nopen, int k, double tres, bool isOpen, int i, int m, int r) {

    /* create Q-matrix */
    DCProgs::t_rmatrix qmatrix(k ,k);
    for (int i=0; i < k; i++){
        for (int j=0; j < k;j++){
            qmatrix(i,j) = Q[(j*k)+i];
        }
    }

    if (tres < 0) {
        std::cerr << "tres cannot be negative" << std::cout;
        return 1;
    }

    if (i > k) {
        std::cerr << "i cannot be greater than the number of states in the Q-matrix" << std::cout;
        return 1;
    }

    if (m < 0) {
        std::cerr << "m cannot be negative" << std::cout;
        return 1;
    }

    if (r < 0) {
        std::cerr << "r cannot be negative" << std::cout;
        return 1;
    }

    DCProgs::t_rmatrix pSurvive;

    try {
        DCProgs::ExactSurvivor survivor (DCProgs::QMatrix(qmatrix , nopen), tres);
        if (isOpen)
            pSurvive = survivor.recursion_af(i, m, r);
        else
            pSurvive = survivor.recursion_fa(i, m, r);


    } catch (const DCProgs::errors::Root& e){
        std::cerr << "[WARN]: Error in dcprogs::dcpExactSurvivorRecursiveMatrix" << std::endl;
        std::cerr << e.what() << std::endl;
        return 1;
    }

    int elements;
    isOpen ? (elements = nopen * nopen) : (elements = ( k - nopen ) * ( k - nopen ));

    for ( int index = 0; index < elements; index++ ) {
        recursiveMatrix[index] = pSurvive(index);
    }
    return 0;
}


extern "C" int dcpCHSOccupancies ( double * phi, double * Q, int nopen, int k, double tres, double tcrit, bool initial, double * dcpOptions) {

    /* create Q-matrix */
    DCProgs::t_rmatrix qmatrix(k ,k);
    for (int i=0; i < k; i++){
        for (int j=0; j < k;j++){
            qmatrix(i,j) = Q[(j*k)+i];
        }
    }

    if (tcrit < 0) {
        std::cerr << "tcrit cannot be negative" << std::endl;
        return 1;
    }

    if (tres < 0) {
        std::cerr << "tres cannot be negative" << std::cout;
        return 1;
    }

    int nmax = dcpOptions[0];
    double xtol = dcpOptions[1];
    double rtol = dcpOptions[2];
    int itermax = dcpOptions[3];
    double lower_bound = dcpOptions[4];
    double upper_bound = dcpOptions[5];

    DCProgs::t_initvec chsvec;

    try {
        DCProgs::MissedEventsG missedeventsG(DCProgs::QMatrix(qmatrix,nopen),tres,nmax, xtol, rtol, itermax, lower_bound, upper_bound);
        if (initial)
            chsvec = DCProgs::CHS_occupancies(missedeventsG , tcrit);
        else
            chsvec = DCProgs::CHS_occupancies(missedeventsG , tcrit , false);
 
    } catch (const DCProgs::errors::Root& e) {
        std::cerr << "[WARN]: Error in scprogs::dcpCHSOccupancies" << std::endl;
        std::cerr << e.what() << std::endl;
        return 1;
    }

    int elements;
    initial ? (elements = nopen ) : (elements = ( k - nopen ) );

    for ( int index = 0; index < elements; index++ ) {
        phi[index] = chsvec(index);
    }
 
    return 0;
}

extern "C" int dcpOccupancies ( double * phi, double * Q, int nopen, int k, double tau, bool initial, double * dcpOptions ) {

    /* create Q-matrix */
    DCProgs::t_rmatrix qmatrix(k ,k);
    for (int i=0; i < k; i++){
        for (int j=0; j < k;j++){
            qmatrix(i,j) = Q[(j*k)+i];
        }
    }

    DCProgs::t_initvec occvec;

    int nmax = dcpOptions[0];
    double xtol = dcpOptions[1];
    double rtol = dcpOptions[2];
    int itermax = dcpOptions[3];
    double lower_bound = dcpOptions[4];
    double upper_bound = dcpOptions[5];

    try {
        DCProgs::MissedEventsG eG(DCProgs::QMatrix(qmatrix,nopen) , tau , nmax , xtol, rtol , itermax , lower_bound, upper_bound);
        occvec = DCProgs::occupancies(eG , initial );

    } catch (const DCProgs::errors::Root& e ) {
        std::cerr << "[WARN] : Error thrown in dcprogs::dcpOccupancies" << std::endl;
        std::cerr << e.what() << std::endl;
        return 1;
    }

    //assign results
    int elements;
    initial ? (elements = nopen ) : (elements = ( k - nopen ));

    for ( int index = 0; index < elements; index++ ) {
        phi[index] = occvec(index);
    }
 
    return 0;
}

extern "C" int dcpExactGXYt( double * GXYt, double * Q, int nopen, int k, double tres, double t, bool isAF, double * dcpOptions) {
    DCProgs::t_rmatrix qmatrix(k ,k);
    for (int i=0; i < k; i++){
        for (int j=0; j < k;j++){
            qmatrix(i,j) = Q[(j*k)+i];
        }
    }

    if (t < 0) {
        std::cerr << "t cannot be negative" << std::endl;
        return 1;
    }

    if (tres < 0) {
        std::cerr << "tres cannot be negative" << std::cout;
        return 1;
    }

    DCProgs::t_rmatrix expT;

    try {

        int nmax = dcpOptions[0];
        double xtol = dcpOptions[1];
        double rtol = dcpOptions[2];
        int itermax = dcpOptions[3];
        double lower_bound = dcpOptions[4];
        double upper_bound = dcpOptions[5];

        DCProgs::MissedEventsG missedeventsG(DCProgs::QMatrix(qmatrix,nopen) , tres, nmax, xtol, rtol, itermax, lower_bound, upper_bound);

        if (isAF) {
            expT = missedeventsG.af(t);
        } else {
            expT = missedeventsG.fa(t);
        }

    } catch (const DCProgs::errors::Root& e) {
        std::cerr << "[WARN]: Error thrown in DCProgs::MissedEvents" << std::endl;
        std::cerr << e.what() << std::endl;
        return 1;
    }

    //return the call
    int index;
    for (index = 0; index < (nopen * (k - nopen)); index++ ){
        GXYt[index] = expT(index);
    }

    return 0;
}

extern "C" int dcpIdealGXYt( double * GXYt, double * Q, int nopen, int k, double t , bool isAF  ) {
     
    /* error handling */
     
    if (false){
        std::stringstream ss;
        ss << "calculates ideal X to Y transition exp(Q_XX * t) * Q_XY" << std::endl;
        ss << "Four arguments expected - Q-matrix ,nopen, t. isAF" << std::endl;
        ss << "Q-matrix - k x k matrix double - the Q matrix" << std::endl;
        ss << "nopen - scalar int - number of open states in mechanism" << std::endl;
        ss << "t - scalar double - duration time" << std::endl;
        ss << "isAF - scalar bool - calculate exp(Q_AF * t) or exp(Q_FA * t)" << std::endl;
        std::cout << ss.str().c_str() << std::endl;
    }   
     
    /* create Q-matrix */
    DCProgs::t_rmatrix qmatrix(k ,k);
    for (int i=0; i < k; i++){
        for (int j=0; j < k;j++){
            qmatrix(i,j) = Q[(j*k)+i];
        }
    }

    if (t < 0) {
        std::cerr << "t cannot be negative" << std::endl;
        return 1;
    }
     
    DCProgs::t_rmatrix expT;
    try {
        // Create missed-events G
       DCProgs::IdealG idealG(DCProgs::QMatrix(qmatrix,nopen));
        if (isAF)
            expT = idealG.af(t);
        else
            expT = idealG.fa(t);
         
    } catch (const DCProgs::errors::Root& e) {
        std::cerr << "[WARN]: DCProgs - Error thrown in DCProgs::IdealG\n";
        std::cerr << e.what() << std::endl; 
        return 1;      
    }
     
    /* return the results */
      
    int nclose = k - nopen;
    int index;
            
    for ( index = 0; index < nopen*nclose; index++ ) {
        GXYt[index] = expT(index);
    }
     
    return 0;   
     
} 

extern "C" int missed_events_likelihood(double* jl_likelihood , double* jl_bursts, size_t interval_length, size_t* burst_lengths, size_t burst_number, double* Q , size_t nopen, size_t k, double tau, double tcrit, bool useChs) {
    // jl_likelihood -> double*: to hold result
    // jl_bursts -> double*: array of all intervals in the burst set (burst_array_length elements)
    // interval_length -> int: length of jl_bursts, number of intervals in all the bursts 
    // burst_number -> size_t: number of bursts
    // burst_lengths -> int*: array of the lengths of all bursts (burst_number elements)
    // Q -> double*: Q-matrix
    // nopen -> int: number of open states in the mechanism
    // k -> int: number of states in the mechanism
    // tau -> double: resolution time for the burst set
    // tcrit -> double: tcritical time for seperating the bursts
    // useChs -> bool: whether to use chs vectors for low concnetration recordings

  DCProgs::t_Bursts bursts;

  size_t interval_count=0;
  for (size_t b = 0; b < burst_number; b++){
      DCProgs::t_Burst dburst;
      size_t elem;
      for ( elem = 0; elem < burst_lengths[b]; elem++){
          const double interval = jl_bursts[interval_count];
          dburst.push_back(interval);
          interval_count++;
      }
      bursts.push_back(dburst);
  }
  
  //check the number of bursts equals burst_number and all the intervals have been accounted for
  if (interval_count != interval_length) {
      std::cout << "[WARN]: interval count does not equal the length of the intervals array" << interval_count << " vs " << interval_length << std::endl;
  }

  if( bursts.size() != burst_number) {
      std::cout << "[WARN]: burst count does not match the number of bursts parsed " << bursts.size() << " vs " << burst_number << std::endl;
  }

  if(useChs == 0)
      tcrit = -tcrit; //use equilibrium vectors
  
  DCProgs::t_rmatrix qmatrix(k ,k);
  for (size_t i=0; i < k; i++){
      for (size_t j=0; j < k;j++){
          qmatrix(i,j) = Q[(j*k)+i];
      }
  }

  try {
      DCProgs::Log10Likelihood likelihood( bursts, nopen, tau, tcrit);
      double lik = likelihood(qmatrix);
      *jl_likelihood = lik * log(10);
  } catch (const DCProgs::errors::Root& ex) {
      std::cerr << "[ERROR]: Error thrown in dcprogs" << ex.what() << std::endl;
      *jl_likelihood = -std::numeric_limits<double>::infinity();
      return 1;
  } catch(...) { 
      std::cerr << "Some error" <<std::endl;
      return 1;
  }
  
  //dcprogs can safely return Inf or NaN so dead with these cases
  if (fabs(*jl_likelihood) == std::numeric_limits<double>::infinity() || std::isnan(*jl_likelihood)){
      *jl_likelihood =  -std::numeric_limits<double>::infinity();
      return 1;
  }
  return 0;
}

extern "C" int dcprogs_julia_likelihood(double* jl_likelihood){
   DCProgs::t_Bursts bursts{
     {0.1, 0.2, 0.1},                  /* 1st burst */
     {0.2},                            /* 2nd burst */
     {0.15, 0.16, 0.18, 0.05, 0.1}     /* 3rd burst */
  };
  
  DCProgs::Log10Likelihood likelihood( bursts, /*nopen=*/2, /*tau=*/1e-2,
                                       /*tcrit=*/DCProgs::quiet_nan );
  DCProgs::t_rmatrix qmatrix(5 ,5);
  qmatrix << -3050,        50,  3000,      0,    0, 
            2./3., -1502./3.,     0,    500,    0,  
               15,         0, -2065,     50, 2000,  
                0,     15000,  4000, -19000,    0,  
                0,         0,    10,      0,  -10;


  //std::cout << likelihood << std::endl;
  *jl_likelihood = likelihood(qmatrix) * log(10);
  //std::cout << *jl_likelihood << std::endl;
  return 0;

}

