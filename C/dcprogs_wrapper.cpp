#include <iostream>
#include "limits.h"
#include "math.h"
#include "likelihood.h"

extern "C" int add(int a, int b) {
    int c = a + b;
    return c;
}

extern "C" int missed_events_likelihood(double* jl_likelihood , double* jl_bursts, size_t interval_length, size_t* burst_lengths, size_t burst_number, double* Q , int nopen, int k, double tau, double tcrit, bool useChs) {
    // jl_likelihood -> double*: to hold result
    // jl_bursts -> double*: array of all intervals in the burst set (burst_array_length elements)
    // interval_length -> size_t: length of jl_bursts, number of intervals in all the bursts 
    // burst_number -> size_t: number of bursts
    // burst_lengths -> size_t*: array of the lengths of all bursts (burst_number elements)
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
    //std::cout << b << " Burst lengths " << burst_lengths[b] << std::endl;
    for ( elem = 0; elem < burst_lengths[b]; elem++){
      const double interval = jl_bursts[interval_count];
      dburst.push_back(interval);
      interval_count++;
      //std::cout << "\t" << interval << std::endl;
    }
    bursts.push_back(dburst);
  }
  std::cout << "interval count " << interval_count << std::endl;
  //check the number of bursts equals burst_number and all the intervals have been accounted for
  if (interval_count != interval_length) {
    std::cout << "[WARN]: interval count does not equal the length of the intervals array" << interval_count << " vs " << interval_length << std::endl;
  }

  if( bursts.size() != burst_number) {
    std::cout << "[WARN]: burst count does not match the number of bursts parsed " << bursts.size() << " vs " << burst_number << std::endl;
  }


  DCProgs::Log10Likelihood likelihood( bursts, nopen, tau,
                                       tcrit);
  
  DCProgs::t_rmatrix qmatrix(k ,k);
  for (int i=0; i < k; i++){
      for (int j=0; j < k;j++){
          qmatrix(i,j) = Q[(j*k)+i];
          //printf("[%d][%d] = %.16f\n",i,j,qmatrix(i,j));
      }
  }

  //std::cout << likelihood << std::endl;
  *jl_likelihood = likelihood(qmatrix) * log(10);
  //std::cout << *jl_likelihood << std::endl;
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

/*
extern "C" int dcprogs_likelihood(const gmcmc_burst_dataset * dataset, double * qmatrix, int k , size_t ldq,int kA, double tau, double tcrit, size_t  set_no, int useChs, double * likelihood) {
    //first, unwrap the data


    DCProgs::t_Bursts bursts;

    size_t m = gmcmc_dataset_get_num_bursts(dataset,set_no);
    double burst_time = 0;
    int interval_count=0;
    //printf("Set %zu has %zu bursts\n", set_no, m);
    for (size_t b = 0; b < m; b++){
        size_t * burst_length = gmcmc_dataset_get_burst_number(dataset,b,set_no);
        //printf("Burst %zu -> %zu elements\n", set_no, *burst_length);

        //put the bursts in the vector format required by dcprogs

        DCProgs::t_Burst dburst;
        size_t elem;
        const double * burst = gmcmc_dataset_get_burst(dataset , b, set_no);
        for ( elem = 0; elem < *burst_length; elem++){
            dburst.push_back(burst[elem]);
            burst_time+=burst[elem];
            interval_count++;
            //printf("za%f\n", burst[elem]);
        }
        bursts.push_back(dburst);
        free(burst_length); //instantiated with malloc in the C code
    }

    //fill in the q_matrix
    //printf("In wrapper: ldq -> %zu \n",ldq);
    DCProgs::t_rmatrix matrix(k , k);
    for (int i=0; i < k; i++){
        for (int j=0; j < k;j++){
            matrix(i,j) = qmatrix[(j*ldq)+i];
            //printf("[%d][%d] = %.16f\n",i,j,matrix(i,j));
        }
    }
    // printf("End of Q\n");
    if (! useChs > 0)
        tcrit=-1;//-tcrit;

    DCProgs::Log10Likelihood dc_likelihood(bursts, kA, tau, tcrit);
    DCProgs::t_real result = 0;
    int error = 0;
    try {
        result = dc_likelihood(matrix);
    }
    catch (std::exception& e) {
        printf("Exception thrown\n");
        printf(e.what());
        printf("Q-matrix in dcprogs\n");
        std::cout << matrix << std::endl;
        error = -1;
    }

    if (fabs(result) == std::numeric_limits<double>::infinity() || std::isnan(result)){
        error=-1;
    }

    *likelihood = result*log(10);
    *printf ("parsed bursts %zu\n", bursts.size());
    printf("t-crit %.10f\n",tcrit);
    printf("t-res %.10f\n",tau);
    printf("kA %i\n",kA);
    printf("Burst time %.16f\n",burst_time);
    printf("Interval count %i\n",interval_count);
    printf("Likelihood %.17f\n",*likelihood);

    return error;
}
*/
