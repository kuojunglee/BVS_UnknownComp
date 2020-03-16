//
//  UnknownCompFMR.h
//  
//
//  Created by Kuo-Jung Lee on 7/25/16.
//
//

//#include <math.h>

#include <stdio.h>
//#include <algorithm>
//#include <assert.h>
#include <cmath>
#include <ctime>    // For time()
//#include <cstdlib>  // For srand() and rand()
//#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <iomanip>
//#include <list>
//#include <limits>
#include <vector>
#include <string>
//#include <sstream>
#include <algorithm>

#include <RcppArmadillo.h>
using namespace std;
using namespace arma;
using namespace Rcpp;



#ifndef UnknownCompFMR_h
#define UnknownCompFMR_h


class UnknownCompFMR{
private:
    int Num_of_iterations, Num_of_iterations_Inside;
    int Max_num_of_groups, Num_of_obs, Num_of_covariates;
    bool Robust, RJMCMC, CompChg;
    colvec Y;
    mat X;
    vec alpha_g;
    double c2, sigma2_a, sigma2_b, nu, tuning_para, lambda;
    vec Group_Info, z;
    cube r_samples;
    mat rho_samples;
    mat z_samples, omega_samples;
    List InitialValues, GivenValues, Data;
    vec dgj, acc_rate;
public:
    UnknownCompFMR(int iNum_of_iterations, int iMax_num_of_groups, colvec cY, mat mX, List list_InitialValues, List list_GivenValues);
    UnknownCompFMR();
    UnknownCompFMR(int iNum_of_iterations, int iNum_of_iterations_Inside, int iMax_num_of_groups,  List list_Data, List list_InitialValues, List list_GivenValues, bool bRobust, bool bRJMCMC);
    vec rdirichlet_cpp(vec alpha_m);
    vec beta_distance(mat beta_mat, int g_index);
    int min_beta_distance_group_index(mat beta_mat, int g_index);
    double Log_Posterior_Calculation(mat r_sample, vec z_sample, vec rho_sample, vec group_info, vec omega_sample);
    mat Calculate_Beta_Est(vec z, mat r_sample, vec omega_sample);
    vec Update_Z(vec z_sample_input, vec rho_sample_input, mat r_sample_input, vec omega_sample);
    vec Updata_rho();
    mat Updata_gamma(vec z_sample_input, mat r_sample_input, vec omega_sample);
    vec Update_weight(mat r_sample, vec z_sample, vec rho_sample, vec omega_sample); 
    void Generating_Eliminating(int iter, vec z_sample, mat r_sample, vec rho_sample, vec omega_sample);
    void Splitting_Merging(int iter, vec z_sample, mat r_sample, vec rho_sample, vec omega_sample);
    //void RemoveLabelSwitching(List PosteriorSamplesForLabelSwitching);
    SEXP MCMC_Procedure();
    //SEXP MCMC_Procedure_Robust();
    
};


#endif /* UnknownCompFMR_h */


