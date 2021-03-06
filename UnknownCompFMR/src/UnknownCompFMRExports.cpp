// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "UnknownCompFMR.h"

RcppExport SEXP WorkUnknownCompFMR(SEXP i_Num_of_iterations, SEXP i_Num_of_iterations_Inside, SEXP i_Max_num_of_groups, SEXP list_Data, SEXP list_InitialValues, SEXP list_GivenValues, SEXP b_Robust, SEXP b_RJMCMC)
{
    List PosteriorSamples;
    
    List lInitialValues(list_InitialValues), lData(list_Data), lGivenValues(list_GivenValues);
    
    int iNum_of_iterations = Rcpp::as<int> (i_Num_of_iterations);
    
    int iNum_of_iterations_Inside = Rcpp::as<int> (i_Num_of_iterations_Inside);
    
    int iMax_num_of_groups = Rcpp::as<int> (i_Max_num_of_groups);
    
    bool bRobust = Rcpp::as<bool> (b_Robust);
    bool bRJMCMC = Rcpp::as<bool> (b_RJMCMC);
    
    UnknownCompFMR DoUnknownCompFMR(iNum_of_iterations, iNum_of_iterations_Inside, iMax_num_of_groups, lData, lInitialValues, lGivenValues, bRobust, bRJMCMC);
    
    PosteriorSamples = DoUnknownCompFMR.MCMC_Procedure();
    
    return (PosteriorSamples);
    
}


