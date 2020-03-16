// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "UnknownCompFMR.h"
RNGScope scope;
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]


// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

UnknownCompFMR::UnknownCompFMR(int iNum_of_iterations, int iMax_num_of_groups, colvec cY, mat mX, List list_InitialValues, List list_GivenValues)
{
    Num_of_iterations = iNum_of_iterations;
    Max_num_of_groups = iMax_num_of_groups;
    Y = cY;
    X = mX;
    Num_of_obs = cY.n_elem;
    Num_of_covariates = mX.n_cols;
    InitialValues = list_InitialValues;
    GivenValues = list_GivenValues;
    //tau2 = as<vec>(GivenValues["tau2.vec"]);
    c2 = as<double>(GivenValues["c2"]);
    Group_Info = as<vec>(GivenValues["group.info"]);
    sigma2_a = as<double>(GivenValues["sigma2.a"]);
    sigma2_b = as<double>(GivenValues["sigma2.b"]);
    lambda = as<double>(GivenValues["lambda"]);
    Rcout<<"lambda = " << lambda << endl;
}

UnknownCompFMR::UnknownCompFMR(int iNum_of_iterations, int iNum_of_iterations_Inside, int iMax_num_of_groups,  List list_Data, List list_InitialValues, List list_GivenValues, bool bRobust, bool bRJMCMC){
    
    CompChg = false;
    
    InitialValues = list_InitialValues;
    Data = list_Data;
    GivenValues = list_GivenValues;
    Rcout << "read data." << endl;
    //Group_Info = as<vec>(InitialValues["group.info"]);
    z = as<vec>(InitialValues["z"]);
    //r_sample = as<mat>(InitialValues["r.sample"]);
    //rho_sample = as<vec>(InitialValues["rho.sample"]);
    Max_num_of_groups = iMax_num_of_groups;
    Num_of_iterations = iNum_of_iterations;
    Num_of_iterations_Inside = iNum_of_iterations_Inside; //50;
    Robust = bRobust;
    RJMCMC = bRJMCMC;
    
    Y = as<colvec>(Data["Y"]);
    X = as<mat>(Data["X"]);
    
    Group_Info.set_size(Max_num_of_groups);
    Group_Info.fill(datum::nan);
    vec group_info_value = unique(z);
    uvec group_info_num = hist(z,group_info_value);
    Group_Info(conv_to< uvec >::from(group_info_value)) = conv_to< vec >::from(group_info_num);
    
    //cout << "z = " << z << endl;
    //cout <<"Group_Info = " << Group_Info << endl;

    
    Num_of_covariates = X.n_cols;
    Num_of_obs = Y.n_elem;
    
    r_samples.set_size(Max_num_of_groups, Num_of_covariates, Num_of_iterations);
    rho_samples.set_size(Max_num_of_groups, Num_of_iterations);
    z_samples.set_size(Num_of_obs, Num_of_iterations);
    
    
    r_samples.fill(datum::nan);
    rho_samples.fill(0);
    
    r_samples.slice(0) = as<mat>(InitialValues["r.sample"]);
    rho_samples.col(0) = as<vec>(InitialValues["rho.sample"]);
    z_samples.col(0) = as<vec>(InitialValues["z"]);
    
    alpha_g = as<vec>(GivenValues["alpha.g"]);
    //tau2 = as<vec>(GivenValues["tau2"]);
    c2 = as<double>(GivenValues["c2"]);
    sigma2_a = as<double>(GivenValues["sigma2.a"]);
    sigma2_b = as<double>(GivenValues["sigma2.b"]);
    
    nu = as<double>(GivenValues["nu"]);
    tuning_para = as<double>(GivenValues["tuning.para"]);
    dgj = as<vec>(GivenValues["dgj"]);
    
    omega_samples.set_size(Num_of_obs, Num_of_iterations);
    omega_samples.col(0) = as<vec>(InitialValues["omega.sample"]);;
    lambda = as<double>(GivenValues["lambda"]);
    
    Rcout << "end read data."<<endl;
}



vec UnknownCompFMR::rdirichlet_cpp(vec alpha_m)
{
    int distribution_size = alpha_m.n_elem;

    vec distribution = zeros(distribution_size);
    
    for (int j = 0; j < distribution_size; j++){
        double cur = Rf_rgamma(alpha_m[j],1.0);
        distribution(j) = cur;
    
    }
    distribution = distribution/accu(distribution);

    return(distribution);
}

double UnknownCompFMR::Log_Posterior_Calculation(mat r_sample, vec z_sample, vec rho_sample, vec group_info, vec omega_sample)
{
    double alpha_g = 1;
    
    double log_post_prob = 0, SSE_g=0;
    
    uvec groups_index_tmp =  find_finite(group_info); //unique(z_sample);
    
    int num_of_groups_tmp = groups_index_tmp.n_elem;
    
    colvec y_g;
    mat r_g;
    int n_g = 0;
    
    mat Omega_g, Omega_g_inv, Weight_inv;
    
    
    mat x_g, Lambda_g, x_g_1;
    double h_g;
    int g;
    
    for(int g_tmp =0; g_tmp<groups_index_tmp.n_elem; g_tmp++){
        //Rcout << "Log_Posterior_Calculation : g_temp =" << g_tmp << endl;
        g = groups_index_tmp(g_tmp);
        y_g = Y.elem( find(z_sample ==g) );
        if(y_g.n_elem!=0){
            r_g = r_sample.row(g);
            n_g = accu(z_sample==g);
            h_g = (n_g+sigma2_a)/2.;
            x_g = X.rows(find(z_sample ==g));
            x_g_1 = x_g.cols(find(r_g>0));
            
            if(y_g.is_empty())
                Weight_inv = zeros(1,1);
            else
                Weight_inv = diagmat(1./omega_sample(find(z_sample ==g)));
            
            
            
            
            if(accu(r_g)>0){
                
                if(y_g.is_empty()){
                    Lambda_g = eye(accu(r_g),accu(r_g))/c2;
                    SSE_g = 0.;
                }
                else{
                    Lambda_g = x_g_1.t()*Weight_inv*x_g_1 + eye(x_g_1.n_cols, x_g_1.n_cols)/c2;
                    SSE_g = as_scalar(y_g.t()*(Weight_inv - Weight_inv*x_g_1*Lambda_g.i()*x_g_1.t()*Weight_inv)*y_g);
                }
                
                
                //if(Robust)
                    //log_post_prob  += (0.5*(nu+1.)+1.)*log(prod(1./omega_sample(find(z_sample ==g)))) - 0.5*accu(nu/omega_sample(find(z_sample ==g))) + 0.5*nu*log(nu*0.5) - lgamma(0.5*nu) +  0.5*sigma2_a*log(0.5*sigma2_b) -lgamma(0.5*sigma2_a) + 0.5*n_g*log(rho_sample(g)/(2.*datum::pi))- 0.5*accu(r_g)*log(c2) -0.5*log(det(Lambda_g))+ lgamma(h_g)-h_g*log(0.5*(sigma2_b+SSE_g))+ (alpha_g+1.)*log(1.*rho_sample(g)) - lgamma(alpha_g) + Num_of_covariates*log(0.5);
                //else
                log_post_prob  += 0.5*sigma2_a*log(0.5*sigma2_b) -lgamma(0.5*sigma2_a) + 0.5*n_g*log(rho_sample(g)*rho_sample(g)/(2.*datum::pi)) - 0.5*accu(r_g)*log(c2) - 0.5*log(det(Lambda_g)) + lgamma(h_g) - h_g*log(0.5*(sigma2_b+SSE_g))+ (alpha_g+1.)*log(1.*rho_sample(g)) - lgamma(alpha_g) + accu(log((dgj(find(r_g>0.))))) + accu( log ( (1.- dgj(find(r_g<1.)) )) );
                
                //Rcout << 1.- dgj(find(r_g<1.)) << endl;
            }
            else{
                SSE_g = as_scalar(y_g.t()*(Weight_inv)*y_g);
                //if(Robust)
                    //log_post_prob  += (0.5*(nu+1.)+1.)*log(prod(1./omega_sample(find(z_sample ==g)))) - 0.5*accu(nu/omega_sample(find(z_sample ==g))) + 0.5*nu*log(nu*0.5) - lgamma(0.5*nu) + 0.5*sigma2_a*log(0.5*sigma2_b) -lgamma(0.5*sigma2_a) + 0.5*n_g*log(rho_sample(g)/(2.*datum::pi))+ lgamma(h_g)-h_g*log(0.5*(sigma2_b+SSE_g))+ (alpha_g+1.)*log(1.*rho_sample(g)) - lgamma(alpha_g) + Num_of_covariates*log(0.5);
                //else
                log_post_prob  += 0.5*sigma2_a*log(0.5*sigma2_b) -lgamma(0.5*sigma2_a) + 0.5*n_g*log(rho_sample(g)*rho_sample(g)/(2.*datum::pi)) + lgamma(h_g)-h_g*log(0.5*(sigma2_b+SSE_g))+ (alpha_g+1.)*log(1.*rho_sample(g)) - lgamma(alpha_g) + accu(log((dgj(find(r_g>0.)))))+ accu( log ( (1.- dgj(find(r_g<1.)) )) );

            }
        }
        else{
            log_post_prob  += 0.5*sigma2_a*log(0.5*sigma2_b) -lgamma(0.5*sigma2_a) + lgamma(0.5*sigma2_a)-0.5*sigma2_a*log(0.5*sigma2_b) + (alpha_g+1.)*log(1.*rho_sample(g)) - lgamma(alpha_g) + accu(log ( (1.- dgj)) ) ;
            //Rcout << 1.- dgj << endl;
        }
        
    }
    
    
    return (log_post_prob  + lgamma(num_of_groups_tmp*alpha_g));
}

/*
double UnknownCompFMR::Log_Posterior_Calculation_Robust(mat r_sample, vec z_sample, vec rho_sample, vec group_info, vec omega_sample)
{
    double alpha_g = 1;
    
    double log_post_prob = 0, SSE_g=0;
    
    uvec groups_index_tmp =  find_finite(group_info); //unique(z_sample);
    
    int num_of_groups_tmp = groups_index_tmp.n_elem;
    
    colvec y_g;
    mat r_g;
    int n_g;
    
    mat Omega_g, Omega_g_inv, Weight_inv;
    
    
    mat x_g, Lambda_g, x_g_1;
    double h_g;
    int g;
    
    for(int g_tmp =0; g_tmp<groups_index_tmp.n_elem; g_tmp++){
        //Rcout << "Log_Posterior_Calculation : g_temp =" << g_tmp << endl;
        g = groups_index_tmp(g_tmp);
        y_g = Y.elem( find(z_sample ==g) );
        if(y_g.n_elem!=0){
            r_g = r_sample.row(g);
            n_g = accu(z_sample==g);
            h_g = (n_g+sigma2_a)/2.;
            x_g = X.rows(find(z_sample ==g));
            x_g_1 = x_g.cols(find(r_g>0));
            
            if(y_g.is_empty())
                Weight_inv = zeros(1,1);
            else
                Weight_inv = diagmat(1./omega_sample(find(z_sample ==g)));

            
            
            
            if(accu(r_g)>0){
            
                if(y_g.is_empty()){
                    Lambda_g = eye(accu(r_g),accu(r_g))/c2;
                    SSE_g = 0.;
                }
                else{
                    Lambda_g = x_g_1.t()*Weight_inv*x_g_1 + eye(x_g_1.n_cols, x_g_1.n_cols)/c2;
                    SSE_g = as_scalar(y_g.t()*(Weight_inv - Weight_inv*x_g_1*Lambda_g.i()*x_g_1.t()*Weight_inv)*y_g);
                }
                
                
                
                log_post_prob  += log(prod(1./omega_sample(find(z_sample ==g)))) + 0.5*sigma2_a*log(0.5*sigma2_b) -lgamma(0.5*sigma2_a) + 0.5*n_g*log(rho_sample(g)/(2.*datum::pi))- 0.5*accu(r_g)*log(c2) -0.5*log(det(Lambda_g))+ lgamma(h_g)-h_g*log(0.5*(sigma2_b+SSE_g))+ (alpha_g+1.)*log(1.*rho_sample(g)) - lgamma(alpha_g) + Num_of_covariates*log(0.5);
            }
            else{
                SSE_g = as_scalar(y_g.t()*(Weight_inv)*y_g);;
                log_post_prob  += log(prod(1./omega_sample(find(z_sample ==g)))) + 0.5*sigma2_a*log(0.5*sigma2_b) -lgamma(0.5*sigma2_a) + 0.5*n_g*log(rho_sample(g)/(2.*datum::pi))+ lgamma(h_g)-h_g*log(0.5*(sigma2_b+SSE_g))+ (alpha_g+1.)*log(1.*rho_sample(g)) - lgamma(alpha_g) + Num_of_covariates*log(0.5);

            }
        }
        else
            log_post_prob  += (alpha_g+1.)*log(1.*rho_sample(g)) - lgamma(alpha_g) + Num_of_covariates*log(0.5);
            
    }
                                                                                                            
    
    return (log_post_prob  + lgamma(num_of_groups_tmp*alpha_g));
}
 */

vec UnknownCompFMR::beta_distance(mat beta_mat, int g_index)
{
    int nbeta=beta_mat.n_rows;
    vec dist = zeros<vec>(nbeta-1);
    int i = 0;
    while(i<nbeta){
        if(i!=g_index)
            dist(i) = accu(square(beta_mat.row(g_index)-beta_mat.row(i)));
        i++;
    }
    //beta.dist.g1 = apply(beta.est[-g1, ], 1, function(x) sum((beta.est[g1, ]-x)^2))
    return dist;
}


int UnknownCompFMR::min_beta_distance_group_index(mat beta_mat, int g_index)
{
    //cout << g_index << "\t" << "beta_mat = " << beta_mat <<endl;
    int nbeta=beta_mat.n_rows;
    int selected_index;
    vec dist = zeros<vec>(nbeta);
    double dist0;
    
    for(int i=0; i<nbeta; i++)
        dist(i) = accu(square(beta_mat.row(i)-beta_mat.row(g_index)));
    
    dist(find_nonfinite(dist)).fill(datum::inf);
    //cout << "dist = " << dist << endl;
    if(g_index==0)
        selected_index = 1;
    else
        selected_index = g_index-1;
    dist0 = dist(selected_index);
    for(int i=0; i<nbeta; i++)
        if((i!=g_index) && (dist(i) < dist0)){
            selected_index = i;
            dist0 = dist(selected_index);
        }
    return selected_index;
}

mat UnknownCompFMR::Calculate_Beta_Est(vec z_sample, mat r_sample, vec omega_sample)
{
    //beta.est.tmp = matrix(NA, max.num.of.groups, num.of.covariates)
    colvec y_g;
    mat r_g;
    mat beta_est_tmp(Max_num_of_groups, X.n_cols);
    beta_est_tmp.fill(datum::nan);
    vec beta_est_tmp1= zeros<vec>(X.n_cols);
    vec groups_index_tmp = unique(z);
    mat x_g, x_g_1, Lambda_g, Weight_inv;
    int g;
    //rownames(beta.est.tmp) = 1:max.num.of.groups
     for(int g_tmp=0; g_tmp<groups_index_tmp.n_elem; g_tmp++){
         g = groups_index_tmp(g_tmp);
         beta_est_tmp1.fill(0);
         y_g = Y.elem( find(z ==g) );
         r_g = r_sample.row(g);
         x_g = X.rows(find(z ==g));
         x_g_1 = x_g.cols(find(r_g>0));
         if(accu(r_g)>0){
             Weight_inv = diagmat(1./omega_sample(find(z_sample ==g)));
             Lambda_g = x_g_1.t()*Weight_inv*x_g_1 + eye(x_g_1.n_cols,x_g_1.n_cols)/c2;
             beta_est_tmp1(find(r_g>0))= Lambda_g.i()*x_g_1.t()*Weight_inv*y_g;
             //cout << "Lambda_g.i()*x_g_1.t()*y_g ＝"<< Lambda_g.i()*x_g_1.t()*y_g << endl;
             beta_est_tmp.row(g) = beta_est_tmp1.t();
         }
         else
             beta_est_tmp.row(g).zeros();
    }
    return(beta_est_tmp);
}


vec UnknownCompFMR::Update_Z(vec z_sample_input, vec rho_sample_input, mat r_sample_input, vec omega_sample)
{
    //Rcout << "Start Updating Z." <<endl;
    //cout <<"Z: Group_Info = " << endl << Group_Info << endl;
    //vec Group_Info_tmp = Group_Info(find_finite(Group_Info));

    //cout <<"Z: Group sum = " << accu(Group_Info_tmp) << endl;
    double prob_z_tmp = 0;
    
    double SSE_g=0;
    
    mat r_sample = r_sample_input; //r_samples.slice(iter);
    vec rho_sample = rho_sample_input; //rho_samples.col(iter);
    
    
    
    vec groups_index_tmp = unique(z_sample_input);
    
    
    vec prob_z = zeros<vec>(groups_index_tmp.n_elem);
    colvec y_g;
    mat r_g;
    int n_g;
    
    mat Omega_g, Omega_g_inv, Weight_inv;
    
    mat x_g, Lambda_g, x_g_1;
    double h_g;
    int g;
    
    uvec group_indice_nozeros = conv_to<uvec>::from(unique(z_sample_input));
    //Rcout << "Update Z: r_sample_input = " << r_sample_input << endl;
    //Rcout << "Start Updating Z. Check 1" <<endl;
    for(int i=0; i<Y.n_elem; i++){
        //Rcout << "Origin z_sample_input(i) = " << z_sample_input(i) << endl;
        for(int g_tmp=0; g_tmp<groups_index_tmp.n_elem; g_tmp++){
            g = groups_index_tmp(g_tmp);
            //cout << "g = " << g << endl;
            z_sample_input(i) = g;
            y_g = Y.elem( find(z_sample_input ==g) );
            r_g = r_sample.row(g);
            n_g = y_g.n_elem; //accu(z_sample_input==g);
            
            h_g = (n_g+sigma2_a)/2.;
            x_g = X.rows(find(z_sample_input ==g));
            //Weight_inv = diagmat(1./omega_sample(find(z_sample_input ==g)));
            if(y_g.is_empty())
                Weight_inv = zeros(1,1);
            else
                Weight_inv = diagmat(1./omega_sample(find(z_sample_input ==g)));
            
            if(accu(r_g)>0){
                x_g_1 = x_g.cols(find(r_g>0));
                
                if(y_g.is_empty()){
                    Lambda_g = eye(accu(r_g),accu(r_g))/c2;
                    SSE_g = 0.;
                }
                else{
                    Lambda_g = x_g_1.t()*Weight_inv*x_g_1 + eye(x_g_1.n_cols, x_g_1.n_cols)/c2;
                    SSE_g = as_scalar(y_g.t()*(Weight_inv - Weight_inv*x_g_1*Lambda_g.i()*x_g_1.t()*Weight_inv)*y_g);
                }
                
                //Lambda_g = x_g_1.t()*x_g_1 + eye(accu(r_g),accu(r_g))/c2;
                
                
                //SSE_g = as_scalar(y_g.t()*(eye(n_g,n_g) - x_g_1*Lambda_g.i()*x_g_1.t())*y_g);
                
                //cout << "SSE_g = " << SSE_g << endl;
                //cout << "0.5*n_g*log(rho_sample(g)/(2.*datum::pi)) = " << 0.5*n_g*log(rho_sample(g)/(2.*datum::pi)) << endl;
                //cout << "r_g = " << r_g << endl;
                //cout << "0.5*accu(r_g) = " << 0.5*accu(r_g) << endl;
                //cout << "0.5*log(det(Lambda_g)) = " << 0.5*log(det(Lambda_g)) << endl;
     
                prob_z_tmp  = 0.5*n_g*log(rho_sample(g)*rho_sample(g)/(2.*datum::pi))- 0.5*accu(r_g)*log(c2)-0.5*log(det(Lambda_g))+ lgamma(h_g)-h_g*log(0.5*(sigma2_b+SSE_g))+ (alpha_g(g)+1.)*log(rho_sample(g)) +  accu(log((dgj(find(r_g>0.)))))+ accu( log ( (1.- dgj(find(r_g<1.)) )) );
            }
            else{
                SSE_g = as_scalar(y_g.t()*(Weight_inv)*y_g);
                prob_z_tmp  = 0.5*n_g*log(rho_sample(g)*rho_sample(g)/(2.*datum::pi))+ lgamma(h_g)-h_g*log(0.5*(sigma2_b+SSE_g))+ accu(log ( (1.- dgj)) );
            }
            //cout << "prob_z_tmp = " << prob_z_tmp << endl;
            prob_z(g_tmp) = prob_z_tmp;
            
        }
        //Rcout <<"prob_z = " << prob_z <<endl;
        prob_z = prob_z - median(prob_z);
        prob_z = exp(prob_z)/accu(exp(prob_z));
        //Rcout << "prob_z = " << prob_z << endl;
        z_sample_input(i) = Rcpp::RcppArmadillo::sample(group_indice_nozeros, 1, true, as<NumericVector>(wrap(prob_z)))(0);
        //Rcout << " z_sample_input(i) = " << z_sample_input(i) << "\t" << group_indice_nozeros << endl;
    }
    //Rcout << "End Update Z" << endl;
    
    return z_sample_input;
}


vec UnknownCompFMR::Updata_rho()
{
    
    uvec group_info_value = find_finite(Group_Info);// unique(z);
    vec group_info_num = Group_Info(group_info_value); //hist(z,group_info_value);
    vec rho_sample = zeros<vec>(Group_Info.n_elem);
    
    //Rcout<< "group_info_num = " <<endl << group_info_num << endl;
    
    //uvec group_info_value_uvec = conv_to<uvec>::from(group_info_value);
    vec group_info_alpha = alpha_g(group_info_value);
    
    rho_sample.fill(0.);
    for(int j=0; j < Num_of_iterations_Inside; j++)
        rho_sample(group_info_value)= rdirichlet_cpp( (group_info_num + group_info_alpha));
    
    
    //Rcout << "Update rho rho_sample = " << rho_sample << endl;
    //Rcout << "Update rho" << endl;
    return rho_sample;
}

mat UnknownCompFMR::Updata_gamma(vec z_sample_input, mat r_sample_input, vec omega_sample)
{
    //Rcout << "Start Updating gamma" << endl;
    uvec group_indice_nozeros = conv_to<uvec>::from(unique(z_sample_input));  //find(Group_Info>0);
    colvec y_g;
    mat r_g, r_g_0, r_g_1, x_g_0, x_g_1;
    int n_g;
    
    double SSE_g_0=0, SSE_g_1=0;
    
    mat Omega_g, Omega_g_inv, r_sample;
    r_sample = r_sample_input;
    
    mat x_g, Lambda_g_0, Lambda_g_1, Weight_inv;
    double h_g;
    
    //mat AA, BB;
   
    double prob_r_tmp_0 = 0, prob_r_tmp_1=0, prob_r_1 = 0;
    int g;
    for(int g_tmp=0; g_tmp<group_indice_nozeros.n_elem; g_tmp++){
        //Rcout << "Updata_gamma : g_temp =" << g_tmp << endl;
        g = group_indice_nozeros(g_tmp);
        y_g = Y.elem( find(z_sample_input ==g) );
        //Rcout << g_tmp << "\t" << y_g << endl;
        r_g = r_sample.row(g);
        n_g = y_g.n_elem; //accu(z_sample_input==g);
        
        h_g = 0.5*(n_g+sigma2_a);
        x_g = X.rows(find(z_sample_input ==g));
        
        if(y_g.is_empty()){
            Weight_inv = zeros(1,1);
            y_g = zeros(1);
        }
        else
            Weight_inv = diagmat(1./omega_sample(find(z_sample_input ==g)));
            
        
        for(int p=0; p<Num_of_covariates; p++){
            r_g(p) = 0;
            r_g_0 = r_g;
            
            r_g(p) = 1;
            r_g_1 = r_g;
            
            x_g_0 = x_g.cols(find(r_g_0>0));
            x_g_1 = x_g.cols(find(r_g_1>0));
            
            if(y_g.is_empty()){
                Lambda_g_1 = eye(accu(r_g_1),accu(r_g_1))/c2;
                SSE_g_1 = 0.;
            }
            else{
                Lambda_g_1 = x_g_1.t()*Weight_inv*x_g_1 + eye(x_g_1.n_cols, x_g_1.n_cols)/c2;
                SSE_g_1 = as_scalar(y_g.t()*(Weight_inv - Weight_inv*x_g_1*Lambda_g_1.i()*x_g_1.t()*Weight_inv)*y_g);
            }
            
            if(accu(r_g_0)>0){

                if(y_g.is_empty()){
                    Lambda_g_0 = eye(accu(r_g_0),accu(r_g_0))/c2;
                    SSE_g_0 = 0.;
                }
                else{
                    Lambda_g_0 = x_g_0.t()*Weight_inv*x_g_0 + eye(x_g_0.n_cols, x_g_0.n_cols)/c2;
                    SSE_g_0 = as_scalar(y_g.t()*(Weight_inv - Weight_inv*x_g_0*Lambda_g_0.i()*x_g_0.t()*Weight_inv)*y_g);
                }
                
                prob_r_tmp_0  =             - 0.5*log(det(Lambda_g_0))-h_g*log(0.5*(sigma2_b+SSE_g_0)) + log(1-dgj(p));
                
                prob_r_tmp_1  = -0.5*log(c2)- 0.5*log(det(Lambda_g_1))-h_g*log(0.5*(sigma2_b+SSE_g_1)) + log(dgj(p));
                
            }
            else{
                
                SSE_g_0 = as_scalar(y_g.t()*(Weight_inv)*y_g);;
                
                prob_r_tmp_0  =                                       -h_g*log(0.5*(sigma2_b+SSE_g_0)) + log(1-dgj(p));
                
                prob_r_tmp_1  = -0.5*log(c2)- 0.5*log(det(Lambda_g_1))-h_g*log(0.5*(sigma2_b+SSE_g_1)) + log(dgj(p));

            }
            
            if(prob_r_tmp_0-prob_r_tmp_1 >100)
                prob_r_1 = 0.;
            else if( prob_r_tmp_0-prob_r_tmp_1 < -100)
                prob_r_1 = 1.;
            else
                prob_r_1 = 1./(1.+exp(prob_r_tmp_0-prob_r_tmp_1));
            
            r_sample(g, p) = Rf_rbinom(1., prob_r_1);
            
        }
    }
    //Rcout << "Update gamma" << endl;
    return r_sample;
    
}


void UnknownCompFMR::Generating_Eliminating(int iter, vec z_sample, mat r_sample, vec rho_sample, vec omega_sample)
{
    //Rcout << "Perform Generating or Eliminating Step." << endl;
    
    //uvec G1_vec = find(Group_Info>0);
    uvec G0_vec = find(Group_Info==0);
    uvec G_vec = find_finite(Group_Info);
    uvec GNaN_vec = find_nonfinite(Group_Info);
    vec Group_Info_tmp;
    
    int G0 = G0_vec.n_elem, G = G_vec.n_elem;
    double bG =0.5, dG = 0.5;
    int G_star;
    
    int new_empty_index, empty_index_to_remove;
    double new_rho;
    mat r_sample_tmp;
    rowvec r_sample_remove;
    //r_sample = r_samples.slice(iter);
    r_sample_tmp = r_sample;
    
    vec rho_sample_tmp;
    double rho_remove;
    //rho_sample = rho_samples.col(iter);
    
    double post_prob_num = 0, post_prob_den = 0;
    double E_prob_num = 0, E_prob_den = 0;
    //double bin_sum_tmp = 0;
    
    bool BIRTH = Rf_rbinom(1., bG)>0 ? true : false;
    if(G==1) BIRTH=true;
    if(G == Max_num_of_groups) BIRTH = false;
    
    if(BIRTH){
        Rcout << "Try generating" << endl;
        G_star = G+1;
        
        if(G==1) bG = 1;
        if(G_star == Max_num_of_groups) dG = 1;
        
        new_empty_index = GNaN_vec(0);
        new_rho = Rf_rbeta(1., G);
        //bin_sum_tmp = 0;


        

        rho_sample_tmp = rho_sample * (1-new_rho);
        rho_sample_tmp(new_empty_index)= new_rho;
        
        r_sample_tmp = r_sample;

        for(int p=0; p<Num_of_covariates; p++){
            r_sample_tmp(new_empty_index, p) = Rf_rbinom(1., 0.5);
            //vec_new_r_sample(p) = Rf_rbinom(1., 0.5);
            //bin_sum_tmp += Rf_dbinom(vec_new_r_sample(p), 1, 0.5, true);
        }
        //r_sample_tmp.row(new_empty_index) = vec_new_r_sample;
        
        //Rcout << "Generating: Group_Info = " <<endl << Group_Info << endl;
        
        Group_Info_tmp = Group_Info;
        Group_Info_tmp(new_empty_index) = 0;
        
        //Rcout << "Generating: Group_Info_tmp = " <<endl << Group_Info_tmp << endl;
       
        post_prob_num = Log_Posterior_Calculation(r_sample_tmp, z, rho_sample_tmp, Group_Info_tmp, omega_sample);
        post_prob_den = Log_Posterior_Calculation(r_sample, z, rho_sample, Group_Info, omega_sample);

        E_prob_num = log(dG) - log(G0+1.) + G*log(1.-new_rho) + Rf_dpois(G_star, lambda, 1);
        E_prob_den = log(bG) + Rf_dbeta(new_rho, 1., G, true) + Num_of_covariates*log(0.5) + Rf_dpois((G_star-1), lambda, 1) ;
        
        //cout <<"Birth (post_prob_num - post_prob_den + E_prob_num-E_prob_den) =" << (post_prob_num - post_prob_den + E_prob_num-E_prob_den) << endl;
        if(log(Rf_runif(0., 1.)) < (post_prob_num - post_prob_den + E_prob_num-E_prob_den) ){
            r_samples.slice(iter+1)  = r_sample_tmp;
            rho_samples.col(iter+1)  = rho_sample_tmp;
            z_samples.col(iter+1) = z = z_sample;
            Group_Info = Group_Info_tmp;
            CompChg = true;
            Rcout<<"A new empty component was born"<<endl;
        }
        else{
            r_samples.slice(iter+1) = r_sample;
            rho_samples.col(iter+1) = rho_sample;
            z_samples.col(iter+1) = z = z_sample;
            CompChg = false;
            Rcout<<"A new empty component was NOT born"<<endl;
        }
    
        //Rcout << "Updated Generating" << endl;
    }
    else{
        Rcout << "Try Eliminating" << endl;
        G_star = G-1;
        
        if(G0>0 && G>1){
            if(G==Max_num_of_groups) dG = 1;
            if(G_star == 1) bG = 1;
        
            
            empty_index_to_remove = Rcpp::RcppArmadillo::sample(G0_vec, 1, true)(0);
            
            //Rcout << "empty_index_to_remove = " << empty_index_to_remove << endl;

            rho_sample_tmp = rho_sample;
            rho_remove = rho_sample(empty_index_to_remove);
            rho_sample_tmp(empty_index_to_remove) = 0;
            rho_sample_tmp = rho_sample_tmp/(1.-rho_remove);
     
            r_sample_tmp = r_sample;
            r_sample_remove = r_sample.row(empty_index_to_remove);
            r_sample_tmp.row(empty_index_to_remove).fill(datum::nan);
            
            //bin_sum_tmp = 0;
            //for(int p=0; p<Num_of_covariates; p++){
            //    bin_sum_tmp += Rf_dbinom(r_sample_remove(p), 1, 0.5, true);
            //}
            
            //Rcout << "Eliminating: Group_Info_tmp = " <<endl << Group_Info << endl;

            Group_Info_tmp = Group_Info;
            Group_Info_tmp(empty_index_to_remove) = datum::nan;
            
            //Rcout << "Eliminating: Group_Info_tmp = " <<endl << Group_Info_tmp << endl;
            
            post_prob_num = Log_Posterior_Calculation(r_sample_tmp, z, rho_sample_tmp, Group_Info_tmp, omega_sample);
            post_prob_den = Log_Posterior_Calculation(r_sample, z, rho_sample, Group_Info, omega_sample);
            
            
            

            E_prob_num = log(bG) + Rf_dbeta(rho_remove, 1., G, true) + Num_of_covariates*log(0.5) + Rf_dpois(G_star, lambda, 1);
            E_prob_den = log(dG) - log(G0) + (G-1)*log(1-rho_remove) + Rf_dpois((G_star+1), lambda, 1);
            
            //cout <<"Death (post_prob_num - post_prob_den + E_prob_num-E_prob_den) =" << (post_prob_num - post_prob_den + E_prob_num-E_prob_den) << endl;
            
            Rcout << "empty remove or not = " << (post_prob_num - post_prob_den + E_prob_num-E_prob_den) << endl;
            Rcout << "post_prob_num =" << post_prob_num << endl;
            Rcout << "post_prob_den =" << post_prob_den << endl;
            Rcout << "E_prob_num =" << E_prob_num << endl;
            Rcout << "E_prob_den =" << E_prob_den << endl;
            
            if(log(Rf_runif(0., 1.)) < (post_prob_num - post_prob_den + E_prob_num-E_prob_den) ){
                r_samples.slice(iter+1) = r_sample_tmp;
                rho_samples.col(iter+1) = rho_sample_tmp;
                z_samples.col(iter+1) = z = z_sample;
                Group_Info(empty_index_to_remove) = datum::nan;
                CompChg = true;
                Rcout << "A new empty component was deleted" << endl;
            }
            else{
                r_samples.slice(iter+1) = r_sample;
                rho_samples.col(iter+1) = rho_sample;
                z_samples.col(iter+1) = z = z_sample;
                CompChg = false;
                Rcout << "A new empty component was NOT deleted" << endl;
            }
        }
        else{
            r_samples.slice(iter+1) = r_sample;
            rho_samples.col(iter+1) = rho_sample;
            z_samples.col(iter+1) = z = z_sample;
            CompChg = false;
            Rcout << "No empty" << endl;
        }
        
        //Rcout << "Updated Eliminating" << endl;
    }
}

void UnknownCompFMR::Splitting_Merging(int iter, vec z_sample, mat r_sample, vec rho_sample, vec omega_sample)
{
    //Rcout << "Perform Splitting and Merging." << endl;
    
    
    uvec G1_vec = find(Group_Info>0);
    uvec G0_vec = find(Group_Info==0);
    uvec G_vec = find_finite(Group_Info);
    uvec GNaN_vec = find_nonfinite(Group_Info);
    vec Group_Info_tmp;
    

    
    int G1 = G1_vec.n_elem, G0 = G0_vec.n_elem, G=G1+G0;
    int G_star, G1_star;
    double bG =0.5, dG = 0.5;
    bool SPLIT;
    rowvec r_sample_g1 = zeros<rowvec>(X.n_cols), r_sample_g2 = zeros<rowvec>(X.n_cols);
    mat r_sample_tmp;
    //r_sample = r_samples.slice(iter);
    
    //cout << "Split&Merge r_sample = " << r_sample << endl;
    
    double P_chos =0., P_alloc =0., u=0.;

    P_chos = 0.; //NULL
    int g1, g2;
    
    //vec z_sample;
    mat beta_est;
    vec sample_integers;
    sample_integers << 0 << 1 << 2;
    
    vec rho_sample_tmp;
    //rho_sample = rho_samples.col(iter);
    
    SPLIT = Rf_rbinom(1., bG)>0 ? true : false;
    if(G1==1) SPLIT=true;
    if(G1 == Max_num_of_groups) SPLIT = false;
    
    int s1, s2;
    
    vec beta_dist_g1, beta_dist_g2;
    vec z_sample_in_g1, z_sample_in_g2;
    uvec z_sample_in_g1_index;
    
    double S_prob_num, S_prob_den, M_prob_num, M_prob_den, post_prob_num, post_prob_den;

    int vt;
    vec z_sample_tmp;
    z_sample_tmp = z_sample;
    
    Rcout << "GNaN_vec = " << GNaN_vec << endl;
    
    if(SPLIT){
        Rcout << "Try Splitting" << endl;
        G_star = G1+1;
        
        if(G1==1) bG = 1;
        if(G_star == Max_num_of_groups) dG = 1;
        
        Rcout << "Check 1" << endl;
        Rcout << "G1_vec =" << G1_vec << endl;
        
        if( GNaN_vec.is_empty()){
            r_samples.slice(iter+1) = r_sample;
            rho_samples.col(iter+1) = rho_sample;
            z_samples.col(iter+1) = z = z_sample;
            CompChg = false;
            return;
        }
        
        g1 = Rcpp::RcppArmadillo::sample(G1_vec, 1, true)(0);
        g2 = GNaN_vec(0);
        
        
        
        //cout << "Check 1" << endl;
        u = Rf_rbeta(2., 2.);
        
        Rcout << "u =" << u << endl;
   
        z_sample_in_g1_index = find(z_sample==g1);
        
        for(int i=0; i<z_sample_in_g1_index.n_elem; i++)
            z_sample_tmp(z_sample_in_g1_index(i)) = (Rf_runif(0., 1.)<u)? g1: g2;
       
        //Rcout << z_sample_tmp << endl;
        
        z_sample_in_g1 = z_sample_tmp(find(z_sample_tmp==g1));
        z_sample_in_g2 = z_sample_tmp(find(z_sample_tmp==g2));
        

        if((z_sample_in_g1.n_elem==0) || (z_sample_in_g2.n_elem==0)){
            r_samples.slice(iter+1) = r_sample;
            rho_samples.col(iter+1) = rho_sample;
            z_samples.col(iter+1) = z = z_sample;
            CompChg = false;
            return;
        }
        
        //cout << "Check 1" << endl;
        
        //r_sample_g1 =r_sample.row(g1);
        //r_sample_g2 =r_sample.row(g2);
        
        for(int p=0; p<Num_of_covariates; p++){
            if(r_sample(g1, p)==0)
                r_sample_g1(p) = r_sample_g2(p)=0;
            else{
                vt = Rcpp::RcppArmadillo::sample(sample_integers, 1, true)(0);
                if(vt==0){ r_sample_g1(p) = r_sample_g2(p)=1; }
                if(vt==1){ r_sample_g1(p) = 1; r_sample_g2(p)=0;}
                if(vt==2){ r_sample_g1(p) = 0; r_sample_g2(p)=1;}
                }
        }
        r_sample_tmp = r_sample;
        r_sample_tmp.row(g1) = r_sample_g1;
        r_sample_tmp.row(g2) = r_sample_g2;
        
        //cout <<"r_sample_tmp = " << r_sample_tmp << endl;
        beta_est = Calculate_Beta_Est(z_sample, r_sample_tmp, omega_sample);
        //cout << "beta_est = " << beta_est << endl;
        
        s1 = min_beta_distance_group_index(beta_est, g1);
        s2 = min_beta_distance_group_index(beta_est, g2);
        //cout << "Check 2" << endl;
        //cout << "g1 = " << g1 << "\t" << "g2 = " << g2 << endl;
        //cout << "s1 = " << s1 << "\t" << "s2 = " << s2 << endl;

        if( (z_sample_in_g1.n_elem>0) && (z_sample_in_g2.n_elem>0)){
            if( (s1 !=g2) && (s2 !=g1)){
                r_samples.slice(iter+1) = r_sample;
                rho_samples.col(iter+1) = rho_sample;
                z_samples.col(iter+1) = z = z_sample;
                CompChg = false;
                return;
            }
            if((s1 ==g2) && (s2 == g1))
                P_chos = 2./(G1+G0+1.);
            if((s1 ==g2) && (s2 != g1))
                P_chos = 1./(G1+G0+1.);
            if((s1 !=g2) && (s2 == g1))
                P_chos = 1./(G1+G0+1.);

        }
        else
            P_chos = 1./( (G1+G0+1.)*G1);
 
        rho_sample_tmp = rho_sample;
        rho_sample_tmp(g1) = rho_sample(g1)*u;
        rho_sample_tmp(g2) = rho_sample(g1)*(1-u);
        
        P_alloc = pow( 1./3., accu(r_sample.row(g1)) );
        
        
        S_prob_num = log(dG) + log(P_chos) + log(rho_sample(g1)) + Rf_dpois(G_star, lambda, 1);
        S_prob_den = log(bG) - log(G1) + (z_sample_in_g1.n_elem)*log(u) + (z_sample_in_g2.n_elem)*log(1-u) + Rf_dbeta(u, 2., 2., true)+ log (P_alloc) + Rf_dpois((G_star-1), lambda, 1);
        
        Group_Info_tmp = Group_Info;
        Group_Info_tmp(g1) = z_sample_in_g1.n_elem;
        Group_Info_tmp(g2) = z_sample_in_g2.n_elem;
        
        post_prob_num = Log_Posterior_Calculation(r_sample_tmp, z_sample_tmp, rho_sample_tmp, Group_Info_tmp, omega_sample);
        post_prob_den = Log_Posterior_Calculation(r_sample, z_sample, rho_sample, Group_Info, omega_sample);
        //cout <<"S_prob_num = " << S_prob_num << "\t" << "S_prob_den = " << S_prob_den << endl;
        
        //cout <<"post_prob_num = " << post_prob_num <<"\t" << "post_prob_den = " << post_prob_den << endl;
        
        //cout <<"Split (post_prob_num - post_prob_den + E_prob_num-E_prob_den) =" << (post_prob_num - post_prob_den + S_prob_num-S_prob_den) << endl;
        
        if(log(Rf_runif(0., 1.)) < (post_prob_num - post_prob_den + S_prob_num-S_prob_den) ){
            r_samples.slice(iter+1) = r_sample_tmp;
            rho_samples.col(iter+1) = rho_sample_tmp;
            z_samples.col(iter+1) = z = z_sample_tmp;
            
            //Group_Info.fill(datum::nan);
            //vec group_info_value = unique(z);
            //uvec group_info_num = hist(z,group_info_value);
            //Group_Info(conv_to< uvec >::from(group_info_value)) = conv_to< vec >::from(group_info_num);
            Group_Info = Group_Info_tmp;
            CompChg = true;
            Rcout << "A component was splitted" << endl;
        }
        else{
            r_samples.slice(iter+1) = r_sample;
            rho_samples.col(iter+1) = rho_sample;
            z_samples.col(iter+1) = z = z_sample;
            CompChg = false;
            /*
            vec group_info_value_check = unique(z);
            uvec group_info_num_check = hist(z,group_info_value_check);
            vec Group_Info_check;
            Group_Info_check.set_size(Max_num_of_groups);
            
            Group_Info_check(conv_to< uvec >::from(group_info_value_check)) = conv_to< vec >::from(group_info_num_check);

            Rcout << "Group_Info_check = " << endl << Group_Info_check << endl;
            */
            
            Rcout << "A component was NOT splitted" << endl;
        }
        
        //Rcout << "Update splitting" << endl;
    }
    else{
        //cout << "Try Merging" << endl;
        G_star = G-1;
        
        if(G==Max_num_of_groups) dG = 1;
        if(G_star == 1) bG = 1;
        
        g1 = Rcpp::RcppArmadillo::sample(G_vec, 1, true)(0);
        
        
        if(!any(z_sample==g1)){
            g2 = Rcpp::RcppArmadillo::sample(G1_vec, 1, true)(0);
            P_chos = 1./( (G1+G0)*G1);
            G1_star = G;
        }
        else{
            beta_est = Calculate_Beta_Est(z_sample, r_sample, omega_sample); //Calculate_Beta_Est(z_sample, r_sample_tmp)
            
            g2 = min_beta_distance_group_index(beta_est, g1);
            
            s2 = min_beta_distance_group_index(beta_est, g2);
            
            if(s2 !=g1)
                P_chos = 2./(G1+G0);
            if(s2 == g1)
                P_chos = 1./(G1+G0);
            G1_star = G1-1;

        }
        
        //Rcout << "g1 = " << g1 << "\t" << "g2 = " << g2 << endl;
        rho_sample_tmp = rho_sample;
        rho_sample_tmp(g1) = rho_sample(g1)+rho_sample(g2);
        rho_sample_tmp(g2) = 0;
        u = rho_sample(g1)/(rho_sample(g1)+rho_sample(g2));
        
        //Rcout << "u = " << u << endl;
        
        r_sample_tmp = r_sample;
        
        for(int i = 0; i<r_sample_tmp.n_cols; i++){
            r_sample_tmp(g1, i) = (r_sample(g1,i)==0 && r_sample(g2,i)==0)? 0 : 1;
        }
        
        //Rcout <<"r_sample_tmp.row(g1) = " << r_sample_tmp.row(g1) << endl;
        //Rcout <<"r_sample_tmp.row(g2) = " << r_sample_tmp.row(g2) << endl;

        //Rcout << "(Num_of_covariates-accu( (r_sample_tmp.row(g1)==0)* (r_sample_tmp.row(g2)==0))) = " << (Num_of_covariates-as_scalar( (r_sample_tmp.row(g1)==0)* (r_sample_tmp.row(g2)==0).t())) << endl;
        
        P_alloc = pow( 1./3., (Num_of_covariates-as_scalar( (r_sample_tmp.row(g1)==0)* (r_sample_tmp.row(g2)==0).t())) );

        z_sample_tmp = z_sample;
        
        z_sample_in_g1 = z_sample_tmp(find(z_sample_tmp==g1));
        z_sample_in_g2 = z_sample_tmp(find(z_sample_tmp==g2));

        
        z_sample_tmp(find(z_sample_tmp==g2)).fill(g1);
        
        //cout << "z_sample =" << z_sample << endl;
        r_sample_tmp.row(g2).fill(datum::nan);
        
        //Rcout << "z_sample = " << z_sample << endl;
        //Rcout << "u = " << u << "\t" << "(z_sample_in_g1.n_elem) = " << (z_sample_in_g1.n_elem) << "\t" << "(z_sample_in_g2.n_elem) = " << (z_sample_in_g2.n_elem) << endl;
        //Rcout << "Rf_dbeta(u, 100., 100., true) = " << Rf_dbeta(u, 2., 2., true) << endl;
        //Rcout << "log (P_alloc) = " << log (P_alloc) << endl;
        
        
        M_prob_num = log(bG) - log(G1_star) + (z_sample_in_g1.n_elem)*log(u) + (z_sample_in_g2.n_elem)*log(1-u)  + Rf_dbeta(u, 2., 2., true)+ log (P_alloc) + Rf_dpois(G_star, lambda, 1);
        M_prob_den = log(dG) + log(P_chos) + log(rho_sample(g1)+rho_sample(g2)) + Rf_dpois( (G_star+1), lambda, 1);
       
        //Rcout<<"G_star = " << G_star << "\tlambda = " << lambda << "\t" << "Rf_dpois(G_star, lambda, 1) = " << Rf_dpois(4, lambda, 1) << "\t" << "Rf_dpois( (G_star+1), lambda, 1) = " << Rf_dpois( (G_star+1), lambda, 1) << endl;
        
        Group_Info_tmp = Group_Info;
        Group_Info_tmp(g1) = z_sample_in_g1.n_elem + z_sample_in_g2.n_elem;
        Group_Info_tmp(g2) = datum::nan;
       
        
        //cout << "r_sample_tmp = " << r_sample_tmp << endl;
        //cout << "rho_sample_tmp = " << rho_sample_tmp << endl;
        
        post_prob_num = Log_Posterior_Calculation(r_sample_tmp, z_sample_tmp, rho_sample_tmp, Group_Info_tmp, omega_sample);
        post_prob_den = Log_Posterior_Calculation(r_sample, z_sample, rho_sample, Group_Info, omega_sample);
        
        //cout << "z_sample = " << unique(z_sample) << "\t" << "z = " << unique(z) << endl;
        //Rcout << "M_prob_num = " << M_prob_num << "\t" << "M_prob_den = " << M_prob_den << endl;
        //Rcout << "post_prob_num = " << post_prob_num << "\t" << "post_prob_den = " << post_prob_den << endl;

        //Rcout <<"Merge (post_prob_num - post_prob_den + E_prob_num-E_prob_den) =" << (post_prob_num - post_prob_den + M_prob_num-M_prob_den) << endl;
        
        if(log(Rf_runif(0., 1.)) < (post_prob_num - post_prob_den + M_prob_num-M_prob_den) ){
            r_samples.slice(iter+1) = r_sample_tmp;
            rho_samples.col(iter+1) = rho_sample_tmp;
            z_samples.col(iter+1) = z = z_sample_tmp;
            Group_Info = Group_Info_tmp;
            CompChg = true;
            Rcout << "Two components were merged." << endl;
        }
        else{
            r_samples.slice(iter+1) = r_sample;
            rho_samples.col(iter+1) = rho_sample;
            z_samples.col(iter+1) = z = z_sample;
            CompChg = false;
            Rcout << "Two components were not merged." << endl;
        }
        //Rcout << "Update merging." << endl;
    }
}

vec UnknownCompFMR::Update_weight(mat r_sample, vec z_sample, vec rho_sample, vec omega_sample)
{
    
    //Rcout << "Start update weight." << endl;
    double log_post_prob_old = 0, log_post_prob_new =0, SSE_g=0;
    
    uvec groups_index_tmp =  find(Group_Info>0); //find_finite(Group_Info); //unique(z_sample);
    
    colvec y_g;
    mat r_g;
    int n_g;
    
    
    mat Omega_g, Omega_g_inv, Weight_inv;
    vec omega_sample_tmp;
    uvec obs_g_index;
    mat x_g, Lambda_g, x_g_1;
    double h_g;
    
    double omega_new, omega_old;
    int g;
    //Rcout << "groups_index_tmp.n_elem = " << groups_index_tmp.n_elem << endl;
    for(int g_tmp =0; g_tmp<groups_index_tmp.n_elem; g_tmp++){
        g = groups_index_tmp(g_tmp);
        obs_g_index = find(z_sample == g);
        y_g = Y.elem( obs_g_index );
        
        //Rcout << "obs_g_index = " << obs_g_index << endl;
        //Rcout <<"y_g.n_elem 1= " << y_g.n_elem << endl;
        if(y_g.n_elem!=0){
            //Rcout <<"y_g.n_elem 2 = " << y_g.n_elem << endl;
            r_g = r_sample.row(g);
            n_g = accu(z_sample==g);
            h_g = 0.5*(n_g+sigma2_a);
            x_g = X.rows(obs_g_index);
            x_g_1 = x_g.cols(find(r_g>0));
            omega_sample_tmp = omega_sample(obs_g_index);
            for(int obs_index = 0; obs_index<obs_g_index.n_elem; obs_index++){
                omega_old = omega_sample_tmp(obs_index);
                Weight_inv = diagmat(1./omega_sample_tmp);
                
                if(accu(r_g)>0){
                    
                    Lambda_g = x_g_1.t()*Weight_inv*x_g_1 + eye(accu(r_g),accu(r_g))/c2;
                    
                    SSE_g = as_scalar(y_g.t()*(Weight_inv - Weight_inv*x_g_1*Lambda_g.i()*x_g_1.t()*Weight_inv)*y_g);
                    
                    log_post_prob_old = -(0.5*(nu+1.)+1.)*log(omega_old) - 0.5*(nu/omega_old) - 0.5*log(det(Lambda_g))-h_g*log(0.5*(sigma2_b+SSE_g));
                }
                else{
                    
                    SSE_g = as_scalar(y_g.t()*(Weight_inv)*y_g);
                    
                    log_post_prob_old =  -(0.5*(nu+1.)+1.)*log(omega_old) - 0.5*(nu/omega_old) -h_g*log(0.5*(sigma2_b+SSE_g));
                    
                }
                do {
                    omega_new = Rf_rnorm(omega_old, tuning_para);
                } while (omega_new<0);
                
                omega_sample_tmp(obs_index) = omega_new;
                Weight_inv = diagmat(1./omega_sample_tmp);
                if(accu(r_g)>0){
                    
                    Lambda_g = x_g_1.t()*Weight_inv*x_g_1 + eye(accu(r_g),accu(r_g))/c2;
                    
                    SSE_g = as_scalar(y_g.t()*(Weight_inv - Weight_inv*x_g_1*Lambda_g.i()*x_g_1.t()*Weight_inv)*y_g);
                
                    log_post_prob_new = -(0.5*(nu+1.)+1.)*log(omega_new) - 0.5*(nu/omega_new) -0.5*log(det(Lambda_g))-h_g*log(0.5*(sigma2_b+SSE_g));
                }
                else{
                    
                    SSE_g = as_scalar(y_g.t()*(Weight_inv)*y_g);
                    
                    log_post_prob_new = -(0.5*(nu+1.)+1.)*log(omega_new) - 0.5*(nu/omega_new)-h_g*log(0.5*(sigma2_b+SSE_g));

                }
                //Rcout << "log_post_prob_new-log_post_prob_old = " << log_post_prob_new-log_post_prob_old << endl;
                if( (log_post_prob_new-log_post_prob_old) < log(Rf_runif(0., 1.)) ){
                    omega_sample_tmp(obs_index) = omega_old;
                    acc_rate(obs_g_index(obs_index))--;
                   
                }
            }
        }
        omega_sample(obs_g_index) = omega_sample_tmp;
    }
    //Rcout << "End update weight." << endl;
    return omega_sample;
}

/*
SEXP UnknownCompFMR::RemoveLabelSwitching(List PosteriorSamplesForLabelSwitching)
{
    mat z_samples_LS;
    cube r_samples_LS;
    mat rho_samples_LS;
    
    int num_of_efficient_iterations, est_num_of_groups;
    
    mat beta_est;
    colvec y_g;
    mat r_g;
    mat x_g, x_g_1, Lambda_g;
    double est_mean, est_variance;
    int q_g;
    double alpha_tmp, beta_tmp;
    
    z_samples_LS = as<mat>(PosteriorSamplesForLabelSwitching["rho.samples"]);
    r_samples_LS = as<cube>(PosteriorSamplesForLabelSwitching["r.samples"]);
    rho_samples_LS = as<mat>(PosteriorSamplesForLabelSwitching["rho.samples"]);
    
    num_of_efficient_iterations = rho_samples_LS.n_cols;
    est_num_of_groups = rho_samples_LS.n_rows;
    
    cube prob_cube(num_of_efficient_iterations, Num_of_obs,  est_num_of_groups);
    
    for(int m=0; m<num_of_efficient_iterations; m++)
        for(int g=0; g<est_num_of_groups; g++){
            y_g = Y.elem( find(z_samples_LS.col(m) ==g) );
            x_g = X.rows(find(z_samples_LS.col(m) ==g));
            r_g = r_samples_LS.slice(m).row(g);
            q_g = accu(r_g);
            x_g_1 = x_g.cols(find(r_g>0));
            Lambda_g = x_g_1.t()*x_g_1 + eye(accu(r_g),accu(r_g))/c2;
            beta_est = Lambda_g.i()*x_g_1.t()*y_g;
            est_mean = as_scalar(x_g_1.t() * beta_est);
            alpha_tmp =0.5*( y_g.n_elem + q_g + sigma2_a);
            beta_tmp = 0.5*(as_scalar((y_g - x_g_1.t() * beta_est).t()*(y_g - x_g_1.t() * beta_est) + beta_est.t()*Lambda_g.i()*beta_est) + sigma2_b);
            est_variance = beta_tmp/(alpha_tmp-1);
            for(int i=0; i<Num_of_obs; i++){
                prob_cube(m, i, g)= Rf_dnorm4(Y(i), est_mean, est_variance, 0);
            }
        }
    
}
*/
SEXP UnknownCompFMR::MCMC_Procedure()
{
    Rcout << "============= FMR: MCMC =============="<< endl;
    
    List PosteriorSamples;
    
    mat r_sample = r_samples.slice(0);
    vec rho_sample = rho_samples.col(0);
    vec z_sample = z;
    vec omega_sample = omega_samples.col(0);
    vec group_info_value = unique(z);
    uvec group_info_num = hist(z,group_info_value);
    //Group_Info(conv_to< uvec >::from(group_info_value)) = conv_to< vec >::from(group_info_num);

    //vec Group_Info_tmp;
    acc_rate.set_size(Num_of_obs);
    acc_rate.fill(Num_of_iterations);

    
    int iter = 0;
    while(iter < Num_of_iterations-1){
        
        Rcout << "iter = " << iter << endl;
        
        group_info_value = unique(z);
        group_info_num = hist(z,group_info_value);
        Group_Info(conv_to< uvec >::from(group_info_value)) = conv_to< vec >::from(group_info_num);

        if (CompChg) {
            for(int iter_inner=0; iter_inner < Num_of_iterations_Inside; iter_inner++){
                r_sample = Updata_gamma(z_sample, r_sample, omega_sample);
                //Rcout << "r_sample" << endl;
                z_sample = Update_Z(z_sample, rho_sample, r_sample, omega_sample);
                //Rcout << "z_sample" << endl;
                rho_sample = Updata_rho();
                //Rcout << "rho_sample" << endl;
            }

        }
        else{
            r_sample = Updata_gamma(z_sample, r_sample, omega_sample);
            z_sample = Update_Z(z_sample, rho_sample, r_sample, omega_sample);
            rho_sample = Updata_rho();
        }
        

        
        Rcout << "Check #1" << endl;
        
        z_samples.col(iter+1) = z = z_sample;
        r_samples.slice(iter+1) = r_sample;
        rho_samples.col(iter+1) = rho_sample;
        Rcout << "Check #2" << endl;
        if(Robust){
            omega_sample = Update_weight(r_sample, z_sample, rho_sample, omega_sample);
            omega_samples.col(iter+1) = omega_sample;
            //Rcout << "Rboust" << endl;
        }
        else{
            omega_sample.ones();
            omega_samples.col(iter+1).ones();
        }
        Rcout << "Check #3" << endl;

        group_info_value = unique(z);
        group_info_num = hist(z,group_info_value);
        Group_Info(conv_to< uvec >::from(group_info_value)) = conv_to< vec >::from(group_info_num);
        Rcout << "Check #4" << endl;
        if(RJMCMC){
            Rcout << "Before Generating_Eliminating " << endl;
            Generating_Eliminating(iter, z_sample, r_sample, rho_sample, omega_sample);
            Rcout << "After Generating_Eliminating " << endl;
            
            z_sample = z = z_samples.col(iter+1);
            r_sample = r_samples.slice(iter+1);
            rho_sample = rho_samples.col(iter+1);
            
            //Rcout << Group_Info << endl;
            //Rcout << "===========================" << endl;
            //Rcout << G1_vec << endl;
            
            Rcout << "Before Splitting_Merging " << endl;
            Splitting_Merging(iter, z_sample, r_sample, rho_sample, omega_sample);
            Rcout << "After Splitting_Merging " << endl;
            
            z_sample = z = z_samples.col(iter+1);
            r_sample = r_samples.slice(iter+1);
            rho_sample = rho_samples.col(iter+1);
            
            //Rcout << Group_Info << endl;
            //Rcout << "===========================" << endl;
            //Rcout << G1_vec << endl;
        }
        
        Rcout << "Group_Info = \n" << Group_Info << endl;
        
        iter++;
    }
    
    
    PosteriorSamples["z"] = z_samples;
    PosteriorSamples["rho.samples"] = rho_samples;
    PosteriorSamples["r.samples"] = r_samples;
    PosteriorSamples["omega.samples"] = omega_samples;
    PosteriorSamples["acceptance.rate"] = acc_rate/Num_of_iterations;
    
    Rcout << "MCMC_Procedure" << endl;
    
    return (PosteriorSamples);
}

