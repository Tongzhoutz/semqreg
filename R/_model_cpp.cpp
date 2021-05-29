////////////   Data created: This file is to do the outer and inner loops for my model
// [[Rcpp::plugins("cpp11")]]

// [[Rcpp::depends(RcppArmadillo,RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppArmadillo.h>
#include <RcppNumerical.h>
using namespace Numer;


typedef Eigen::Map<Eigen::MatrixXd> MapMat;
typedef Eigen::Map<Eigen::VectorXd> MapVec;

class QuantileReg: public MFuncGrad
{
private:
  const MapMat X;
  const MapVec Y;
  const double tau;
public:
  QuantileReg(const MapMat x_, const MapVec y_, const double tau_) : X(x_), Y(y_),tau(tau_) {}
  
  double f_grad(Constvec& beta, Refvec grad)
  {
    // Negative log likelihood
    //   sum(log(1 + exp(X * beta))) - y' * X * beta
    
    Eigen::VectorXd xbeta = X * beta;
    //Eigen::VectorXd bool_v = (Y.array() < xbeta.array()).cast<double>().array();
    Eigen::VectorXd vv = tau -  (Y.array() < xbeta.array()).cast<double>().array();
    Eigen::VectorXd vvv = Y.array() - xbeta.array();
    Eigen::VectorXd v4 = vvv.array()*vv.array();
    const double f = v4.sum()/Y.rows();
    
    grad.noalias() = (-X.transpose() * vv)/Y.rows();
    
    return f;
  }
};

// [[Rcpp::export]]
Rcpp::NumericVector q_reg(Rcpp::NumericMatrix x, Rcpp::NumericVector y,double tau)
{
  const MapMat xx = Rcpp::as<MapMat>(x);
  const MapVec yy = Rcpp::as<MapVec>(y);
  
  // Negative log likelihood
  QuantileReg nll(xx, yy,tau);
  // Initial guess
  Eigen::VectorXd beta(xx.cols());
  beta.setRandom();
  
  double fopt;
  int maxite = 10000000;
    double eps_f = 1e-8;
    double eps_g = 1e-8;
  int status = optim_lbfgs(nll, beta, fopt,maxite,eps_f,eps_g);
  if(status < 0)
    Rcpp::stop("fail to converge");
  
  return Rcpp::wrap(beta);
}

//[[Rcpp::export]]
arma::vec if_elses(arma::vec min, arma::vec max, arma::colvec V1, arma::colvec V2){
  return V1 % (min <= max) + V2 % (min > max);
}


// [[Rcpp::export]]
arma::mat  posterior_cpp2(arma::mat Mat_draw,    // input: U
                         arma::vec utoty_nor,    
                         arma::vec uc_resid,
                         arma::mat Herm_age,
                         arma::mat regressor_V1,arma::mat regressor_U1,
                         arma::mat Mat_init_Vt, arma::mat Mat_init_V1,arma::mat Mat_init_Ut, arma::mat Mat_init_U1,
                         double b1_ut,double bL_ut,double b1_u1,double bL_u1,double b1_vt, double bL_vt,double b1_v1,double bL_v1,double a,double b, double c, double d, double e, arma::vec tauList, double meanY, double stdY
                       ) {
  double tau_diff = tauList(1) - tauList(0);
  int  N = 792;
  int T = 6;
  int L = tauList.n_elem; 
  //  conditions 
  arma::vec condition_1;   // the first interval
  arma::vec condition_Ntau;  // the last interval

  // vectorization of Mat_draw
 arma::vec Mat_vec_15 = arma::vectorise( Mat_draw.cols(0,4)  );
 arma::vec Mat_vec_26 = arma::vectorise( Mat_draw.cols(1,5) );
  
  // Hermite polynomials for age variable in DT
     
    // Density of V_t | V_{t-1} 
    int  k3 = 2; 
     // Vectorize 
  arma::vec M_vec = arma::vectorise(Mat_draw); 
  arma::vec Vect = utoty_nor - M_vec;   // y - u 
  arma::mat Herm_Vect = arma::ones( N*T, (k3+1)); // normalize Vect
  Herm_Vect.col(1) = (Vect-meanY)/stdY;
  Herm_Vect.col(2) = pow( (Vect-meanY)/stdY,2.0) - arma::ones(N*T);

            // regressors
  arma::mat regressor_Vt_updated( N*(T-1), 2*(k3+1), arma::fill::zeros  ); // k3 = 2
  for (int i = 0 ; i < (k3+1)  ; ++i){
    for (int j = 0; j < 2 ; ++j){
            regressor_Vt_updated.col( 2*i + j )= Herm_Vect.col(i).head(N*(T-1)) % Herm_age.col(j).tail(N*(T-1));
    }
  }
             // prepare 
  arma::mat predict_vt_updated = regressor_Vt_updated * Mat_init_Vt;
  arma::mat mat_diff_vt = arma::diff(predict_vt_updated , 1, 1);
// if (  mat_diff_vt.min() < 0  ){ Rcpp::stop("Negative mat_diff_vt"); }
    arma::umat vt_infinite = find(mat_diff_vt == 0);
    mat_diff_vt.elem(vt_infinite).fill(10000000);
  arma::vec dens2_vt( N*(T-1), arma::fill::zeros  );
  arma::vec Vect_dep = Vect.tail(N*(T-1));
  arma::vec condition;

  for ( int i = 0 ; i < L-1 ; ++i ){
            
            condition = arma::conv_to<arma::vec>::from( Vect_dep > predict_vt_updated.col(i) && Vect_dep <= predict_vt_updated.col(i+1)   );

           dens2_vt = dens2_vt +  condition % (tau_diff/mat_diff_vt.col(i)); 
  } 

       // tail intervals: v_t

       condition_1 = arma::conv_to<arma::vec>::from( Vect_dep <= predict_vt_updated.col(0)  ); 
       
       condition_Ntau = arma::conv_to<arma::vec>::from( Vect_dep > predict_vt_updated.col(10)  ); 

       dens2_vt = dens2_vt + tauList(0) * b1_vt * arma::trunc_exp( b1_vt* (Vect_dep - predict_vt_updated.col(0)  )) % condition_1 + (1 - tauList(10)) * bL_vt * arma::trunc_exp( -bL_vt * (Vect_dep - predict_vt_updated.col(10))  ) % condition_Ntau;
    

       // v_1

             // prepare 
  arma::mat predict_v1_updated = regressor_V1 * Mat_init_V1;
  arma::mat mat_diff_v1 = arma::diff(predict_v1_updated , 1, 1);
 //if (  mat_diff_v1.min() < 0  ){ Rcpp::stop("Negative mat_diff_v1"); }
    arma::umat v1_infinite = find(mat_diff_v1 == 0);
    mat_diff_v1.elem(v1_infinite).fill(10000000000);
  arma::vec dens2_v1( N, arma::fill::zeros  );
  arma::vec Vect_dep_v1 = Vect.head(N);


  for (int i= 0 ; i < L-1 ; ++i){

      condition = arma::conv_to<arma::vec>::from( Vect_dep_v1 > predict_v1_updated.col(i) && Vect_dep_v1 <= predict_v1_updated.col(i+1)   );

           dens2_v1 = dens2_v1 +  condition % (tau_diff/mat_diff_v1.col(i)); 

  }

       // tail intervals: v_1

      condition_1 = arma::conv_to<arma::vec>::from( Vect_dep_v1 <= predict_v1_updated.col(0)  ); 
       
       condition_Ntau = arma::conv_to<arma::vec>::from( Vect_dep_v1 > predict_v1_updated.col(10)  ); 

       dens2_v1 = dens2_v1 + tauList(0) * b1_v1 * arma::trunc_exp( b1_v1* (Vect_dep_v1 - predict_v1_updated.col(0)  )) % condition_1 + (1 - tauList(10)) * bL_v1 * arma::trunc_exp( -bL_v1 * (Vect_dep_v1 - predict_v1_updated.col(10))  ) % condition_Ntau;
    


    // Prior: U_t | U_{t-1}, AGE_t

            // Prepare
        arma::mat regressor_ut_updated( N*(T-1), 12, arma::fill::zeros  );
        arma::mat Herm_ut(N*(T-1),4,arma::fill::zeros);
        Herm_ut.col(0) = arma::ones(N*(T-1)); 
        Herm_ut.col(1) = (Mat_vec_15 - meanY )/stdY;
        Herm_ut.col(2) = arma::pow( (Mat_vec_15-meanY)/stdY   , 2.0   ) - arma::ones(N*(T-1));
        Herm_ut.col(3) = arma::pow( (Mat_vec_15-meanY)/stdY , 3.0  ) - 3 * (Mat_vec_15-meanY)/stdY;

        for (int i = 0 ; i < 4 ; ++i){
          for (int j = 0 ; j < 3 ; ++j){

             regressor_ut_updated.col( 3*i + j ) = Herm_ut.col(i) % Herm_age.col(j).tail(N*(T-1));

          }
        }

        arma::mat predict_ut_updated = regressor_ut_updated * Mat_init_Ut;
        arma::mat mat_diff_ut = arma::diff(predict_ut_updated , 1, 1);
 //if (  mat_diff_ut.min() < 0  ){ Rcpp::stop("Negative mat_diff_ut"); }
    arma::umat ut_infinite = find(mat_diff_ut == 0);
    mat_diff_ut.elem(ut_infinite).fill(1000000);
        arma::vec dens1_ut(N*(T-1), arma::fill::zeros);

  for ( int i = 0 ; i < L-1 ; ++i ){
            
            condition = arma::conv_to<arma::vec>::from( Mat_vec_26 > predict_ut_updated.col(i) && Mat_vec_26 <= predict_ut_updated.col(i+1)   );

           dens1_ut = dens1_ut +  condition % (tau_diff/mat_diff_ut.col(i)); 
  } 

             // tail : U_t

       condition_1 = arma::conv_to<arma::vec>::from( Mat_vec_26 <= predict_ut_updated.col(0)  ); 
       
       condition_Ntau = arma::conv_to<arma::vec>::from( Mat_vec_26 > predict_ut_updated.col(10)  ); 

       dens1_ut = dens1_ut + tauList(0) * b1_ut * arma::trunc_exp( b1_ut* (Mat_vec_26 - predict_ut_updated.col(0)  )) % condition_1 + (1 - tauList(10)) * bL_ut * arma::trunc_exp( -bL_ut * (Mat_vec_26 - predict_ut_updated.col(10))  ) % condition_Ntau;
    
     //  arma::umat exp_infinite = find( Mat_vec_26 == 0);
      // mat_diff_ut.elem(ut_infinite).fill(1000000);
    


       // Prior: U_1
    arma::mat predict_u1_updated = regressor_U1 * Mat_init_U1;
    arma::mat mat_diff_u1 = diff(predict_u1_updated,1,1);
 //if (  mat_diff_u1.min() < 0  ){ Rcpp::stop("Negative mat_diff_u1"); }
    arma::umat u1_infinite = find(mat_diff_u1 == 0);
    mat_diff_u1.elem(u1_infinite).fill(100000000000);
    arma::vec dens1_u1(N,arma::fill::zeros);
    arma::vec M_1 = Mat_draw.col(0);

   for (int i= 0 ; i < L-1 ; ++i){

      condition = arma::conv_to<arma::vec>::from( M_1  > predict_u1_updated.col(i) && M_1 <= predict_u1_updated.col(i+1)   );

           dens1_u1 = dens1_u1 +  condition % (tau_diff/mat_diff_u1.col(i)); 

  }

       // tail intervals: U_1

      condition_1 = arma::conv_to<arma::vec>::from( M_1 <= predict_u1_updated.col(0)  ); 
       
       condition_Ntau = arma::conv_to<arma::vec>::from( M_1 > predict_u1_updated.col(10)  ); 

       dens1_u1 = dens1_u1 + tauList(0) * b1_u1 * arma::trunc_exp( b1_u1* (M_1 - predict_u1_updated.col(0)  )) % condition_1 + (1 - tauList(10)) * bL_u1 * arma::trunc_exp( -bL_u1 * (M_1 - predict_u1_updated.col(10))  ) % condition_Ntau;
    

       // Prior: consumption
       arma::vec dens_cmp =d + e/(1 + arma::trunc_exp( b * M_vec + c*Vect - a) );
       //arma::vec dens_cmp(N*T,arma::fill::ones);,
       arma::vec dens_tot = dens1_u1 % dens2_v1 % dens_cmp.tail(N);
       
       for ( int t = 0 ; t < 5 ; ++t ){
              
         dens_tot = dens_tot % dens2_vt.subvec( t*N , t*N+N-1) % dens1_ut.subvec(t*N,t*N+N-1) % dens_cmp.subvec(t*N, t*N+N-1) ;

       }
      // arma::vec dens_ret(N,arma::fill::ones);
      // for ( int t = 0 ; t < 5 ; ++t ){
       //dens_ret_vt = dens_ret % dens2_vt.subvec( t*N , t*N+N-1);
     //  }
     //  arma::vec dens_ret_ut(N,arma::fill::ones);
       
      // for ( int t = 1 ; t < 2 ; ++t ){
       //dens_ret_ut = dens_ret_ut % dens1_ut.subvec( t*N , t*N+N-1);
     //  }
    return dens_tot; 
}

//[[Rcpp::export]]

Rcpp::List outerLoop(int maxiter, int draws, int N, int T,
                     arma::vec tauList,
                     double var_prop1, double var_prop2, double var_prop3,
                     double var_prop4, double var_prop5, double var_prop6,
                     arma::vec utoty_nor,arma::mat init, arma::vec uc_resid,
                     arma::mat Herm_age,arma::mat regressor_V1,
                     arma::mat regressor_U1, arma::mat Mat_init_Vt, arma::mat Mat_init_V1, arma::mat Mat_init_Ut, arma::mat Mat_init_U1,double b1_ut, double bL_ut,
                     double b1_u1, double bL_u1, double b1_vt, double bL_vt, double b1_v1, double bL_v1, double a, double b,double c,double d, double e, double meanY, double stdY){

 // Initialize matrices: Nu_chain, Obj_chain, acc, acceptrate

  arma::vec onesVec= arma::ones(N);
  arma::vec zerosVec = arma::zeros(N);
  arma::mat Matdraws = init;
  arma::vec likelihood(maxiter, arma::fill::zeros);
  arma::vec logg;

  // cubes,matrices and vectors used in the E-step
  arma::cube Nu_chain(N,draws,T,arma::fill::zeros);
  arma::cube acc(N,draws,T,arma::fill::zeros);
  // matrices and vectors in the M-step:
  arma::vec Vt_2T(N*(T-1),arma::fill::zeros);
  arma::vec U1(N,arma::fill::zeros); 
  arma::vec V1(N,arma::fill::zeros); 
  arma::vec Vt_15(N*(T-1),arma::fill::zeros);
  arma::mat Udraw_lag(N*(T-1),12,arma::fill::zeros); 
  arma::mat Vdraw_lag(N*(T-1),6,arma::fill::zeros);
  arma::mat Herm_ut_regressor(N*(T-1),4,arma::fill::ones);
  arma::mat Herm_vt_regressor(N*(T-1),3,arma::fill::ones);

  arma::vec Ut_2T(N*(T-1),arma::fill::zeros);
  arma::vec Ut_lag(N*(T-1),arma::fill::zeros);
  arma::mat coef_Ut(12,11,arma::fill::zeros);
  arma::mat coef_U1(4,11,arma::fill::zeros);
  arma::mat coef_Vt(6,11,arma::fill::zeros);
  arma::mat coef_V1(3,11,arma::fill::zeros);
  arma::cube Resqnew_Ut(12,11,maxiter, arma::fill::zeros);
  arma::cube Resqnew_U1(4,11,maxiter, arma::fill::zeros);
  arma::cube Resqnew_Vt(6,11,maxiter, arma::fill::zeros);
  arma::cube Resqnew_V1(3,11,maxiter, arma::fill::zeros);

  // updates for b1_ and bL_
  arma::vec Vect_b1_ut(N*(T-1),arma::fill::zeros);
  arma::vec Vect_bL_ut(N*(T-1),arma::fill::zeros);
  arma::vec Vect_b1_u1(N,arma::fill::zeros);
  arma::vec Vect_bL_u1(N,arma::fill::zeros);
  arma::vec Vect_b1_vt(N*(T-1),arma::fill::zeros);
  arma::vec Vect_bL_vt(N*(T-1),arma::fill::zeros);
  arma::vec Vect_b1_v1(N,arma::fill::zeros);
  arma::vec Vect_bL_v1(N,arma::fill::zeros);
  arma::uvec ids_utb1_num ; 
  arma::uvec ids_utbL_num ; 
  arma::uvec ids_u1b1_num ; 
  arma::uvec ids_u1bL_num ; 
  arma::uvec ids_vtb1_num ; 
  arma::uvec ids_vtbL_num ; 
  arma::uvec ids_v1b1_num ; 
  arma::uvec ids_v1bL_num ; 
  arma::vec zeros_15u_b1 = arma::zeros(N*(T-1));
  arma::vec zeros_15u_bL = arma::zeros(N*(T-1));
  arma::vec zeros_1u_b1 = arma::zeros(N);
  arma::vec zeros_1u_bL = arma::zeros(N);
  arma::vec zeros_15v_b1 = arma::zeros(N*(T-1));
  arma::vec zeros_15v_bL = arma::zeros(N*(T-1));
  arma::vec zeros_1v_b1 = arma::zeros(N);
  arma::vec zeros_1v_bL = arma::zeros(N);
  arma::vec denom_15_b1(N*(T-1),arma::fill::zeros);
  arma::vec denom_15_bL(N*(T-1),arma::fill::zeros);
  arma::vec denom_1_b1(N,arma::fill::zeros);
  arma::vec denom_1_bL(N,arma::fill::zeros);
  for (int t = 0 ; t < T ; ++t){
    Nu_chain.slice(t) = Matdraws.col(t) * arma::ones(1,draws);
  }

  arma::mat Obj_chain(N,draws,arma::fill::zeros); 
  Obj_chain.col(0) =  posterior_cpp2(Matdraws,utoty_nor,uc_resid,Herm_age,regressor_V1,regressor_U1,Mat_init_Vt,Mat_init_V1,Mat_init_Ut,Mat_init_U1,b1_ut,bL_ut,b1_u1,bL_u1,b1_vt,bL_vt,b1_v1,bL_v1,a,b,c,d,e,tauList,meanY,stdY);

  arma::mat acceptrate(draws,T,arma::fill::zeros);
  arma::mat mean_accept(maxiter,T,arma::fill::zeros);
  arma::mat last_Ut;
  arma::mat mat_b(maxiter,8,arma::fill::zeros);

  // matrices and vectors used in the inner loop: Metropolis-Hastings Algorithm
  arma::vec sd_mh = {sqrt(var_prop1),sqrt(var_prop2),sqrt(var_prop3),sqrt(var_prop4),
  sqrt(var_prop5), sqrt(var_prop6)};
  arma::vec newObj(N,arma::fill::zeros);
  arma::vec r(N,arma::fill::zeros);
  arma::vec prob(N,arma::fill::zeros);

  for (int iter = 0 ; iter < maxiter ; ++iter){
    Rcpp::Rcout << iter << "\n" ; 
    for (int j  = 1 ; j < draws ; ++j){
             Matdraws.zeros(); 
      for (int t = 0 ; t < T; ++t){
             Matdraws.col(t) = Nu_chain.slice(t).col(j-1);
      }
          
          // U_1

          Matdraws.col(0) = Nu_chain.slice(0).col(j-1) + sd_mh(0) * arma::randn(N);
          newObj = posterior_cpp2(Matdraws,utoty_nor,uc_resid,Herm_age,regressor_V1,regressor_U1,Mat_init_Vt,Mat_init_V1,Mat_init_Ut,Mat_init_U1,b1_ut,bL_ut,b1_u1,bL_u1,b1_vt,bL_vt,b1_v1,bL_v1,a,b,c,d,e,tauList,meanY,stdY);
          r = arma::min(onesVec,newObj/Obj_chain.col(j-1));
          prob = arma::randu(N);
          acc.slice(0).col(j) = arma::conv_to<arma::vec>::from( prob <= r );
          Obj_chain.col(j) = acc.slice(0).col(j) % newObj + ( onesVec - acc.slice(0).col(j)  ) % Obj_chain.col(j-1);
          Nu_chain.slice(0).col(j) = acc.slice(0).col(j) % Matdraws.col(0) + (onesVec - acc.slice(0).col(j)) % Nu_chain.slice(0).col(j-1);
          Matdraws.col(0) = Nu_chain.slice(0).col(j);


          // U_2 - U_6

          for ( int t = 1 ; t < T ; ++t ){
               Matdraws.col(t) = Nu_chain.slice(t).col(j-1) + sd_mh(t) * arma::randn(N);
               newObj = posterior_cpp2(Matdraws,utoty_nor,uc_resid,Herm_age,regressor_V1,regressor_U1,Mat_init_Vt,Mat_init_V1,Mat_init_Ut,Mat_init_U1,b1_ut,bL_ut,b1_u1,bL_u1,b1_vt,bL_vt,b1_v1,bL_v1,a,b,c,d,e,tauList,meanY,stdY);
          r = arma::min(onesVec,newObj/Obj_chain.col(j));
          prob = arma::randu(N);
          acc.slice(t).col(j) = arma::conv_to<arma::vec>::from( prob <= r );
          Obj_chain.col(j) = acc.slice(t).col(j) % newObj + ( onesVec - acc.slice(t).col(j) )% Obj_chain.col(j);
          Nu_chain.slice(t).col(j) = acc.slice(t).col(j) % Matdraws.col(t) + (onesVec - acc.slice(t).col(j)) % Nu_chain.slice(t).col(j-1);
          Matdraws.col(t) = Nu_chain.slice(t).col(j);
          } 
          
          for ( int k = 0 ; k < T ; ++k ){
              acceptrate( j-1 ,k ) = arma::mean( acc.slice(k).col(j-1)  );

          }

    } 
    for (int t = 0 ; t < T; ++t){
         mean_accept(iter,t) = arma::mean(acceptrate.col(t) );
      }
    Rcpp::Rcout <<  mean_accept(iter,0) << "\n" ; 
    Rcpp::Rcout <<  mean_accept(iter,1) << "\n" ; 
    Rcpp::Rcout <<  mean_accept(iter,2) << "\n" ; 
    Rcpp::Rcout <<  mean_accept(iter,3) << "\n" ; 
    Rcpp::Rcout <<  mean_accept(iter,4) << "\n" ; 
    Rcpp::Rcout <<  mean_accept(iter,5) << "\n" ; 
// E-step is over 
      
          // Keep the last draw 
    arma::mat last_Ut1 = join_rows(Nu_chain.slice(0).col(draws-1),Nu_chain.slice(1).col(draws-1),Nu_chain.slice(2).col(draws-1),Nu_chain.slice(3).col(draws-1));

    last_Ut = join_rows(last_Ut1, Nu_chain.slice(4).col(draws-1),Nu_chain.slice(5).col(draws-1));
//  M-Step

         //  U_t | U_{t-1}

     Ut_2T = arma::vectorise( last_Ut.cols(1,5)  ); //dependent variable
     Ut_lag = arma::vectorise(last_Ut.cols(0,4));

    Herm_ut_regressor.col(1) = ( Ut_lag - meanY ) / stdY;
    Herm_ut_regressor.col(2) = arma::pow( (Ut_lag-meanY)/stdY   , 2.0 ) - arma::ones(N*(T-1));
    Herm_ut_regressor.col(3) = arma::pow( (Ut_lag-meanY)/stdY   , 3.0 ) - 3*(Ut_lag - meanY  )/stdY;

    for (int k1 = 0 ; k1 < 4 ; ++k1){
      for (int k2 = 0 ; k2 < 3; ++k2){
            Udraw_lag.col(3*k1 + k2) = Herm_ut_regressor.col(k1) % Herm_age.col(k2).rows((N),6*N-1);
      }
    }

      //  U_1

     U1 = last_Ut.col(0);

    // V_t 
    Vt_2T = utoty_nor.tail(N*(T-1)) - vectorise(last_Ut.cols(1,5));
    Vt_15 = utoty_nor.head(N*(T-1)) - vectorise(last_Ut.cols(0,4));

    Herm_vt_regressor.col(1) = (Vt_15-meanY)/stdY;
    Herm_vt_regressor.col(2) = arma::pow((Vt_15-meanY)/stdY,2.0) - arma::ones(N*(T-1));

    for (int k1 = 0 ; k1 < 3 ; ++k1){
      for (int k2 = 0 ; k2 < 2; ++k2){
            Vdraw_lag.col(2*k1 + k2) = Herm_vt_regressor.col(k1) % Herm_age.col(k2).rows((N),6*N-1);
      }
    }
    // V_1
    V1 = utoty_nor.head(N) - last_Ut.col(0);

    for (int i = 0 ; i < 11 ; ++i){
     Resqnew_Ut.slice(iter).col(i) = Rcpp::as<arma::vec>(q_reg(Rcpp::wrap(Udraw_lag),Rcpp::wrap(Ut_2T),tauList(i)));
     Resqnew_U1.slice(iter).col(i) = Rcpp::as<arma::vec>(q_reg(Rcpp::wrap(regressor_U1),Rcpp::wrap(U1),tauList(i)));
     Resqnew_Vt.slice(iter).col(i) = Rcpp::as<arma::vec>(q_reg(Rcpp::wrap(Vdraw_lag),Rcpp::wrap(Vt_2T),tauList(i)));
     Resqnew_V1.slice(iter).col(i) = Rcpp::as<arma::vec>(q_reg(Rcpp::wrap(regressor_V1),Rcpp::wrap(V1),tauList(i)));
    }

    // Normalization
   
     arma::vec norm_Ut = mean(Resqnew_Ut.slice(iter),1);
            Resqnew_Ut.slice(iter) = Resqnew_Ut.slice(iter) - norm_Ut*arma::ones(1,11);
            Resqnew_Ut.slice(iter).row(0) =Resqnew_Ut.slice(iter).row(0) - ((1-tauList.tail(1))/bL_ut - tauList.head(1)/b1_ut )*arma::ones(1,11);
           Resqnew_Ut.slice(iter).row(3) = Resqnew_Ut.slice(iter).row(3)  + stdY;
            //Resqnew_Ut.slice(iter).row(3) = Resqnew_Ut.slice(iter).row(3) + 1.5;
     arma::vec norm_Vt = mean(Resqnew_Vt.slice(iter),1);
            Resqnew_Vt.slice(iter) = Resqnew_Vt.slice(iter) - norm_Vt*arma::ones(1,11);
            Resqnew_Vt.slice(iter).row(0) =Resqnew_Vt.slice(iter).row(0) - ((1-tauList.tail(1))/bL_vt - tauList.head(1)/b1_vt )*arma::ones(1,11);
        

  
    // tail parameters: b1_ and bL_
  
    Vect_b1_ut = Ut_2T - Udraw_lag * Resqnew_Ut.slice(iter).col(0);
    Vect_bL_ut = Ut_2T - Udraw_lag * Resqnew_Ut.slice(iter).col(10);
    ids_utb1_num = find(Vect_b1_ut <= 0);
    ids_utbL_num = find(Vect_bL_ut >= 0);
    zeros_15u_b1.elem(ids_utb1_num).fill(1.0); 
    zeros_15u_bL.elem(ids_utbL_num).fill(1.0);
    denom_15_b1 = Vect_b1_ut % ( Vect_b1_ut <= 0 );
    denom_15_bL = Vect_bL_ut % ( Vect_bL_ut >= 0 );
    b1_ut = - arma::sum(zeros_15u_b1) / arma::sum( denom_15_b1 ) ;
    bL_ut =   arma::sum(zeros_15u_bL) / arma::sum( denom_15_bL ) ; 


    Vect_b1_u1 = U1 - regressor_U1 * Resqnew_U1.slice(iter).col(0);
    Vect_bL_u1 = U1 - regressor_U1 * Resqnew_U1.slice(iter).col(10);
    ids_u1b1_num = find(Vect_b1_u1 <= 0);
    ids_u1bL_num = find(Vect_bL_u1 >= 0);
    zeros_1u_b1.elem(ids_u1b1_num).fill(1.0); 
    zeros_1u_bL.elem(ids_u1bL_num).fill(1.0);
    denom_1_b1 = Vect_b1_u1 % ( Vect_b1_u1 <= 0 );
    denom_1_bL = Vect_bL_u1 % ( Vect_bL_u1 >= 0 );
    b1_u1 = - arma::sum(zeros_1u_b1) / arma::sum( denom_1_b1 ) ;
    bL_u1 =   arma::sum(zeros_1u_bL) / arma::sum( denom_1_bL ) ; 

    Vect_b1_vt = Vt_2T - Vdraw_lag * Resqnew_Vt.slice(iter).col(0);
    Vect_bL_vt = Vt_2T - Vdraw_lag * Resqnew_Vt.slice(iter).col(10);
    ids_vtb1_num = find(Vect_b1_vt <= 0);
    ids_vtbL_num = find(Vect_bL_vt >= 0);
    zeros_15v_b1.elem(ids_vtb1_num).fill(1.0); 
    zeros_15v_bL.elem(ids_vtbL_num).fill(1.0);
    denom_15_b1 = Vect_b1_vt % ( Vect_b1_vt <= 0 );
    denom_15_bL = Vect_bL_vt % ( Vect_bL_vt >= 0 );
    b1_vt = - arma::sum(zeros_15v_b1) / arma::sum( denom_15_b1 ) ;
    bL_vt =   arma::sum(zeros_15v_bL) / arma::sum( denom_15_bL ) ; 

    Vect_b1_v1 = V1 - regressor_V1 * Resqnew_V1.slice(iter).col(0);
    Vect_bL_v1 = V1 - regressor_V1 * Resqnew_V1.slice(iter).col(10);
    ids_v1b1_num = find(Vect_b1_v1 <= 0);
    ids_v1bL_num = find(Vect_bL_v1 >= 0);
    zeros_1v_b1.elem(ids_v1b1_num).fill(1.0); 
    zeros_1v_bL.elem(ids_v1bL_num).fill(1.0);
    denom_1_b1 = Vect_b1_v1 % ( Vect_b1_v1 <= 0 );
    denom_1_bL = Vect_bL_v1 % ( Vect_bL_v1 >= 0 );
    b1_v1 = - arma::sum(zeros_1v_b1) / arma::sum( denom_1_b1 ) ;
    bL_v1 =   arma::sum(zeros_1v_bL) / arma::sum( denom_1_bL ) ; 

    // Update parameters to start the next iteration:

   Mat_init_Ut = Resqnew_Ut.slice(iter);
   Mat_init_U1 = Resqnew_U1.slice(iter);
   Mat_init_Vt = Resqnew_Vt.slice(iter);
   Mat_init_V1 = Resqnew_V1.slice(iter);

   // store tail parameters:

   mat_b(iter,0) = b1_ut; 
   mat_b(iter,1) = bL_ut;
   mat_b(iter,2) = b1_u1;
   mat_b(iter,3) = bL_u1; 
   mat_b(iter,4) = b1_vt;
   mat_b(iter,5) = bL_vt;
   mat_b(iter,6) = b1_v1;
   mat_b(iter,7) = bL_v1;

   // likelihood
   logg = arma::log(posterior_cpp2(last_Ut,utoty_nor,uc_resid,
                                   Herm_age,regressor_V1,regressor_U1,Mat_init_Vt,
                                   Mat_init_V1,Mat_init_Ut,Mat_init_U1,
                                   b1_ut,bL_ut,b1_u1,bL_u1,b1_vt,bL_vt,b1_v1,bL_v1,a,b,c,d,e,tauList,meanY,stdY));
  likelihood(iter) = arma::mean(logg);
 //if ( Rcpp::NumericVector::is_na(likelihood(iter))){ Rcpp::stop("Negative likelihood!"); }
   
   // reshape Obj_chain, Nu_chain, acc, acceptrate

   Obj_chain.col(0) = Obj_chain.col(draws-1);
   Obj_chain.cols(1,draws-1).zeros();

   for (int t = 0 ; t < T ; ++t){
         Nu_chain.slice(t).col(0) = Nu_chain.slice(t).col(draws-1);
         Nu_chain.slice(t).cols(1,draws-1).zeros();
         acc.slice(t).zeros();
   }

    Rcpp::Rcout << likelihood(iter) << "\n" ; 
  }
  Rcpp::List L = Rcpp::List::create(last_Ut,likelihood,logg,Resqnew_Ut, Resqnew_U1, 
                                    Resqnew_Vt, Resqnew_V1,mat_b,last_Ut,Obj_chain,Nu_chain,
                                    mean_accept);
 //arma::vec try1 = {b1_ut, bL_ut, b1_u1, bL_u1, b1_vt, bL_vt, b1_v1, bL_v1};
return L ;



}


//[[Rcpp::export]]

Rcpp::List innerLoop(int draws, int N, int T,
                     arma::vec tauList,
                     double var_prop1, double var_prop2, double var_prop3,
                     double var_prop4, double var_prop5, double var_prop6,
                     arma::mat Nu_chain1,arma::mat Nu_chain2,arma::mat Nu_chain3,
                     arma::mat Nu_chain4,arma::mat Nu_chain5,arma::mat Nu_chain6,
                     arma::mat Obj_chain,arma::mat acc1,arma::mat acc2,arma::mat acc3,
                     arma::mat acc4,arma::mat acc5,arma::mat acc6, 
                      arma::vec utoty_nor,arma::vec uc_resid,arma::mat Herm_age,
                    arma::mat regressor_V1,arma::mat regressor_U1,arma::mat Mat_init_Vt, 
                    arma::mat Mat_init_V1,arma::mat Mat_init_Ut, arma::mat Mat_init_U1,
                    double b1_ut,double bL_ut,double b1_u1,double bL_u1,double b1_vt, 
                  double bL_vt,double b1_v1,double bL_v1,double a,double b, double c,double d, double e,
                   double meanY, double stdY){

 // Initialize matrices: Nu_chain, Obj_chain, acc, acceptrate
  arma::mat Matdraws(N,T, arma::fill::zeros);
  arma::vec onesVec= arma::ones(N);
  arma::vec zerosVec = arma::zeros(N);

  // cubes,matrices and vectors used in the E-step
  arma::cube Nu_chain(N,draws,T,arma::fill::zeros);
  arma::cube acc(N,draws,T,arma::fill::zeros);

  Nu_chain.slice(0) = Nu_chain1;
  Nu_chain.slice(1) = Nu_chain2;
  Nu_chain.slice(2) = Nu_chain3;
  Nu_chain.slice(3) = Nu_chain4;
  Nu_chain.slice(4) = Nu_chain5;
  Nu_chain.slice(5) = Nu_chain6;
  acc.slice(0) = acc1;
  acc.slice(1) = acc2;
  acc.slice(2) = acc3;
  acc.slice(3) = acc4;
  acc.slice(4) = acc5;
  acc.slice(5) = acc6;

  //arma::mat Obj_chain(N,draws,arma::fill::zeros); 

  //arma::mat randnn = arma::join_rows(rand2,rand3,rand4,rand5); 
  //arma::mat randnn2 = arma::join_rows(randnn,rand6);
  arma::mat acceptrate(draws,T,arma::fill::zeros);
  
  arma::mat last_Ut(N,T,arma::fill::zeros);

  // matrices and vectors used in the inner loop: Metropolis-Hastings Algorithm
  arma::vec sd_mh = {sqrt(var_prop1),sqrt(var_prop2),sqrt(var_prop3),sqrt(var_prop4),
  sqrt(var_prop5), sqrt(var_prop6)};
  arma::vec newObj(N,arma::fill::zeros);
  arma::vec r(N,arma::fill::zeros);
  arma::vec prob(N,arma::fill::zeros);

    for (int j  = 1 ; j < draws; ++j){
      Matdraws.zeros();
      for (int t = 0 ; t < T; ++t){
             Matdraws.col(t) = Nu_chain.slice(t).col(j-1);
      }
          
          // U_1

          Matdraws.col(0) = Nu_chain.slice(0).col(j-1) + sd_mh(0) * arma::randn(N);
          newObj = posterior_cpp2(Matdraws,utoty_nor,uc_resid,Herm_age,regressor_V1,regressor_U1,Mat_init_Vt,Mat_init_V1,Mat_init_Ut,Mat_init_U1,b1_ut,bL_ut,b1_u1,bL_u1,b1_vt,bL_vt,b1_v1,bL_v1,a,b,c,d,e,tauList,meanY,stdY);
          if (  newObj.has_nan() == TRUE  ){ Rcpp::stop("newObj NA in 1"); }
          
          r = arma::min(onesVec,newObj/Obj_chain.col(j-1));
          prob = arma::randu(N);
        // prob = probl.col(0); 
          acc.slice(0).col(j) = arma::conv_to<arma::vec>::from( prob <= r );
          Obj_chain.col(j) = acc.slice(0).col(j) % newObj + ( onesVec - acc.slice(0).col(j)  ) % Obj_chain.col(j-1);
          Nu_chain.slice(0).col(j) = acc.slice(0).col(j) % Matdraws.col(0) + (onesVec - acc.slice(0).col(j)) % Nu_chain.slice(0).col(j-1);
          Matdraws.col(0) = Nu_chain.slice(0).col(j);


          // U_2 - U_6

          for ( int t = 1 ; t < T ; ++t ){
               Matdraws.col(t) = Nu_chain.slice(t).col(j-1) + sd_mh(t) *arma::randn(N);
               newObj = posterior_cpp2(Matdraws,utoty_nor,uc_resid,Herm_age,regressor_V1,regressor_U1,Mat_init_Vt,Mat_init_V1,Mat_init_Ut,Mat_init_U1,b1_ut,bL_ut,b1_u1,bL_u1,b1_vt,bL_vt,b1_v1,bL_v1,a,b,c,d,e,tauList,meanY,stdY);
          if (  newObj.has_nan() == TRUE  ){ Rcpp::stop("newObj NA in t"); }
          r = arma::min(onesVec,newObj/Obj_chain.col(j));
          prob = arma::randu(N);
         // prob = probl.col(t);
          acc.slice(t).col(j) = arma::conv_to<arma::vec>::from( prob <= r );
          Obj_chain.col(j) = acc.slice(t).col(j) % newObj + ( onesVec - acc.slice(t).col(j) )% Obj_chain.col(j);
          Nu_chain.slice(t).col(j) = acc.slice(t).col(j) % Matdraws.col(t) + (onesVec - acc.slice(t).col(j)) % Nu_chain.slice(t).col(j-1);
          Matdraws.col(t) = Nu_chain.slice(t).col(j);
          } 
          
          for ( int t = 0 ; t < T ; ++t ){
              acceptrate( j ,t ) = arma::mean( acc.slice(t).col(j)  );

          }

    } 
// E-step is over 
      
          // Keep the last draw 
    arma::mat last_Ut1 = join_rows(Nu_chain.slice(0).col(draws-1),Nu_chain.slice(1).col(draws-1),Nu_chain.slice(2).col(draws-1),Nu_chain.slice(3).col(draws-1));

    last_Ut = join_rows(last_Ut1, Nu_chain.slice(4).col(draws-1),Nu_chain.slice(5).col(draws-1));
    Obj_chain.col(0) = Obj_chain.col(draws-1);
    Obj_chain.cols(1,draws-1).zeros();

   for (int t = 0 ; t < T ; ++t){
         Nu_chain.slice(t).col(0) = Nu_chain.slice(t).col(draws-1);
         Nu_chain.slice(t).cols(1,draws-1).zeros();
         acc.slice(t).zeros();
   }


  Rcpp::List L = Rcpp::List::create(last_Ut,Nu_chain, Obj_chain,acceptrate );
 //arma::vec try1 = {b1_ut, bL_ut, b1_u1, bL_u1, b1_vt, bL_vt, b1_v1, bL_v1};
return L ;



   }





