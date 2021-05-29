// [[Rcpp::depends(RcppEigen,RcppArmadillo)]]
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
    Eigen::VectorXd bool_v = (Y.array() < xbeta.array()).cast<double>();
    Eigen::VectorXd vv = tau -  bool_v.array();
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
  int status = optim_lbfgs(nll, beta, fopt);
  if(status < 0)
    Rcpp::stop("fail to converge");
  
  return Rcpp::wrap(beta);
}


//SEXP rq_cpp(Rcpp::NumericMatrix X, Rcpp::NumericVector y, Rcpp::NumericVector tauList  ){
//     arma::mat coef_rq_mat;
//  for (int i = 0; i < 11 ; ++i){
      
//      arma::vec output = Rcpp::as<arma::colvec>( q_reg(X,y,tauList(i)) );
//      coef_rq_mat.col(i) = output;
//    }
//  return Rcpp::wrap(coef_rq_mat);
//  }
  
  
//Rcpp::NumericMatrix rq_cpp(Rcpp::NumericMatrix X, Rcpp::NumericVector y,   Rcpp::NumericVector tauList ){
//  Rcpp::NumericMatrix coef_rq( X.cols()  ,11  );
//  for (int i = 0 ; i < 11 ; ++i){
//     Rcpp::NumericVector col_qreg=q_reg(X,y,tauList(i));
//     Rcpp::NumericMatrix::Column col = coef_rq( _ , i )   ;
//     coef_rq( _ , i ) = 3;
//     col = col_qreg;
//  }
//  return Rcpp::wrap(coef_rq);
//}



  
  