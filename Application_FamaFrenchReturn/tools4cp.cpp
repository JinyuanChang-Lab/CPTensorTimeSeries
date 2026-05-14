#include <RcppEigen.h>
#include <Rcpp.h>
#include <iostream>
#include <algorithm>
// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
Eigen::MatrixXd sigmak(Eigen::MatrixXd Y, Eigen::MatrixXd Y_mean, int k ,int n){
  double nn = n;
  Y = Y-Y_mean.replicate(1,n);
  //Rcout<<Y_mean.replicate(1,n);
  Eigen::MatrixXd Cov_lagk = Y.rightCols(n-k)*Y.leftCols(n-k).transpose()/nn;
  return Cov_lagk;
}
// [[Rcpp::export]]
Eigen::MatrixXd thresh_C(Eigen::MatrixXd mat, double delta){
  //double threshold = lambda * sqrt(log(p) / n);
  for (int i = 0; i < mat.rows(); i++) {
    for (int j = 0; j < mat.cols(); j++) {
      if ( abs(mat(i, j)) < delta) {
        mat(i, j) = 0;
      }
    }
  }
  return mat;
}


// [[Rcpp::export]]
SEXP MatMult(Eigen::MatrixXd A, Eigen::MatrixXd B){
  Eigen::MatrixXd C = A * B;
  return Rcpp::wrap(C);
}


// [[Rcpp::export]]
Eigen::VectorXd minor_P(Eigen::MatrixXd Wr, Eigen::MatrixXd Ws, int d1, int d2){
  
  Eigen::VectorXd P(d1*d1*d2*d2);
  int sss = 0;
  for(int l = 0; l < d2;l++){
    for(int k = 0; k < d2;k++){
      for(int j = 0; j < d1;j++){
        for(int i = 0; i < d1;i++){
          P(sss) = Wr(i,k) * Ws(j,l) + Ws(i,k) * Wr(j,l) - Wr(i,l) * Ws(j,k) - Ws(i,l) * Wr(j,k);
          sss++;
        }
      }
    }
  }
  
  return P;
  
}


// [[Rcpp::export]]
Eigen::VectorXd minor_P_vector(Eigen::VectorXd Wr, Eigen::VectorXd Ws, int d1, int d2){
  
  Eigen::VectorXd P(d1*d1*d2*d2);
  int sss = 0;
  for(int l = 0; l < d2;l++){
    for(int k = 0; k < d2;k++){
      for(int j = 0; j < d1;j++){
        for(int i = 0; i < d1;i++){
          P(sss) = Wr(i,k) * Ws(j,l) + Ws(i,k) * Wr(j,l) - Wr(i,l) * Ws(j,k) - Ws(i,l) * Wr(j,k);
          sss++;
        }
      }
    }
  }
  
  return P;
  
}

// [[Rcpp::export]]
Eigen::MatrixXd Vech2Mat(Eigen::VectorXd P, int d){
  
  Eigen::MatrixXd M  = Eigen::MatrixXd::Zero(d,d);
  Eigen::MatrixXd Md = Eigen::MatrixXd::Zero(d,d);
  int k = 0;
  for(int j = 0; j < d;j++){
    for(int i = j; i < d;i++){
      M(i,j) = P(k);
      k++;
    }
  }
  
  Md.diagonal() = M.diagonal();
  M = M + M.transpose() - Md;
  
  return M;
  
}

// [[Rcpp::export]]
Eigen::MatrixXd Vech2Mat_new(Eigen::VectorXd P, int d){
  
  Eigen::MatrixXd M  = Eigen::MatrixXd::Zero(d,d);
  
  int k = 0;
  for(int j = 0; j < d;j++){
    for(int i = j; i < d;i++){
      M(i,j) = P(k)/2;
      k++;
    }
  }
  
  M = M + M.transpose();
  
  return M;
  
}





