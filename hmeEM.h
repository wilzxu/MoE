#pragma once

#ifndef _hmeEM_H_ 
#define _hmeEM_H_

#include <cstdlib>
#include <iomanip>
#include <cmath>
#include <vector>
#include <fstream>
#include <ctime>
#include <limits>
#include <Eigen/Dense>
#include <boost/random/uniform_real.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions/normal.hpp>
#include <Eigen/QR>    
using namespace std;
using namespace Eigen;
using namespace boost;


// enabled when compile with c++11
// struct MyException : public std::exception {
//   const char * what () const throw () {
//     return "C++ Exception";
//   }
// };

class hmeEM {
public:
  MatrixXd X; 
  MatrixXd y;
  VectorXi label;
  int j;
  int p;
  int k;
  int mode;

  // constructor 1: 
  // ME simulation mode;
  hmeEM(int pp, int pj, int pmode); 

  // constructor 2 (to be built):
  // HME simulation mode;
  hmeEM(int pp, int pj, int pk, int pmode);


  // constructor 3:
  // real data mode:
  hmeEM(MatrixXd& pX, MatrixXd& py, int pj, int pmode);
  //------------------------getters----------------------------//
  MatrixXd get_BETA(){return BETA;} 
  MatrixXd get_THETA(){return THETA;} 
  MatrixXd get_SIGMA2(){return SIGMA2;} 
  MatrixXd get_H(){return H;} 

private:
  int n; // sample size
  double loglik;
  MatrixXd BETA;
  MatrixXd THETA;
  MatrixXd GAMMA;
  MatrixXd SIGMA2;
  MatrixXd all;

  MatrixXd H; //posterierS=s;
  MatrixXd H_top; //posterierS=s;

  
  MatrixXd MU; //X * BETA;
  MatrixXd PI; //inv.logit(X * THETA);
  MatrixXd XI; //inv.logit(X * GAMMA);
  MatrixXd ETA; //X * THETA;
  MatrixXd ETA_top; //X * GAMMA;
  MatrixXd PROB; // p(y|X, BETA);
  double epsilon; //tolerance;
  int iteration; 
  bool convergence; 
  int max_iter; 

  void parameter_initilaizer(int pp, int pj);
  int fit_hmeEM(int pmode);
  void Estep();
  void Mstep(int pmode);
  MatrixXd qNR(MatrixXd H, MatrixXd PI, MatrixXd PROB, MatrixXd THETA, double loglik);
  void check_convergence();
  void output_parameters();
  void output_simulation();
  void output_log();
  double evaluateLoglik(MatrixXd &H, MatrixXd &THETA );
  // void svdDecomp(MatrixXd &X){
  //   JacobiSVD<MatrixXd> svd(X, ComputeThinU|ComputeThinV);
  //   this->V = svd.matrixV();
  //   this->D = svd.singularValues()    cout << "Please choose the optimizer: "<< endl;
;
  // }
  // MatrixXd V ; //SVD of X
  // VectorXd D ; //SVD of X
  MatrixXd clamped_softmax(MatrixXd &X, double lower);
  MatrixXd clamped_linear(MatrixXd &X, double upper, double lower);
  void hessianApprox(MatrixXd & C);
  double lineSearch(MatrixXd & GRADIENT , MatrixXd &THETA, Map<Eigen::MatrixXd> & S,  double & target_previous);
  double lineSearch(MatrixXd & GRADIENT , MatrixXd &THETA, double & target_previous);
  vector<double> range(double lower, double upper, int number);
  MatrixXd gradientDescent(MatrixXd H, MatrixXd PI, MatrixXd PROB, MatrixXd THETA, double loglik);
  //------------prepare for pseudo inverse of diagonal matrix (in VectorXd format)--------------//
  VectorXd v_pinv(VectorXd v){
    for (int i = 0; i< v.size(); i++){
      if (v(i) > 1e-10) {
        v(i) = v(i);
      } else {
        v(i) = 0;
      }
    }
    return v;
}



  
};



void hmeEM::parameter_initilaizer(int pp, int pj){
  double sigma = 0.3;
  random::mt19937 rng;
  rng.seed(61515);
  random::normal_distribution<> norm(0, sigma);
  random::gamma_distribution<> gamma(1, 1);
  uniform_real<> r(0, 1);
  double z;
  variate_generator<random::mt19937&, random::normal_distribution<> > die(rng, norm);
  variate_generator<random::mt19937&, uniform_real<> > die_ratio(rng, r);
  variate_generator<random::mt19937&, random::gamma_distribution<> > die_gamma(rng, gamma);

  //-----------------------BETA---------------------------//
  BETA.resize(pp + 1, pj);

  for (int i = 0; i < pp + 1; i++){
    for (int j = 0; j < pj ; j++){
      BETA(i, j) = die();
    }
  }
  //----------------------THETA--------------------------//
  THETA.resize(pp + 1, pj);

  for (int i = 0; i < pp + 1; i++){
    for (int j = 0; j < pj ; j++){
      THETA(i, j) = die();
    }
  }

    
  //----------------------SIGMA2--------------------------//
  SIGMA2.resize(1, pj);
  for (int i = 0; i < pj ; i ++){
    SIGMA2(0, i) = die_gamma();
  }
  
 
  all.resize(BETA.rows() + THETA.rows() + SIGMA2.rows(), pj);
  all << BETA,
    THETA,
    SIGMA2;

  loglik = (PI.array() * PROB.array()).rowwise().sum().log().sum();
  }

void hmeEM::output_simulation(){
    //---------------output X (intercept included)----------------//
    ofstream myfile;
    myfile.open("X.txt");
    myfile << X <<endl;
    myfile.close();

    //----------output y and the label of data generating model-----//
    myfile.open("Y.txt");
    for (int i = 0; i < n; i++){
      myfile << y(i, 0) << "\t" << label(i) << endl;
    }
    myfile.close();
  } 

void hmeEM::output_parameters(){
  ofstream myfile;
  //-----------------------output model parameters-----------------//
  string suffix;
  switch(mode){
  case 1:
    {
      suffix = "gd";
      break;
    }
  case 2:
    {
      suffix = "vN";
      break;
    } 
  case 3:
    {
      suffix = "qN";
      break;
    }    
  }
  string base("output_parameters_");
  myfile.open(base + suffix + ".txt");
  myfile << all <<endl;
  myfile.close();
  //-----------------------output final H--------------------------//
  base = "output_H_";
  myfile.open(base + suffix + ".txt");
  myfile << H << endl;
  myfile.close();
  //-----------------------output log-------------------------------//
  hmeEM::output_log();
}

void hmeEM::output_log(){
  Matrix<double,Dynamic,Dynamic,RowMajor> tmp(all);
  Map<RowVectorXd> all_vector(tmp.data(), tmp.size());
  
  ofstream myfile;
  //-----------------------output log------------------------------//
  string suffix;
  switch(mode){
  case 1:
    {
      suffix = "gd";
      break;
    }
  case 2:
    {
      suffix = "vN";
      break;
    } 
    case 3:
      {
        suffix = "qN";
        break;
      }
  }
  string base("output_log_");
  myfile.open(base + suffix + ".txt", fstream::app);
  myfile << all_vector <<"\t" << loglik << endl;
  myfile.close();
}



int hmeEM::fit_hmeEM(int pmode){
    
    while (!convergence & (iteration < max_iter)){
    iteration++;
    cout << "******************iteration******************* " << endl;
    cout << "iteration : " << iteration << endl;  
    hmeEM::Estep();
    cout << "Estep" << endl;
    hmeEM::Mstep(pmode);
    cout << "Mstep" << endl;
    hmeEM::check_convergence();
    cout << "convergence : " << convergence << endl;
    cout << "//-------------v8--------------//"<< endl;
  }
 
  return 0;
};

void hmeEM::Estep(){
  cout << "loglik = " << loglik <<endl;
 
  // try{
    MU = X * BETA;
    //------------------------PROB----------------------------//
    for (int i = 0; i < n; i++){
      for (int e = 0; e < j; e++){
        math::normal_distribution<> prob(MU(i, e), pow(SIGMA2(0, e), 0.5));
        PROB(i, e) = pdf(prob, y(i));
      }
    }
    //-----------------------PI-------------------------------//
    ETA = X * THETA;
    PI = ETA.array().exp();
    hmeEM::clamped_softmax(PI, -15);
    //--------------------------H-------------------------------//
    H = PI.array() * PROB.array(); 
    for (int i = 0; i < n; i++){
      double m = H.row(i).maxCoeff();
      double log_m = log(m);
      H.row(i) = exp(H.row(i).array().log()-log_m) / exp(log(H.row(i).sum()) - log_m);
      // for (int e =0; e < j; e++){
      //   if (H(i, e) < 1e-10){H(i, e) = 0; }
      // }
      // H.row(i) = H.row(i)/H.row(i).sum();
    
    }
    //------------------------loglik---------------------------//
    loglik = (PI.array() * PROB.array()).rowwise().sum().log().sum();

  
    // auto a = H.array().isNaN();
    // bool f = (a.array() != 0).any();
    // if (f){
  //     throw MyException();
  // }
// }
    
// catch(MyException& e)
//   {
//     cout << "NaN exits, ending model " << endl;
//     abort();
     
// }
    hmeEM::output_parameters();
    hmeEM::output_log();

}

void hmeEM::Mstep(int pmode){

  MatrixXd X_tilta;
  MatrixXd y_tilta;
  VectorXd res;
  VectorXd v;
  for (int e = 0 ; e < j ; e++){
    //-------------------update SIGMA2--------------------------//
    res = y.col(0).array() - MU.col(e).array();
    double tmp = H.col(e).sum();
    double tmp2 = res.transpose() * H.col(e).asDiagonal() * res ;
    SIGMA2(0, e) = (double) tmp2/tmp;
    SIGMA2(0, e) = (SIGMA2(0, e) < 1e-10) ? 1e-10: SIGMA2(0, e);
    //-------------------update BETA---------------------------//
    v = H.col(e).array().pow(0.5);
    X_tilta = v.asDiagonal() * X;
    y_tilta = v.asDiagonal() * y;
    JacobiSVD<MatrixXd> svd(X_tilta,ComputeThinU|ComputeThinV);
    BETA.col(e) = svd.solve(y_tilta);
  }
    //-------------------update THETA--------------------------//
    switch(pmode)
      {
        
      case 1 : 
        {
    // --------------------gradient descent-------------------------//
          THETA = gradientDescent(H, PI, PROB, THETA, loglik);
          break;
        }
      case 2 :
        {
          for (int e = 0; e < j; e++){
    // --------------variant NR solved by SVD-----------------------//
            v = PI.col(e).array() * (X.col(0).array() - PI.col(e).array());
            v = v.array().pow(0.5);
            X_tilta = v.asDiagonal() * X;
            y_tilta = v.asDiagonal() * y;
            v = v_pinv(v);
            JacobiSVD<MatrixXd> svd_2(X_tilta,ComputeThinU|ComputeThinV);
            THETA.col(e) = svd_2.solve(X_tilta * THETA.col(e) + v.asDiagonal()* (H.col(e) - PI.col(e)));
          } 
     
          break;
        }

        // --------------------quasi Newton-------------------------//
      case 3:
        {
        THETA = qNR(H, PI, PROB, THETA, loglik);
        break;
        }
      }    


}

MatrixXd hmeEM::qNR(MatrixXd H, MatrixXd PI, MatrixXd PROB, MatrixXd THETA, double loglik){
  //Initialize C;  
  MatrixXd C = MatrixXd::Zero( (p + 1) * j, (p + 1) * j);
  hessianApprox(C);
  
  // MatrixXd C = MatrixXd::Identity( (p + 1) * j, (p + 1) * j);
  MatrixXd PI_new(n , j);
  MatrixXd THETA_new(p + 1, j);
 
  VectorXd s((p + 1) * j);
  VectorXd y((p + 1) * j);
  MatrixXd ETA_new(p + 1, j);

  
  double target_previous = -evaluateLoglik(H, THETA);
  
  for (int iter = 0; iter <5; iter++){
    //---------------------compute THETA_new-----------------------//
    // cout << "inner interation : " << iter << endl;
   
      MatrixXd GRADIENT_previous = - X.transpose() *  (H.array() - PI.array()).matrix();
      Map<VectorXd> gradient_previous(GRADIENT_previous.data(), GRADIENT_previous.size());
      clamped_linear(C, 1e6, -1e6);
      s = - C * gradient_previous;
      Map<MatrixXd> S(s.data(), THETA.rows(), THETA.cols());
      double alpha = lineSearch(GRADIENT_previous, THETA, S, target_previous);

      cout << iter << " inner alpha chosen :" << alpha << endl;
      if (alpha == 0) break;
      s = s.array() * alpha;
       //-------------------- makesure no non-----------------------//
     
      THETA_new = THETA.array() + S.array() * alpha;
      clamped_linear(THETA_new, 10, -10);//on log scale
      PI_new = (X * THETA_new).array().exp();
      clamped_softmax(PI_new, -20);

      MatrixXd GRADIENT_new = -X.transpose() *  (H.array() - PI_new.array()).matrix();
      Map<VectorXd> gradient_new(GRADIENT_new.data(), GRADIENT_new.size());
      //---------------------compute y, s-------------------------//
      y = gradient_new - gradient_previous;
      // if (y.array().sum() == 0) break;
      //---------------------update C------------------------// 
      MatrixXd c_left = MatrixXd::Identity((p + 1) * j, (p + 1) * j).array() - (s * y.transpose()).array()/(y.transpose() * s);
      MatrixXd c_right = MatrixXd::Identity((p + 1) * j, (p + 1) * j).array() - (y * s.transpose()).array()/(y.transpose() * s);
      C = (c_left * C * c_right).array() + (s * s.transpose()).array()/(y.transpose() * s);

      target_previous = -evaluateLoglik(H, THETA_new);
      THETA = THETA_new;
      PI = PI_new;
  }
  
  
  return THETA_new;
}
  

//line search for qNR
double hmeEM::lineSearch(MatrixXd & GRADIENT, MatrixXd &THETA, Map<Eigen::MatrixXd> & S,  double & target_previous){
  double alpha = 1; 
  double target_new = DBL_MAX;
  for (int i = 0; i < 100; i++) {
    MatrixXd THETA_new = THETA.array() + alpha * S.array();
    clamped_linear(THETA_new, 10, -10);
    target_new = -evaluateLoglik(H, THETA_new);
    double delta =   1/2 * alpha * GRADIENT.squaredNorm();
    // double delta = 0;
    if (target_new < target_previous - delta){
      // target_previous = target_new;
      return alpha; 
    }
    alpha *= 0.5;
  }
  cout << "target previous" << target_previous << endl;
  cout << "target new "<< target_new << endl;
    
  return 0;
      
}


vector<double> hmeEM::range(double lower, double upper, int number){
  vector<double> alpha;
  for (int i =0; i < number; i++)
    alpha.push_back(lower + i/(double)number * (upper - lower));
  return alpha;
}


MatrixXd hmeEM::gradientDescent(MatrixXd H, MatrixXd PI, MatrixXd PROB, MatrixXd THETA, double loglik){

  MatrixXd THETA_new(p + 1, j);
  MatrixXd PI_new(n, j);
  double decay = 1;
  double target_previous = -evaluateLoglik(H, THETA);
  for (int iter = 0; iter <5; iter++){
    //---------------------compute THETA_new-----------------------//
    // cout << "inner interation : " << iter << endl;
   
      MatrixXd GRADIENT_previous =  -X.transpose() *  (H.array() - PI.array()).matrix();

      double eta = decay * lineSearch(GRADIENT_previous, THETA, target_previous);

      cout << iter << "inner eta chosen : " << eta << endl;
      //-------------------- makesure no non-----------------------//
     
      THETA_new = THETA.array() - GRADIENT_previous.array() * eta;
      PI_new = (X * THETA_new).array().exp();
      clamped_softmax(PI_new, -15);
      PI = PI_new;
      THETA = THETA_new;
      target_previous = -evaluateLoglik(H, THETA_new);
      cout << "eta = " << eta << endl;
      if (eta == 0) {
        return THETA;
      }
      // if (iter % 100 ==0) decay *= 0.5;
  }
  return THETA;
}

//line search for gradient descent
double hmeEM::lineSearch(MatrixXd & GRADIENT, MatrixXd & THETA, double& target_previous){
 
  double eta = 1e-2;
  double target_new = DBL_MAX;
  int i = 0;
  while (i < 100){
    MatrixXd THETA_new = THETA.array() - eta * GRADIENT.array();
    clamped_linear(THETA_new, 10, -10);
    target_new = -evaluateLoglik(H, THETA_new);
    double delta = 1/2 * fabs(eta) * GRADIENT.squaredNorm();
    if (target_new < target_previous - delta){
      target_previous = target_new;
      return eta; 
      break;
    }
    i++;
    eta *= 0.5;
  }
  cout << i << "target previous" << target_previous << endl;
  cout << "target new "<< target_new << endl;
    
  return 0;
}

void hmeEM::check_convergence(){
  MatrixXd all_previous = all;
  all << BETA,
    THETA,
    SIGMA2;
      
  cout << "norm : " <<(all_previous - all).norm() << endl;
  if ((all_previous-all).norm() < epsilon  ){
      convergence = true;
    };
}


double hmeEM::evaluateLoglik(MatrixXd &H, MatrixXd &THETA ){
  MatrixXd PI_tmp = exp((X * THETA).array());
  clamped_softmax(PI_tmp, -10);
  double loglik_new = (H.array() * PI_tmp.array().log()).sum();
  return loglik_new;
};

//time limiting 38.18%
MatrixXd hmeEM::clamped_softmax(MatrixXd &X, double lower){
  for (int i = 0; i < X.rows(); i++){
    double m = X.row(i).maxCoeff();
    double log_m = log(m);
    VectorXd tmp = X.row(i).array().log()-log_m;
    for (int j = 0 ; j< X.cols(); j++){
      tmp(j) = (tmp(j) < lower) ? lower: tmp(j);
    }
    X.row(i) = exp(tmp.array())/exp(tmp.array()).sum();
  }
  return X;
}

MatrixXd hmeEM::clamped_linear(MatrixXd &X, double upper, double lower){
  for (int i = 0; i < X.rows(); i++){
    VectorXd tmp = X.row(i).array();
    for (int j = 0 ; j< X.cols(); j++){
      tmp(j) = (tmp(j) < lower) ? lower: tmp(j);
      tmp(j) = (tmp(j) > upper) ? upper: tmp(j);
    }
    X.row(i) = tmp;
  }
  return X;
}

void hmeEM::hessianApprox(MatrixXd & C){
  
  MatrixXd tmp = PI.array() * (1 - PI.array());
  // tmp = tmp.block(0, 0, p + 1, j);
  // MatrixXd diag(p + 1 , p + 1);
  for (int e = 0; e < j; e ++){
  //   cout << "ok" << endl;
  //   diag = v_pinv(D.array().pow(2) * tmp.col(e).array());
  //   cout << "ok" << endl;
  // C.block(e * (p + 1), e * (p + 1), p + 1, p + 1) = V * diag * V.transpose();
    MatrixXd A = X.transpose() * (tmp.col(e).asDiagonal()) * X;
    MatrixXd pinvA =  A.completeOrthogonalDecomposition().pseudoInverse();
    C.block(e * (p + 1), e * (p + 1), p + 1, p + 1) = pinvA;
  } 
  cout << "c = " << C.determinant() << endl;
}


#endif
