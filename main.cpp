 // *******************************************************************************
 // * Optimizers for Mixture of Experts Model
 // * This c++ library was created to provide optimization tools for Mixture of
 // * Experts model, available methods: gradient descent, quasi-Newton Method, variant
 // * of Newton's method
 // *
 // * Version 1.2.0
 // * Dec 10, 2018
 // *
 // *
 // *  Dependencies:Boost, Eigen
 // *
 // * hmeEM: Defines the entry point for the hmeEM application
 // * To compile:
 // * g++ -std=c++11 -O -o main -I /path/to/dependencies main.cpp
 // *
 // * To use hmeEM, please use one of following commands:
 // *
 // * ./main X.txt y.txt (inputs your own X, y, see X_NO.txt, Y_NO.txt for example)
 // * ./main (simulate toy X, y)
 // *******************************************************************************/




#include "hmeEM.cpp"
#include "hmeEM.h"
// #include "hmeEM_qNRv5_gd.h"
// #include "hmeEM_copy.h"
// #include "hmeEMold.h"
#include "matrixLoader.h"
#include <time.h>

int main(int argc, char* argv[]){
    cout << "###################################################################" <<endl;
    cout << "Welcome to our optimizatation tool for Mixtures of Experts Model! " << endl;
  // reading in X, y, j;
    cout <<  "##################################################################" <<endl;


  if (argc > 1){
    // real data mode;
      cout << "Please enter the number of experts: "<< endl;
      int j;
      cin >> j;
      cout << "Please choose the optimizer: "<< endl;
      cout << "1. gradient descent" <<endl;
      cout << "2. variant Newton's method" <<endl;
      cout << "3. quasi-Newton's method" <<endl;
      int mode;
      cin >> mode;
      MatrixXd X = matrixLoader<MatrixXd>(argv[1]);
      MatrixXd y = matrixLoader<MatrixXd>(argv[2]);
      cout << "Initiating..." << endl;
      // hmeEM a = hmeEM(X, y, j, mode);
      // cout << "BETA = " << a.get_BETA() << endl;
      // cout << "SIGMA = " << a.get_SIGMA2().array().pow(0.5)<<endl;
    } else {

   
    // ME simulation mode;
    int p = 1; // only support p = 1 in toy example;
    bool intercept = true; // assuming always true;
    int j;
    cout << "Please enter the number of experts: " << endl;
    cout << "ME simulation mode supports up to j = 5" << endl;
    cout << "j = " ;
    cin >> j;
    cout << "Please choose the optimizer: "<< endl;
    cout << "1. gradient descent" <<endl;
    cout << "2. variant Newton's method" <<endl;
    cout << "3. quasi-Newton's method" <<endl;
    int mode;
    cin >> mode;
    cout << "Simulate data dimension of X: p = " << p << endl;
    cout << "Simulate data with intercept: " << intercept << endl;

    cout << "Initiating ME..." << endl;
    hmeEM a= hmeEM(p, j, mode);
    // hmeEM a= hmeEM(p, j);
    cout << "BETA = " << a.get_BETA() << endl;
    cout << "SIGMA = " << a.get_SIGMA2().array().pow(0.5)<<endl;
  }
  
  cout << "end" << endl;

  return 0;
  
}
