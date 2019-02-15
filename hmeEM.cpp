#include "hmeEM.h"


// constructor 1:
// ME simulation;
hmeEM::hmeEM(int pp, int pj, int pmode){
  cout << "ME simulation constructor called...\n";
  // this toy example is taken from the reference Jordan (1993);
  // support up to 5 linear experts with pp = 1
  mode = pmode;
  epsilon = 1e-3;
  convergence = false;
  iteration = 0;
  max_iter = 1000;
  
  j = pj;
  p = pp;
  n = pj * 200;

  X.resize(n , pp + 1);
  y.resize(n, 1);
  label.resize(n);
 
  H.resize(n, j);
  MU.resize(n, j);
  PI.resize(n, j);
  ETA.resize(n, j);
  PROB.resize(n, j);

  //---------------------initilaizer--------------------------------//
  hmeEM::parameter_initilaizer(p, j);
  
  //---------------------clean up existing log----------------------//
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
  cout << "suffix is " << suffix << endl;
  string base("output_log_");
  ofstream existing_log;
  existing_log.open(base + suffix + ".txt", ofstream::out|ofstream::trunc);
  existing_log.close();

  // bool log_exist = std::filesystem::exists('output_log.txt');
  // if (log_exist){remove('output_log.txt');}

  MatrixXd simulation_parameters(6, 4);
  // b0, b1, xl, xu, lambda(mix) = {0.4, 0.6}, sigma
  simulation_parameters(0, 0) = 0.4; //b0
  simulation_parameters(0, 1) = 0.8; //b1
  simulation_parameters(0, 2) = -1.0; // xl
  simulation_parameters(0, 3) = 2.0 ; //xu

  simulation_parameters(1, 0) = 2.4; //b0
  simulation_parameters(1, 1) = -0.5; //b1
  simulation_parameters(1, 2) = 1.0; // xl
  simulation_parameters(1, 3) = 4.0 ; //xu

  simulation_parameters(2, 0) = -1.2; //b0
  simulation_parameters(2, 1) = 0.4; //b1
  simulation_parameters(2, 2) = 3.0; // xl
  simulation_parameters(2, 3) = 6.0 ; //xu

  simulation_parameters(3, 0) = 1.2; //b0
  simulation_parameters(3, 1) = -0.1; //b1
  simulation_parameters(3, 2) = 5.0; // xl
  simulation_parameters(3, 3) = 8.0 ; //xu

  simulation_parameters(4, 0) = -1.4; //b0
  simulation_parameters(4, 1) = 0.1; //b1
  simulation_parameters(4, 2) = 7.0; // xl
  simulation_parameters(4, 3) = 10.0 ; //xu

  //this one is dummy, only to ease boundary
  simulation_parameters(5, 0) = -1.4; //b0
  simulation_parameters(5, 1) = 0.1; //b0
  simulation_parameters(5, 2) = 9.0; // xl
  simulation_parameters(5, 3) = 12.0 ; //xu
  

  double xL = simulation_parameters(0, 2); //xl of first
  double xU;



  if (j == 5) {
    xU = simulation_parameters( 4, 3) ; // xu of jth experts
  } else {
    xU = simulation_parameters( j - 1 , 3) - 1 ;
  }

  for (int i =0; i < n; i++){
    X(i, 0) = 1;
    X(i, 1) = (xU - xL)/double(n) * i + xL;
  }
  
  double sigma = 0.3;
  double mix = 0.4;
  random::mt19937 rng;
  random::normal_distribution<> norm(0, sigma);
  random::gamma_distribution<> gamma(1, 1);
  uniform_real<> r(0, 1);
  double z;
  variate_generator<random::mt19937&, random::normal_distribution<> > die(rng, norm);
  variate_generator<random::mt19937&, uniform_real<> > die_ratio(rng, r);
  variate_generator<random::mt19937&, random::gamma_distribution<> > die_gamma(rng, gamma);

  for (int i = 0; i < n; i++){
    rng.seed(i);
    for (int j = 0; j < pj; j ++){
      //-----------------single generating model -------------------// 
      
      if (j  == 0){
        if( X(i, 1) <= simulation_parameters(j + 1, 2)){
          z = die();
          y(i, 0) = simulation_parameters(j, 0) + simulation_parameters(j, 1) * X(i, 1) + z;
          label(i) = j;
        } else if (  X(i, 1)<= simulation_parameters(j , 3)) {
          z = die();
          mix = die_ratio();
          //--------------- generate from jth model -----------------// 
          if (mix < 0.4){
            y(i, 0) = simulation_parameters(j, 0) + simulation_parameters(j, 1) * X(i, 1) + z;
            label(i) = j ;
           //--------------- generate from j + 1 th model -----------// 
          } else{ 
            y(i, 0) = simulation_parameters(j + 1, 0) + simulation_parameters(j + 1, 1) * X(i, 1) + z;
            label(i) = j + 1;
        }
          
          }
         
          } else if (  simulation_parameters(j-1, 3) <= X(i, 1)  & X(i, 1) < simulation_parameters(j + 1, 2)){
        
        z = die();
        y(i, 0) = simulation_parameters(j, 0) + simulation_parameters(j, 1) * X(i, 1) + z;
        label(i) = j;
      }
      //---------------mixture generating model -----------------// 
      else if (simulation_parameters(j + 1, 2) <= X(i, 1) & X(i, 1) < simulation_parameters(j , 3) ){
       
        z = die();
        mix = die_ratio();
        //--------------- generate from jth model -----------------// 
        if (mix < 0.4){
          y(i, 0) = simulation_parameters(j, 0) + simulation_parameters(j, 1) * X(i, 1) + z;
          label(i) = j ;
        //--------------- generate from j + 1 th model -----------// 
            } else{ 
          y(i, 0) = simulation_parameters(j + 1, 0) + simulation_parameters(j + 1, 1) * X(i, 1) + z;
          label(i) = j + 1;
        }
     }
    }
  }
  //--------------------output simulation data-----------------//
  hmeEM::output_simulation();
 

  //----------------------execute EM--------------------------//

  cout << "Timing starts";
  clock_t t1,t2;
  t1=clock();
  fit_hmeEM(pmode);
  t2=clock();
  float diff = ((float)t2-(float)t1)/CLOCKS_PER_SEC;
  cout << "Running time: " << diff << endl;
    }


// constructor 2:
// input X, y
hmeEM::hmeEM(MatrixXd & pX, MatrixXd & py, int pj, int pmode){
  cout << "Data read in " << endl;
  mode = pmode;
  epsilon = 1e-3;
  convergence = false;
  iteration = 0;
  max_iter = 1000;

  j = pj;
  n = pX.rows();
  p = pX.cols();
  X.resize(n, p + 1);
  y = py;
 
  VectorXd ones = VectorXd::Ones(n);
  X << ones, pX;


  H.resize(n, j);
  MU.resize(n, j);
  PI.resize(n, j);
  ETA.resize(n, j);
  PROB.resize(n, j);

  //---------------------initilaizer--------------------------------//
  hmeEM::parameter_initilaizer(p, j);
  
  //---------------------clean up existing log----------------------//
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
  cout << "suffix is " << suffix << endl;
  string base("output_log_");
  ofstream existing_log;
  existing_log.open(base + suffix + ".txt", ofstream::out|ofstream::trunc);
  existing_log.close();
  //---------------
  cout << "Timing starts";
  clock_t t1,t2;
  t1=clock();
  fit_hmeEM(pmode);
  t2=clock();
  float diff = ((float)t2-(float)t1)/CLOCKS_PER_SEC;
  cout << "Running time: " << diff << endl;
  /////------execute EM-----------------------------//
  
}
