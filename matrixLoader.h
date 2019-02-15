#pragma once

#ifndef _matrixLoader_H_
#define _matrixLoader_H_


#include <cstdlib>
#include <iomanip>
#include <cmath>
#include <vector>
#include <fstream>
#include <ctime>
#include <Eigen/Dense>


using namespace std;
using namespace Eigen;

template<class M>
M matrixLoader(const std::string & path) {
  std::ifstream indata;
  // NOTE: ifstream.open doesn't take std::string, convert it to cstyle string! 
  indata.open(path.c_str());
  std::string line;
  std::vector<double> values;
  int row_num = 0;
  while (std::getline(indata, line)) {
    std::stringstream lineStream(line);
    std::string cell;
    while (std::getline(lineStream, cell, ' ')) {
      // values.push_back(std::stod(cell));
      // NOTE: std::stod need c++11 above compiler, here have to convert to string then atof
      values.push_back(atof(cell.c_str()));
    }
    ++row_num;
  }
  int col_num = values.size()/row_num;
  return Map<const Matrix<typename M::Scalar, M::RowsAtCompileTime, M::ColsAtCompileTime, RowMajor> >(values.data(), row_num, col_num);
};

#endif
