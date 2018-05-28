#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
	
	VectorXd rmse(4);
  rmse << 0, 0, 0, 0;
  
  //check the validity of the following inputs:
  //* the estimations vector size should not be zero
  //* the estimations vector size should equal ground truth vector size
  if(estimations.size() != ground_truth.size()
			|| estimations.size() == 0) {
	  cout << "Invalid estimations or ground_truth data" << endl;
	  return rmse;
  }
  
  //accumulate square residuals
  for(unsigned int i = 0; i < estimations.size(); ++i) {
	  
	  VectorXd residual = estimations[i] - ground_truth[i];
	  
	  //coefficient-wise multiplication
	  residual = residual.array() * residual.array();
	  rmse += residual;
  }
  
  //calculate the mean
  rmse = rmse / estimations.size();
  
  //calculate the squared root
  rmse = rmse.array().sqrt();
  
  //return the result
  return rmse;
}