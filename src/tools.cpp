#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  // ... your code here

  if(estimations.size() != ground_truth.size()
  		|| estimations.size() == 0){
  	cout << "Invalid estimation or ground_truth data" << endl;
  	return rmse;
  }

  //accumulate squared residuals
  for(int i=0; i < estimations.size(); ++i){
        // ... your code here
        if (estimations[i].cols() == ground_truth[i].cols() && estimations[i].rows() == ground_truth[i].rows()) {
            VectorXd residual = estimations[i] - ground_truth[i];
            //std::cout << "residual " << residual << std::endl;
            residual = residual.array().square();
            //std::cout << "residual2 " << residual << std::endl;
            rmse += residual;
        }
  }

  //std::cout << rmse << std::endl;

  //calculate the mean
  // ... your code here
    rmse /= (estimations.size());

  //calculate the squared root
  // ... your code here
    rmse = rmse.array().sqrt();

  //return the result
  return rmse;
}
