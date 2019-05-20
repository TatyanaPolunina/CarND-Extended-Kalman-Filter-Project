#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) 
{

  VectorXd rmse = VectorXd::Zero(estimations[0].size());  
  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if (estimations.size() != ground_truth.size()
      || estimations.size() == 0) {
    std::cout << "Invalid estimation or ground_truth data" << std::endl;
    return rmse;
  }
  
  for (unsigned int i=0; i < estimations.size(); ++i) {
    VectorXd diff = estimations[i] - ground_truth[i];
    // coefficient-wise multiplication
    diff = diff.array()*diff.array();
    rmse += diff;
  }

  // calculate the mean
  rmse = rmse/estimations.size();
  // calculate the squared root
  rmse = rmse.array().sqrt();
  // return the result
  return rmse; 
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  MatrixXd Hj(3,4);
  if (x_state.size() != 4) {
    std::cout << "CalculateJacobian () - Error - Incorrect input" << std::endl;
    return Hj;
  }
  // recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  float sqr_sum = px * px + py * py;
  // check division by zero
  if (fabs(sqr_sum) < 0.0001) {
    std::cout << "CalculateJacobian () - Error - Division by Zero" << std::endl;
    return Hj;
  }

  float sqrt_sum = std::sqrt(sqr_sum);
  float sum32 = sqrt_sum * sqr_sum;
  
  // compute the Jacobian matrix
  Hj << (px/sqrt_sum), (py/sqrt_sum), 0, 0,
      -(py/sqr_sum), (px/sqr_sum), 0, 0,
      py*(vx*py - vy*px)/sum32, px*(px*vy - py*vx)/sum32, px/sqrt_sum, py/sqrt_sum;
  return Hj;
}

VectorXd Tools::CalculateH(const VectorXd& x) {
  VectorXd hx = VectorXd::Zero(3);
  
  if (x.size() != 4) {
    std::cout << "incorrect x_ size";
    return hx;
  }
  // recover state parameters
  float px = x(0);
  float py = x(1);
  float vx = x(2);
  float vy = x(3);
  float fi = std::atan(py/px);
  while (fi > M_PI)
    fi -= 2 * M_PI;
  while (fi < -M_PI)
    fi += 2 * M_PI; 
  float sqrt_sum = std::sqrt(px * px + py * py);
  hx << sqrt_sum, fi, (px*vx + py * vy) / sqrt_sum;
  return hx;
}

VectorXd Tools::PolarToCasterian(const VectorXd& polar) {
  VectorXd res = VectorXd::Zero(4);
  if (polar.size() < 3) {
    std::cout << "incorrect polars";
    return res;
  }
  res(0) = polar(0) * cos (polar(1));
  res(1) = polar(0) * sin (polar(1));
  return res;
}