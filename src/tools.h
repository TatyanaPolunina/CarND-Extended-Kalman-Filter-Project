#ifndef TOOLS_H_
#define TOOLS_H_

#include <vector>
#include "Eigen/Dense"

class Tools {
 public:

  /**
   * A helper method to calculate RMSE.
   */
  static Eigen::VectorXd CalculateRMSE(const std::vector<Eigen::VectorXd> &estimations, 
                                const std::vector<Eigen::VectorXd> &ground_truth);

  /**
   * A helper method to calculate Jacobians.
   */
  static Eigen::MatrixXd CalculateJacobian(const Eigen::VectorXd& x_state);

  /**
  * calculate converstion to polar
  */
  static Eigen::VectorXd CalculateH(const Eigen::VectorXd& x);
  
  static Eigen::VectorXd PolarToCasterian(const Eigen::VectorXd& polar);
  
  static void NormalizePolarVector(Eigen::VectorXd& polar);
};

#endif  // TOOLS_H_
