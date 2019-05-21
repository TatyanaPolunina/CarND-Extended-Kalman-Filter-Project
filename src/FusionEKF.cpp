#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  P_ = MatrixXd(4, 4);
  Hj_ = MatrixXd::Zero(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;


   H_laser_ << 1, 0, 0, 0,
            0, 1, 0, 0;
   
    P_ << 1, 0, 0, 0,
          0, 1, 0, 0,
          0, 0, 1000, 0,
          0, 0, 0, 1000;
   noise_ax = 9;
   noise_ay = 9;
}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    /**
     *  Create the covariance matrix.     
     */
     std::cout << " First measurement ";
     // the initial transition matrix F_
    MatrixXd F_(4, 4);
    F_ << 1, 0, 1, 0,
            0, 1, 0, 1,
            0, 0, 1, 0,
            0, 0, 0, 1;
	MatrixXd Q = MatrixXd::Zero(4,4);
    // first measurement
    cout << "EKF: " << endl;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // Convert radar from polar to cartesian coordinates 
      //         and initialize state.
      VectorXd x = Tools::PolarToCasterian(measurement_pack.raw_measurements_);      
      previous_timestamp_ = measurement_pack.timestamp_;
      ekf_.Init(x, P_, F_, H_laser_, R_radar_, Q);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      VectorXd x(4);
      x << measurement_pack.raw_measurements_[0], 
              measurement_pack.raw_measurements_[1], 
              0, 
              0;

      previous_timestamp_ = measurement_pack.timestamp_;
      ekf_.Init(x, P_ , F_, H_laser_, R_laser_, Q);
    }

    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */

   // compute the time elapsed between the current and previous measurements
  // dt - expressed in seconds
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  // Update the state transition matrix F  and  noise covariance matrix
  // according to the new elapsed time.
  ekf_.UpdateStateMatrices(dt, noise_ax, noise_ay);
  ekf_.Predict();

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    ekf_.UpdateEKF(measurement_pack.raw_measurements_, R_radar_);
  } else {
	ekf_.Update(measurement_pack.raw_measurements_, R_laser_);
  }

  // print the output
  cout << "x_ = " << ekf_.GetEstimatedPos() << endl;
  cout << "P_ = " << ekf_.GetCovariance() << endl;
}
