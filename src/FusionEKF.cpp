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
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  // measurement matrix - laser
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;
  
  // measurement matrix - radar (Hj)
  // not constant
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    // Initializing object covariance matrix
    ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_ << 1, 0, 0, 0,
               0, 1, 0, 0,
               0, 0, 1000, 0,
               0, 0, 0, 1000;

    // Initialize state with first measurement
    ekf_.x_ = VectorXd(4);
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      double rho = measurement_pack.raw_measurements_[0]; // range
      double phi = measurement_pack.raw_measurements_[1]; // bearing
      double rho_dot = measurement_pack.raw_measurements_[2]; // velocity of rho
      
      // Coordinates convertion from polar to cartesian
      double x = rho * cos(phi);
      if ( x < 0.0001 )  x = 0.0001;
      double y = rho * sin(phi);
      if ( y < 0.0001 )  y = 0.0001;
      double vx = rho_dot * cos(phi);
      double vy = rho_dot * sin(phi);
      
      ekf_.x_ << x, y, vx , vy;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */
  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  
  // Update state transition matrix
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, dt, 0,
             0, 1, 0, dt,
             0, 0, 1, 0,
             0, 0, 0, 1;
  
  // Update noise covariance matrix with time difference
  double noise_ax = 9.0;  // Noise values are given
  double noise_ay = 9.0;  // Noise values are given
  double dt2 = dt * dt;
  double dt3 = dt2 * dt;
  double dt4 = dt3 * dt;
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << dt4 / 4 * noise_ax, 0, dt3 / 2 * noise_ax, 0,
	           0, dt4 / 4 * noise_ay, 0, dt3 / 2 * noise_ay,
	           dt3 / 2 * noise_ax, 0, dt2 * noise_ax, 0,
 	           0, dt3 / 2 * noise_ay, 0, dt2 * noise_ay;
  
  ekf_.Predict();

  /**
   * Update
   */
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    ekf_.R_ = R_radar_;
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
  }
  else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
    ekf_.R_ = R_laser_;
    ekf_.H_ = H_laser_;
  }
  ekf_.Update(measurement_pack.raw_measurements_);

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
