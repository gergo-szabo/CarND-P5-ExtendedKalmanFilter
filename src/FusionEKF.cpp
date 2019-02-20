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
    MatrixXd P_init(4, 4);
    P_init << 1, 0, 0, 0,
              0, 1, 0, 0,
              0, 0, 1000, 0,
              0, 0, 0, 1000;
	  
    // Initializing state transition matrix
    // values will be updated in prediction step
    MatrixXd F_init(4, 4);
    F_init << 1, 0, 0, 0,
	      0, 1, 0, 0,
              0, 0, 1, 0,
              0, 0, 0, 1;
	  
    // Initializing state transition matrix
    // values will be updated in prediction step
    MatrixXd Q_init(4, 4);
    Q_init << 0, 0, 0, 0,
	      0, 0, 0, 0,
              0, 0, 0, 0,
              0, 0, 0, 0;
	  
    // Initialize Kalman filter with first measurement (state)
    VectorXd x_init(4);
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      float rho = measurement_pack.raw_measurements_[0]; // range
      float phi = measurement_pack.raw_measurements_[1]; // bearing
      float rho_dot = measurement_pack.raw_measurements_[2]; // velocity of rho
      
      // Coordinates convertion from polar to cartesian
      float px = rho * cos(phi);
      float py = rho * sin(phi);
      float vx = rho_dot * cos(phi);
      float vy = rho_dot * sin(phi);
      x_init << px, py, vx , vy;  
      ekf_.Init(x_init, P_init, F_init, H_laser_, Hj_, R_laser_, R_radar_, Q_init);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      x_init << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
      ekf_.Init(x_init, P_init, F_init, H_laser_, Hj_, R_laser_, R_radar_, Q_init);
    }
	
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  
  // Update state transition matrix
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, dt, 0,
             0, 1, 0, dt,
             0, 0, 1, 0,
             0, 0, 0, 1;
	
  // Update noise covariance matrix with time difference
  float noise_ax = 9.0;  // Noise values are given
  float noise_ay = 9.0;  // Noise values are given
  float dt2 = dt * dt;
  float dt3 = dt2 * dt;
  float dt4 = dt3 * dt;
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << dt4/4*noise_ax, 0,              dt3/2*noise_ax, 0,
	     0,              dt4/4*noise_ay, 0,              dt3/2*noise_ay,
	     dt3/2*noise_ax, 0,              dt2*noise_ax,   0,
 	     0,              dt3/2*noise_ay, 0,              dt2*noise_ay;
  
  ekf_.Predict();

  /**
   * Update
   */
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  }
  else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
    ekf_.Update(measurement_pack.raw_measurements_);
  }
	
  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
