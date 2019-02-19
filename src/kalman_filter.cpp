#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  // Kalman Filter variables
  x_ = x_in;  // object state
  P_ = P_in;  // object covariance matrix
  F_ = F_in;  // state transition matrix
  H_ = H_in; 	// measurement matrix
  R_ = R_in; 	// measurement covariance matrix
  I_ = MatrixXd::Identity(x_.size(), x_.size()); // identity matrix
  Q_ = Q_in; 	// process covariance matrix
}

void KalmanFilter::Predict() {
  // x = F*x + B*u + v
  // B*u = 0 because of no control (B: control input matrix; u: control vector)
  // v = 0 because noise is expected to be gaussian distribution with zero mean
  x_ = F_ * x_ ;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  // TODO: Implement Sanity check
  // H(2,4) and not H(3,4)(Jacobian)
  
  VectorXd y = z - H_ * x_;
  
  // Rest of the update cycle shared with the EKF
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K =  P_ * Ht * Si;
  
  // New state
  x_ = x_ + (K * y);
  P_ = (I_ - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  // TODO: Implement Sanity check
  // H(3,4)(Jacobian) and not H(2,4)
  
  // For radar measurements, the functions that map the x vector [px, py, vx, vy] to polar coordinates are non-linear.
  // y = z - h(x') instead of y = z - H_ * x_
  double px = x_(0);
  double py = x_(1);
  double vx = x_(2);
  double vy = x_(3);
  double h1 = sqrt(px*px + py*py);
  double h2 = atan2(py, px);
  double h3 = (px*vx + py*vy) / h1;  
  VectorXd h = VectorXd(3);
  h << h1, h2, h3;
  VectorXd y = z - h;
  
  // Normalizing angles
  while ( y(1) > M_PI || y(1) < -M_PI ) {
    if ( y(1) > M_PI )  y(1) -= M_PI;
    else                y(1) += M_PI;
  }
    
  // Rest of the update cycle shared with the normal KF
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K =  P_ * Ht * Si;
  
  // New state
  x_ = x_ + (K * y);
  P_ = (I_ - K * H_) * P_;
}
