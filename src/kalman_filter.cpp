#include "kalman_filter.h"
#define PI 3.14159265

using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &Hj_in, MatrixXd &R_l_in,
                        MatrixXd &R_r_in, MatrixXd &Q_in) {
  // Kalman Filter variables
  x_ = x_in;  // object state
  P_ = P_in;  // object covariance matrix
  F_ = F_in;  // state transition matrix
  H_ = H_in;  // measurement matrix, normal KF
  Hj_ = Hj_in;  // measurement matrix, extended KF
  R_l_ = R_l_in;  // measurement covariance matrix, laser
  R_r_ = R_r_in;  // measurement covariance matrix, radar
  I_ = MatrixXd::Identity(x_.size(), x_.size());  // identity matrix
  Q_ = Q_in;  // process covariance matrix
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
  VectorXd y = z - H_ * x_;
  
  // Rest of the update cycle shared with the EKF
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_l_;
  MatrixXd Si = S.inverse();
  MatrixXd K =  P_ * Ht * Si;
  
  // New state
  x_ = x_ + (K * y);
  P_ = (I_ - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) { 
  // For radar measurements, the functions that map the x vector [px, py, vx, vy] to polar coordinates are non-linear.
  // y = z - h(x') instead of y = z - H_ * x_
  double px = x_(0);
  double py = x_(1);
  double vx = x_(2);
  double vy = x_(3);
  double h1 = sqrt(px*px + py*py);
  // Check to avoid dividing by zero
  if (h1 < 0.00001)  return;
  double h2 = atan2(py, px);
  double h3 = (px*vx + py*vy) / h1;  
  VectorXd h = VectorXd(3);
  h << h1, h2, h3;
  VectorXd y = z - h;
  
  // Normalizing angles
  if( y(1) > PI )   y(1) -= 2.f*PI;
  if( y(1) < -PI )  y(1) += 2.f*PI;
    
  // Rest of the update cycle shared with the normal KF
  // H must be the Jacobian matrix in this case!
  Hj_ = tools.CalculateJacobian(x_);
  MatrixXd Hjt = Hj_.transpose();
  MatrixXd S = Hj_ * P_ * Hjt + R_r_;
  MatrixXd Si = S.inverse();
  MatrixXd K =  P_ * Hjt * Si;
  
  // New state
  x_ = x_ + (K * y);
  P_ = (I_ - K * Hj_) * P_;
}
