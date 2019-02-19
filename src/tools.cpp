#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
	VectorXd rmse(4);
	rmse << 0,0,0,0;
	
	// Sanity check for inputs:
	//  * the ground truth vector size should not be zero
	//  * the estimation vector size should not be zero
	//  * the estimation vector size should equal ground truth vector size
	if (ground_truth.size() == 0) {
		std::cout << "ERROR: CalculateRMSE() - The ground-truth vector is empty" << std::endl;
		return rmse;
 	}
	if (estimations.size() == 0) {
		std::cout << "ERROR: CalculateRMSE() - The estimations vector is empty" << std::endl;
		return rmse;
 	}
	if (estimations.size() != ground_truth.size()) {
		std::cout << "ERROR: CalculateRMSE() - The ground-truth and estimations vectors are not the same size." << std::endl;
		return rmse;
	}
	
	// Calculate RMSE
	for (unsigned int i=0; i < estimations.size(); ++i) {
		VectorXd residual = estimations[i] - ground_truth[i];
		// coefficient-wise multiplication
		residual = residual.array()*residual.array();
		rmse += residual;
	}
	rmse = rmse/estimations.size();
	rmse = rmse.array().sqrt();

	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
	MatrixXd Hj(3,4);
	
	// Sanity check for inputs:
	if ( x_state.size() != 4 ) {
		std::cout << "ERROR: CalculateJacobian() - The state vector size is not 4." << std::endl;
		return Hj;
	}
	
	// recover state parameters
	double px = x_state(0);
	double py = x_state(1);
	double vx = x_state(2);
	double vy = x_state(3);

	// pre-compute a set of terms to avoid repeated calculation
	double c1 = px*px+py*py;
	double c2 = sqrt(c1);
	double c3 = (c1*c2);

	// check division by zero
	if (fabs(c1) < 0.0001) {
		std::cout << "ERROR: CalculateJacobian() - Division by Zero" << std::endl;
		return Hj;
	}

	// compute the Jacobian matrix
	Hj << (px/c2), (py/c2), 0, 0,
		  -(py/c1), (px/c1), 0, 0,
		  py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

	return Hj;
}
