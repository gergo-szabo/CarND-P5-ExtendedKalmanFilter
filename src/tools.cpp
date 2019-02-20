#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() : residual(4), lastSum(4), rmse(4) {
	lastSum << 0, 0, 0, 0;
}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
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
	residual = estimations[estimations.size()-1] - ground_truth[estimations.size()-1];
	residual = residual.array()*residual.array();	// coefficient-wise multiplication
	lastSum += residual;
	rmse = lastSum/estimations.size();
	rmse = rmse.array().sqrt();

	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
	float px0 = x_state(0);
	float py0 = x_state(1);
	float vx0 = x_state(2);
	float vy0 = x_state(3);
   
	// pre-compute a set of terms to avoid repeated calculation
	float c1 = px0*px0+py0*py0;
	float c2 = sqrt(c1);
	float c3 = (c1*c2);
	
	// check before division by zero
	if (c1 == 0) {
		std::cout << "ERROR: CalculateJacobian() - Division by Zero" << std::endl;
		return Hj;
	}
	
	// compute the Jacobian matrix
	MatrixXd Hj(3,4);
	Hj << (px0/c2), (py0/c2), 0, 0,
		  -(py0/c1), (px0/c1), 0, 0,
		  py0*(vx0*py0 - vy0*px0)/c3, px0*(px0*vy0 - py0*vx0)/c3, px0/c2, py0/c2;
  
	return Hj;
}
