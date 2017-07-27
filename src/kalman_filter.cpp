/**********************************************
* Self-Driving Car Nano-degree - Udacity
*  Created on: July 12, 2017
*      Author: Munir Jojo-Verge
**********************************************/

#include "kalman_filter.h"
#include <iostream>
#define epsilon 0.0001 // A small number

using Eigen::MatrixXd;
using Eigen::VectorXd;
const long double PI = 3.141592653589793238L;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;

  // Let's create only once an Identity Matrix that we will use on the Update step
  int x_size = x_.size();
  I_ = MatrixXd::Identity(x_size, x_size);
 
}

void KalmanFilter::Predict() {
	/*
	* KF Prediction step
	*/
	x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;

}

void KalmanFilter::Update(const VectorXd &z) {
	/*
	* KF Measurement update step
	*/
	VectorXd y = z - H_ * x_;
	//new updated state x_ and covariance matrix P_
	StateUpdate(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
	/*
	* EKF Measurement update step
	* Algorithm extracted from Lesson 5 - 20. EKF Algorithm Generalization
	* Remeber that we use EKF over the RADAR because the model is not linear and we have to linearize
	*/
	float px = x_[0];
	float py = x_[1];
	float vx = x_[2];
	float vy = x_[3];
	float rho; // The range, ρ, is the distance to the pedestrian
	float rho_dot; // The range rate, ​ρ​˙, is the projection of the velocity, v
	float theta; // φ is the angle between ρ and the x direction

	// maps x_ from Cartesian coordinates to polar coordinates
	rho = sqrt(px * px + py * py); 
   
    // Since we are dividing by rho to calculate rho_dot, we have to make sure is not zero
	if (rho < epsilon) {
		std::cerr << "Error while converting vector x_ to polar coordinates: Division by Zero" << std::endl;
	}
	else {		
		rho_dot = (px * vx + py * vy) / rho; 
	}

    // Since we are dividing by px to calculate arcTangent, we have to make sure is not zero
	
	if (fabs(px) < epsilon) {
		std::cerr << "Error while converting vector x_ to polar coordinates: Division by Zero" << std::endl;
	}
	else {
		theta = atan2(py, px); 
	}
		
	Eigen::VectorXd hx = Eigen::VectorXd(3);
	hx << rho, theta, rho_dot;
	Eigen::VectorXd y = z - hx;	
	if (y(1)>PI) {
		y(1) -= 2 * PI;
	}
	if (y(1)<-PI) {
		y(1) += 2 * PI;
	}

	//new updated state x_ and covariance matrix P_
	StateUpdate(y);
}

// Common STATE Update. Valid for KF and EKF
void KalmanFilter::StateUpdate(const VectorXd &y) {
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd K = P_ * Ht * Si;
	
	// New state
	x_ = x_ + (K * y);	
	P_ = (I_ - K * H_) * P_;
}
