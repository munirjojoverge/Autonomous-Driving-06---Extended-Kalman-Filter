#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#define epsilon 0.0001 // A small number

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  /*****************************************************************************
     initializing matrices
  ****************************************************************************/
  
  //create a 4D state vector, we don't know yet the values of the x state
  ekf_.x_ = VectorXd(4);  
  
  H_laser_ = MatrixXd(2, 4);
  H_laser_ << 1, 0, 0, 0,
	          0, 1, 0, 0;
  H_radar_ = MatrixXd(3, 4); // This will hold the Jacobian of x_ every update step

  
  //measurement covariance matrix - laser
  R_laser_ = MatrixXd(2, 2);
  R_laser_ << 0.0225, 0,
              0,      0.0225;

  //measurement covariance matrix - radar
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << 0.09, 0,      0,
              0,    0.0009, 0,
              0,    0,      0.09;
    
  //state covariance matrix P
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
			 0, 1, 0, 0,
			 0, 0, 1000, 0,
	         0, 0, 0, 1000;

  //the initial transition matrix F_
  ekf_.F_ = MatrixXd(4, 4);

  // noise covariance matrix Q
  ekf_.Q_ = MatrixXd(4, 4);

  //set the acceleration noise components
  noise_ax_ = 7;
  noise_ay_ = 7;
  
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

/**
* Process Measurement.
*/
void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
	
	Eigen::VectorXd z = measurement_pack.raw_measurements_;
		
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */	  
	float px;
	float py;
	float vx;
	float vy;

	//cout << "Z(0) = " << z(0) << endl;
	//cout << "Z(1) = " << z(1) << endl;
	

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
		float rho     = z(0);
		float phi     = z(1);
		float rho_dot = z(3);

		px = rho * cos(phi);
		py = rho * sin(phi);
		vx = rho_dot * cos(phi);
		vy = rho_dot * sin(phi);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
		px = z(0);
		py = z(1);
		vx = 0.0;
		vy = 0.0;
    }
	// Assign the initial state
	ekf_.x_ << px, py, vx, vy;
	// Deal with the special case where we start at 0, 0 and is a Radar measurment. This way we don't wait for the Jacobian to plot the error.
	if (fabs(ekf_.x_(0)) < epsilon && fabs(ekf_.x_(1)) < epsilon) {
		ekf_.x_(0) = epsilon;
		ekf_.x_(1) = epsilon;
	}
	//cout << "Initial State = " << ekf_.x_ << endl;

	previous_timestamp_ = measurement_pack.timestamp_;
    
	// done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
   //compute the time elapsed between the current and previous measurements
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;

  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;

  //Modify the F matrix so that the time is integrated
  ekf_.F_ << 1, 0, dt, 0,
		     0, 1, 0, dt,
	         0, 0, 1, 0,
	         0, 0, 0, 1;
  //cout << "F = " << ekf_.F_ << endl;

  //set the process covariance matrix Q  
  ekf_.Q_ << dt_4 / 4 * noise_ax_, 0, dt_3 / 2 * noise_ax_, 0,
	         0, dt_4 / 4 * noise_ay_, 0, dt_3 / 2 * noise_ay_,
	         dt_3 / 2 * noise_ax_, 0, dt_2*noise_ax_, 0,
	         0, dt_3 / 2 * noise_ay_, 0, dt_2*noise_ay_;

  //predict
  ekf_.Predict();
  //cout << "State After Prediction = " << ekf_.x_ << endl;

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
	  H_radar_ = tools.CalculateJacobian(ekf_.x_);	  
	  ekf_.Init(ekf_.x_, ekf_.P_, ekf_.F_, H_radar_, R_radar_, ekf_.Q_);
	  ekf_.UpdateEKF(z);
  } else {
    // Laser updates	  
	  ekf_.Init(ekf_.x_, ekf_.P_, ekf_.F_, H_laser_, R_laser_, ekf_.Q_);
	  ekf_.Update(z);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
