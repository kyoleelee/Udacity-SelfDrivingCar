#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);
	x_.setZero();

  // initial covariance matrix
  P_ = MatrixXd(5, 5);
	P_.setZero();

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1;   //need to change, try 1

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1;    //need to change try 1
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...

  */

	is_initialized_ = false;
	time_us_ = 0;

	n_x_ = 5;
	n_aug_ = 7;
	n_aug_col_ = 2*n_aug_ + 1;
	lambda_ = 3 - n_aug_;

	weights_ = VectorXd(n_aug_col_);
  double weight_0 = lambda_ / (lambda_ + n_aug_);
  double weight_x = 0.5 / (lambda_ + n_aug_);
  weights_(0) = weight_0;
  for (int i = 1; i < n_aug_col_; i++) { // 2n+1 weights
    weights_(i) = weight_x;
	}
	
	Xsig_pred_ = MatrixXd(n_x_, n_aug_col_);
	Xsig_pred_.setZero();


  // Initial NIS values for each sensor type
  NIS_laser_ = 0.0;
	NIS_radar_ = 0.0;


}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  
  // Reset filter if timestamp resets (restart dataset) or
  // has a gap > 1000 sec (switch dataset)
  if ((meas_package.timestamp_ < time_us_) ||
      (abs(time_us_ - meas_package.timestamp_) > 1000000000.0)) {
    is_initialized_ = false;
    NIS_laser_ = 0.0;
    NIS_radar_ = 0.0;
	}

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    time_us_ = meas_package.timestamp_ ;
	
    // Initialize Covariance Matrix P in lecture: "What to Expect from the Project"
    P_ <<
					1.0, 0.0, 0.0, 0.0, 0.0,
					0.0, 1.0, 0.0, 0.0, 0.0,
					0.0, 0.0, 1.0, 0.0, 0.0,
					0.0, 0.0, 0.0, 1.0, 0.0,
					0.0, 0.0, 0.0, 0.0, 1.0;
	
    // Initialize state vector
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      double rho = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];
			x_ << (rho * cos(phi)), (rho * sin(phi)), 0, 0, 0;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
			//set the state with the initial location and zero velocity
			x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
		  }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }
  
  /*****************************************************************************
  *  Prediction
   ****************************************************************************/

  double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0; 
  time_us_ = meas_package.timestamp_;
	
	//UKF predict with current dt
	Prediction(delta_t);
  

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    UpdateRadar(meas_package);
  } else {
    // Laser updates
	UpdateLidar(meas_package);
  }
}
  

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
	

 	//Step 1 - Choose representative augmented sigma points
  
  
  //create augmented mean vector, state covariance, and sigma point matrix
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.setZero();
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.setZero();
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_aug_col_);
  Xsig_aug.setZero();
  
  //set augmented mean state
  x_aug.head(n_x_) = x_;
  x_aug(n_x_) = 0;
  x_aug(n_x_+1) = 0;
  
  //set augmented covariance matrix
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_x_, n_x_) = std_a_*std_a_; // acceleration covariance
  P_aug(n_x_+1, n_x_+1) = std_yawdd_*std_yawdd_; // yaw_dot_dot covariance
  
  //calculate square root of P_aug
  MatrixXd P_aug_sqrt = P_aug.llt().matrixL();
  
  //set augmented sigma point matrix Xsig_aug
  Xsig_aug.col(0)  = x_aug;

  for (int i = 0; i< n_aug_; i++) {
    Xsig_aug.col(i+1) = x_aug + sqrt(lambda_+n_aug_) * P_aug_sqrt.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * P_aug_sqrt.col(i);
  }
  

	//Step 2 - Predict sigma points to current time step


  //loop through each augmented sigma point column to set Xsig_pred_
  for (int i = 0; i < n_aug_col_; i++) {
    
    // Extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);
    
    // Calculate base predicted state values
    double px_p, py_p;
    
    // Switch motion equation to avoid division by zero from yawd
    if (fabs(yawd) > 0.001) {
      // Driving in a curve
      px_p = p_x + v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
      py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
    }
    else {
      // Driving straight
      px_p = p_x + v * delta_t * cos(yaw);
      py_p = p_y + v * delta_t * sin(yaw);
    }
    double v_p = v; // constant velocity
    double yaw_p = yaw + yawd * delta_t;
    double yawd_p = yawd; // constant turn rate
    
    // Add process noise to predicted state values
    px_p = px_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
    py_p = py_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);
    v_p = v_p + nu_a * delta_t;
    yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t * delta_t;
    yawd_p = yawd_p + nu_yawdd * delta_t;
    
    // Write predicted sigma point into Xsig_pred_ column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }
  

  //Step 3 - Calculate predicted mean state and covariance

  
  // Set predicted state mean from weighted predicted sigma points
  x_.fill(0.0);

  for (int i = 0; i < n_aug_col_; i++) { // iterate over sigma points
    x_ = x_+ weights_(i) * Xsig_pred_.col(i);
  }
  
  // Set predicted state covariance matrix from weighted predicted sigma points
  P_.fill(0.0);

  for (int i = 0; i < n_aug_col_; i++) { // iterate over sigma points
    // State difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    // Angle normalization between [-pi, pi] for yaw
    while (x_diff(3) > M_PI) { x_diff(3) -= 2.*M_PI; }
    while (x_diff(3) < -M_PI) { x_diff(3) += 2.*M_PI; }
    
    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
	}

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
	
	VectorXd z = meas_package.raw_measurements_; //incoming lidar measurement
	
	MatrixXd H_ = MatrixXd(2, 5);
	H_ << 1, 0, 0, 0, 0,
				0, 1, 0, 0, 0;

	MatrixXd R_ = MatrixXd(2, 2);
	R_ << std_laspx_*std_laspx_, 0, 0, std_laspy_*std_laspy_;

	VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;

  //calculate LASER NIS for consistency check
	NIS_laser_ = y.transpose() * Si * y;
 	cout << "LASER NIS: " << NIS_laser_ << endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  
	
	//Step 1 - Transform sigma points to measurement space
	

	//get covariance matrix S
	int n_z = 3;
	MatrixXd Zsig = MatrixXd(n_z, n_aug_col_);
	Zsig.setZero();

	//transform predicted sigma points to measurement space
	for (int i = 0; i < n_aug_col_; i++) {  //2n+1 simga points
		//extract values for better readibility
		double px = Xsig_pred_(0,i);
		double py = Xsig_pred_(1,i);
		double v  = Xsig_pred_(2,i);
		double yaw = Xsig_pred_(3,i);
		double vx = cos(yaw)*v;
		double vy = sin(yaw)*v;
		
		//measurement model
		Zsig(0,i) = sqrt(px*px + py*py);												//r
		Zsig(1,i) = atan2(py,px);   														//phi
		Zsig(2,i) = (px*vx + py*vy) / sqrt(px*px + py*py);   		//r_dot
	}


	//Step 2 - Calculate predicted measurement mean, covariance


	//predicted measurement mean from weighted predicted measurement 
	VectorXd z_pred = VectorXd(n_z);
	z_pred.fill(0.0);
	for (int i=0; i < n_aug_col_; i++) {
		z_pred = z_pred + weights_(i) * Zsig.col(i);
	}

	//measurement covariance matrix S
	MatrixXd S = MatrixXd(n_z,n_z);
	S.fill(0.0);
	for (int i = 0; i < n_aug_col_; i++) {  //2n+1 simga points
			
		//residual
		VectorXd z_diff = Zsig.col(i) - z_pred;

		//angle normalization
		while (z_diff(1)> M_PI) { z_diff(1)-=2.*M_PI; }
		while (z_diff(1)<-M_PI) { z_diff(1)+=2.*M_PI; }
		S = S + weights_(i) * z_diff * z_diff.transpose();
	}

	//add measurement noise covariance matrix
	MatrixXd R = MatrixXd(n_z,n_z);
	R <<    (std_radr_ * std_radr_), 0, 0,
					0, (std_radphi_ * std_radphi_), 0,
					0, 0, (std_radrd_ * std_radrd_);
	S = S + R;


	//Step 3 - Calculate Kalman Gain using cross-covariance


	//create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd(n_x_, n_z);
	Tc.fill(0.0);

	for (int i = 0; i < n_aug_col_; i++) {  //2n+1 simga points
		//residual
		VectorXd z_diff = Zsig.col(i) - z_pred;

		//angle normalization
		while (z_diff(1)> M_PI) { z_diff(1)-=2.*M_PI; }
		while (z_diff(1)<-M_PI) { z_diff(1)+=2.*M_PI; }

		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;

		//angle normalization
		while (x_diff(3)> M_PI) { x_diff(3)-=2.*M_PI; }
		while (x_diff(3)<-M_PI) { x_diff(3)+=2.*M_PI; }

		//cross-covariance matrix Tc
		Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
	}

	//calculate Kalman gain K;
	MatrixXd K = Tc * S.inverse();


	//Step 4 - Update state and covariance with Kalman gain


	VectorXd z = meas_package.raw_measurements_; //incoming radar measurement

	//residual between actual and predicted measurement
	VectorXd z_diff = z - z_pred;

	//angle normalization
	while (z_diff(1)> M_PI) { z_diff(1)-=2.*M_PI; }
	while (z_diff(1)<-M_PI) { z_diff(1)+=2.*M_PI; }

	//update state mean and covariance matrix
	x_ = x_ + K * z_diff;
	P_ = P_ - K*S*K.transpose();  
  
  //calculate RADAR NIS for consistency check
  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
	cout << "RADAR NIS: " << NIS_radar_ << endl;

  
}
