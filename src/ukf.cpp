#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

double unwrap2pi(double angle) {
  double x = cos(angle);
  double y = sin(angle);
  return atan2f(y,x);
}

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
  x_.fill(0.0);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  P_.fill(0.0);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = .5; // .5

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = .3; // 0.3 being ~ pi/4 /2

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

  is_initialized_ = false;

  n_x_ = 5;

  n_aug_ = 7;

  lambda_ = 3 - n_aug_;

  P_ = MatrixXd(n_x_, n_x_);
  P_.fill(0.0);

  weights_ = VectorXd(2*n_aug_ +1);
  weights_.fill(0.0);

  //create matrix with predicted sigma points as columns
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  Xsig_pred_.fill(0.0);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  if (!is_initialized_) {

    is_initialized_ = true;

    if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      double x,y,rho;
      x = meas_package.raw_measurements_(0);
      y = meas_package.raw_measurements_(1);
      rho = atan2f(y,x);
      x_(0) = x;
      x_(1) = y;
    } else {
      // RADAR measurement
      double rho, phi, x, y;
      rho = meas_package.raw_measurements_(0);
      phi = meas_package.raw_measurements_(1);
      x = cos(phi) * rho;
      y = sin(phi) * rho;
      x_(0) = x;
      x_(1) = y;
    }

    // Initialize the covariance matrix
    P_(0,0) = 1.;
    P_(1,1) = 1.;
    P_(2,2) = 6.;
    P_(3,3) = 3.;
    P_(4,4) = 6.;

    time_us_ = meas_package.timestamp_;

    return;
  }

  double delta_t = (meas_package.timestamp_ - time_us_ ) / 1000000.;
  time_us_ = meas_package.timestamp_;

  // std::cout << "Predicting next pose" << std::endl;
  Prediction(delta_t);

  // std::cout << "Predicted next pose" << std::endl;
  if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    // std::cout << "Updating lidar" << std::endl;
    if (use_laser_) UpdateLidar(meas_package);
  } else {
    // std::cout << "Updating radar" << std::endl;
    if (use_radar_) UpdateRadar(meas_package);
    // std::cout << "Updated radar" << std::endl;
  }

  x_(3) = unwrap2pi(x_(3));


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


  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.fill(0.0);
  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  MatrixXd P_noise = MatrixXd(n_aug_-n_x_, n_aug_-n_x_);
  P_noise << std_a_ * std_a_, 0, 0, std_yawdd_ * std_yawdd_;

  //create augmented mean state
  x_aug.head(n_x_) = x_;
  //create augmented covariance matrix
  P_aug.topLeftCorner(n_x_,n_x_) = P_;
  P_aug.bottomRightCorner(n_aug_-n_x_, n_aug_-n_x_) = P_noise;
  //create square root matrix
  MatrixXd A  = P_aug.llt().matrixL();

  // std::cout << "A matrix \n" << A << std::endl;

  //create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  for (int col_n = 0; col_n < n_aug_; col_n++) {
      Xsig_aug.col(1+col_n) = x_aug + sqrt(lambda_ + n_aug_) * A.col(col_n);
      Xsig_aug.col(1+n_aug_+col_n) = x_aug - sqrt(lambda_ + n_aug_) * A.col(col_n);
      // Xsig_aug(3,1+col_n) = unwrap2pi(Xsig_aug(3,1+col_n));
      // Xsig_aug(3,1+n_aug_+col_n) = unwrap2pi(Xsig_aug(3,1+n_aug_+col_n));
  }

  //predict sigma points
  for (int i = 0; i < Xsig_aug.cols(); i++ ) {
      double x,y,v,phi,phi_dot,epsilon_acc, epsilon_ang_acc, delta_t2;
      delta_t2 = delta_t * delta_t;
      x = Xsig_aug(0,i);
      y = Xsig_aug(1,i);
      v = Xsig_aug(2,i);
      phi = Xsig_aug(3,i);
      phi_dot = Xsig_aug(4,i);
      epsilon_acc = Xsig_aug(5,i);
      epsilon_ang_acc = Xsig_aug(6,i);

      Xsig_pred_(0,i) = x + 0.5 * delta_t2 * epsilon_acc * cos(phi);
      Xsig_pred_(1,i) = y + 0.5 * delta_t2 * epsilon_acc * sin(phi);

      if (!(fabs(phi_dot) < 0.001)) {
        Xsig_pred_(0,i) += v / phi_dot * (sin(phi + phi_dot * delta_t) - sin(phi));
        Xsig_pred_(1,i) += v / phi_dot * (-cos(phi + phi_dot * delta_t) + cos(phi));
      } else {
        Xsig_pred_(0,i) += v * delta_t * cos(phi);
        Xsig_pred_(1,i) += v * delta_t * sin(phi);
      }

      Xsig_pred_(2,i) = v + epsilon_acc * delta_t;
      Xsig_pred_(3,i) = phi + phi_dot * delta_t + 0.5 * epsilon_ang_acc * delta_t2;
      Xsig_pred_(4,i) = phi_dot + epsilon_ang_acc * delta_t;

  }

  //set weights
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i < 2*n_aug_+1; i++) weights_(i) = 1 / (2 * (lambda_ + n_aug_) );
  //predict state mean
  for (int i = 0; i < n_x_; i++) x_(i) = Xsig_pred_.row(i) * weights_;
  x_(3) = unwrap2pi(x_(3));

  // std::cout << "Projected forward Xsig_pred_\n " <<  Xsig_pred_ << std::endl;


  //predict state covariance matrix
  P_.fill(0.0); // reset P
  for (int i = 0; i < 2*n_aug_+1; i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    x_diff(3) = unwrap2pi(x_diff(3));
    P_ += x_diff * x_diff.transpose() * weights_(i);
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
  int n_z = 2;
  MatrixXd S = MatrixXd(n_z,n_z);
  // fill in the measurement device noise
  S << std_laspx_ * std_laspx_, 0,
       0, std_laspy_ * std_laspy_;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  Zsig.fill(0.0);
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  VectorXd z_diff = VectorXd(n_z);
  z_diff.fill(0.0);

  //transform sigma points into measurement space
  for (int i=0; i < 2 * n_aug_ + 1; i++) {
     double x,y, v, phi;
     x = Xsig_pred_(0,i);
     y = Xsig_pred_(1,i);

     Zsig(0,i) = x;
     Zsig(1,i) = y;
  }
  //calculate mean predicted measurement
  for (int i = 0; i < n_z; i++) {
   z_pred(i) = Zsig.row(i) * weights_;
  }
  //calculate innovation covariance matrix S
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
   z_diff = z_pred - Zsig.col(i);
   S += z_diff * z_diff.transpose() * weights_(i);
  }

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.setZero(n_x_, n_z);

  //calculate cross correlation matrix
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
      VectorXd x_diff = (Xsig_pred_.col(i) - x_);
      x_diff(3) = unwrap2pi(x_diff(3));
      z_diff = Zsig.col(i) - z_pred;
      Tc += weights_(i) * x_diff * z_diff.transpose() ;
  }

  // std::cout << "computing kalman update" << std::endl;

  //calculate Kalman gain K;
  MatrixXd S_inv = S.inverse();
  MatrixXd K = Tc * S_inv;

  //update state mean and covariance matrix
  z_diff = meas_package.raw_measurements_ - z_pred;

  // std::cout << "NIS value: " << z_diff.transpose() * S_inv * z_diff << std::endl;

  x_ += K * z_diff;

  P_ += -(K * S * K.transpose());

  // std::cout << "x: \n" << x_ << std::endl;
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

  int n_z = 3;
  MatrixXd S = MatrixXd(n_z,n_z);
  // fill in the measurement device noise
  S << std_radr_ * std_radr_, 0 , 0,
       0, std_radphi_ * std_radphi_, 0,
       0, 0, std_radrd_ * std_radrd_;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  VectorXd z_diff = VectorXd(n_z);

  //transform sigma points into measurement space
  for (int i=0; i < 2 * n_aug_ + 1; i++) {
     double x,y, v, phi;
     x = Xsig_pred_(0,i);
     y = Xsig_pred_(1,i);
     v = Xsig_pred_(2,i);
     phi = Xsig_pred_(3,i);

     Zsig(0,i) = sqrt(x*x + y*y);
     Zsig(1,i) = atan2f(y,x);
     Zsig(2,i) = 1./Zsig(0,i) * (x * cos(phi) * v + y * sin(phi) * v);
  }
  //calculate mean predicted measurement
  for (int i = 0; i < 3; i++) {
   z_pred(i) = Zsig.row(i) * weights_;
  }
  //calculate innovation covariance matrix S

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
   z_diff = z_pred - Zsig.col(i);
   z_diff(1) = unwrap2pi(z_diff(1));
   S += z_diff * z_diff.transpose() * weights_(i);
  }

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.setZero(n_x_, n_z);

  //calculate cross correlation matrix
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
      VectorXd x_diff = (Xsig_pred_.col(i) - x_);
      z_diff = Zsig.col(i) - z_pred;
      // wrap around -pi and pi
      x_diff(3) = unwrap2pi(x_diff(3));
      z_diff(1) = unwrap2pi(z_diff(1));

      Tc += weights_(i) * x_diff * z_diff.transpose() ;
  }

  //calculate Kalman gain K;
  MatrixXd S_inv = S.inverse();
  MatrixXd K = Tc * S_inv;

  //update state mean and covariance matrix
  z_diff = meas_package.raw_measurements_ - z_pred;
  z_diff(1) = unwrap2pi(z_diff(1));

  // std::cout << "NIS value: " << z_diff.transpose() * S_inv * z_diff << std::endl;

  x_ += K * z_diff;

  P_ += -(K * S * K.transpose());


}
