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

  /**
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   */
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;
  ekf_.H_ = H_laser_;


}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  VectorXd z = measurement_pack.raw_measurements_;
  if (!is_initialized_) {
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    // first measurement
    ekf_.x_ = VectorXd(4);
    ekf_.P_ = MatrixXd::Zero(4, 4);
    previous_timestamp_ = measurement_pack.timestamp_;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      ekf_.x_(0) = z(0) * std::cos(z(1));
      ekf_.x_(1) = z(0) * std::sin(z(1));
      ekf_.x_(2) = 0;//z(2) * std::cos(z(1));
      ekf_.x_(3) = 0;//z(2) * std::sin(z(1));
      // We assume the uncertainty is proportional to the sensor's distance
      // uncertainty
      ekf_.P_(0, 0) = 1; //R_radar_(0, 0);
      ekf_.P_(1, 1) = 1; //_radar_(0, 0);
      ekf_.P_(2, 2) = 3; // Prior - root of the average velocity of a bicycle 
      ekf_.P_(3, 3) = 3; // on the road, in m/s

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      ekf_.x_(0) = z(0);
      ekf_.x_(1) = z(1);
      ekf_.x_(2) = 0; // We assume a velocity of uncertain magnitude in any
      ekf_.x_(3) = 0; // direction
      ekf_.P_(0, 0) = 1; //R_laser_(0, 0);
      ekf_.P_(1, 1) = 1; //R_laser_(1, 1);
      ekf_.P_(2, 2) = 3; //R_laser_(0, 0);
      ekf_.P_(3, 3) = 3; //R_laser_(1, 1);
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */
  double dt = double(measurement_pack.timestamp_ - previous_timestamp_)/1e6;
  if (dt > 0) {
    ekf_.F_ = MatrixXd(4, 4);
    ekf_.F_ << 1, 0, dt, 0,
                0, 1, 0,  dt,
                0, 0, 1,  0,
                0, 0, 0,  1;
    double dt4 = std::pow(dt, 4);
    double dt3 = std::pow(dt, 3);
    double dt2 = std::pow(dt, 2);
    double s_a = 9; // Acceleration noise variance
    ekf_.Q_ = MatrixXd(4, 4);
    ekf_.Q_ << dt4/4, 0, dt3/2, 0,
                0, dt4/4, 0, dt3/2,
                dt3/2, 0, dt2, 0,
                0, dt3/2, 0, dt2;
    ekf_.Q_ *= s_a;

    /**
    * TODO: Update the state transition matrix F according to the new elapsed time.
    * Time is measured in seconds.
    * TODO: Update the process noise covariance matrix.
    * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
    */

    ekf_.Predict();
  }

  /**
   * Update
   */

  /**
   * TODO:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // TODO: Radar updates
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(z);

  } else {
    // TODO: Laser updates
    ekf_.R_ = R_laser_;
    ekf_.Update(z);

  }
}
