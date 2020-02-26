#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

double normalize_angle(double theta) {
    while (theta > M_PI) 
        theta -= 2*M_PI;
    while (theta < -M_PI)
        theta += 2*M_PI;
    return theta;
}
double angle_difference (double a, double b) {
    return normalize_angle(normalize_angle(a) - normalize_angle(b));
}

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
}

void KalmanFilter::Predict() {
  x_ = F_*x_;
  P_ = F_*P_*F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  Eigen::VectorXd y = z - H_*x_;
  Eigen::MatrixXd S = H_*P_*H_.transpose() + R_;
  Eigen::MatrixXd K = P_*H_.transpose()*S.inverse();
  x_ = x_ + K*y;
  P_ = (Eigen::MatrixXd::Identity(4, 4) - K*H_)*P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  double x  = x_(0);
  double y  = x_(1);
  double vx = x_(2);
  double vy = x_(3);
  double r2 = x*x + y*y;
  if (r2 < 0.001) {
      std::cout << "Range too small to update - div by zero\n";
      return;
  }
  double r = std::sqrt(r2);
  double r3 = r2*r;
  Eigen::MatrixXd H(3, 4);
  H << x/r , y/r, 0, 0,
       -y/r2, x/r2, 0, 0,
       y*(vx*y - vy*x)/r3, x*(x*vy - y*vx)/r3, x/r, y/r;
  Eigen::VectorXd z_p(3);
  z_p << r, normalize_angle(std::atan2(y, x)), (x*vx + y*vy)/r;
  // Wrap angle
  Eigen::VectorXd innov = z - z_p;
  innov(1) = angle_difference(z(1), z_p(1));
  Eigen::MatrixXd S = H*P_*H.transpose() + R_;
  Eigen::MatrixXd K = P_*H.transpose()*S.inverse();
  x_ = x_ + K*innov;
  P_ = P_ - K*H*P_;
}
