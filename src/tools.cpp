#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd error_sum(4);
  error_sum << 0, 0, 0, 0;
  if (estimations.size() != ground_truth.size()) {
      std::cerr << "Estimations and ground truth must have the same size\n";
      return error_sum;
  }
  if (estimations.empty()) {
      std::cerr << "Estimations cannot be empty\n";
      return error_sum;
  }
  for (auto it = estimations.begin(), jt = ground_truth.begin();
       it != estimations.end(); ++it, ++jt) {
      error_sum += (*jt - *it).array().pow(2).matrix();
  }
  return (error_sum/estimations.size()).array().sqrt();
}

//MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
//  /**
//   * TODO:
//   * Calculate a Jacobian here.
//   */
//}
