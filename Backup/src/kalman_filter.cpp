#include <math.h>
#include "kalman_filter.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;


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
  //predict the state
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
  
}

void KalmanFilter::Update(const VectorXd &z) {
  //update the state by using Kalman Filter equations
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  UpdateProcess(y);
  
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  //update the state by using Extended Kalman Filter 
  
  //recover state parameters
  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);

  VectorXd hx(3);

  hx(0) = sqrt(px*px + py*py);
  hx(1) = atan2(py,px);


  //check for low radial distance values
  if(fabs(hx(0))>1.0e-4) {
    hx(2) = (px*vx + py*vy)/hx(0);
  }else {
    hx(2) = (px*vx + py*vy)*1.0e+4;
  }
  
  VectorXd y = z - hx;

  //normalize angle between -pi and pi
  float pi = atan2(1,1)*4;
  if(y(1) < -pi) {
    while(y(1) < -pi){
    y(1) = y(1) + 2*pi;
    }
  }else if(y(1) > pi) {
    while(y(1) > pi){
    y(1) = y(1) - 2*pi;  
    }
  }
  UpdateProcess(y);
 
}

void KalmanFilter::UpdateProcess(const VectorXd &y) {
  
  
  MatrixXd Ht = H_.transpose();
  MatrixXd PHt = P_ * Ht;
  MatrixXd S = H_ * PHt + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K = PHt * Si;

  //update state
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  
  //update state covariance matrix
  P_ = (I - K * H_) * P_;
}
