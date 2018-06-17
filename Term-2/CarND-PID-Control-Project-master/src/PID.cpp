#include "PID.h"

using namespace std;

/*
* TODO: Complete the PID class.
*/

PID::PID() {}

PID::~PID() {}

void PID::Init(double Kp_In, double Ki_In, double Kd_In) {
	Kp = Kp_In;
	Ki = Ki_In;
	Kd = Kd_In;
	
	p_error = 0.0;
	i_error = 0.0;
	d_error = 0.0;
}

void PID::UpdateError(double cte) {
	double pre_cte = p_error;

	p_error = cte;
	i_error += cte;
	d_error = cte - pre_cte;
}

double PID::TotalError() {
	return 0;	
}

double PID::OutputSteerAng() {
  return -Kp*p_error - Ki*i_error - Kd*d_error;
}

double PID::OutputThrottle(double max_thro){
  return max_thro - Kp*p_error - Ki*i_error - Kd*d_error;
}
