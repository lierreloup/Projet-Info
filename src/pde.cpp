#ifndef __PDE_CPP
#define __PDE_CPP

#include "../include/pde.h"
#include <math.h>
#include <iostream>

BlackScholesPDE::BlackScholesPDE(VanillaOption* _option) : option(_option) {}

// Diffusion coefficient
double BlackScholesPDE::diff_coeff(double t, double x) const {
  double vol = option->sigma;
  return 0.5*vol*vol*x*x;  // \frac{1}{2} \sigma^2 S^2
}

// Convection coefficient
double BlackScholesPDE::conv_coeff(double t, double x) const {
  return (option->r)*x;  // rS
}

// Zero-term coefficient
double BlackScholesPDE::zero_coeff(double t, double x) const {
  return -(option->r);  // -r
}

// Source coefficient
double BlackScholesPDE::source_coeff(double t, double x) const {
  return 0.0;
}

// Left boundary-condition 
double BlackScholesPDE::boundary_left(double t, double x) const {
  return this->option->option_price_for_0_spot(t);
  //return 0.0;  // Specifically for a CALL option
}

// Right boundary-condition 
double BlackScholesPDE::boundary_right(double t, double x) const {
  std::cout << "time to mat " << t << " spot " << x << " boundary right " << this->option->option_price_for_big_spot(t, x) << std::endl;
  return this->option->option_price_for_big_spot(t, x);
  // This is via Put-Call Parity and works for a call option
  //return (x-(option->K)*exp(-(option->r)*((option->T)-t))); 
}

// Initial condition 
double BlackScholesPDE::init_cond(double x) const {
  double res = this->option->option_price_at_maturity(x);
  //std::cout << "spot " << x << " payoff " << res << std::endl;
  return res;
  //return option->pay_off->operator()(x);
}

#endif