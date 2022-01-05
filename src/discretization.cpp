#ifndef __DISCRETIZATION_CPP
#define __DISCRETIZATION_CPP

#include "../include/discretization.h"
#include <math.h>
#include <iostream>

// UNIFORM DISCRETIZATION

std::vector<double> get_uniform_x_grid(size_t M, double dx) {
  std::vector<double> x_values(M+1);

  for (size_t m=0; m <= M; m++) {
    double cur_spot = static_cast<double>(m) * dx;
    x_values.at(m) = cur_spot;
  }

  return x_values;
}

std::vector<double> UniformDiscretization::get_grid() const {
    return get_uniform_x_grid(
      this->M
      , this->x_dom / this->M
    );
}

// NON UNIFORM DISCRETIZATION

double compute_delta_eta(double M, double x_dom, double K, double c) {
  return (
    asinh((x_dom - K) / c) - asinh(-K / c)
  ) / (
    M
  );
}

double compute_eta(double i, double K, double c, double delta_eta) {
  return asinh(-K / c) + i * delta_eta;
}

std::vector<double> get_non_uniform_grid(size_t M, double x_dom, double K, double c) {
  //std::cout << "M " << M << " x_dom " << x_dom << " K " << K << " c " << c << '\n';
  std::vector<double> x_values(M+1);
  double delta_eta = compute_delta_eta(static_cast<double>(M), x_dom, K, c);

  for (size_t i = 0; i <= M; i++) {
    double eta_i = compute_eta(static_cast<double> (i), K, c, delta_eta);
    x_values.at(i) = K + c * sinh(eta_i);
  }

  return x_values;
}

std::vector<double> NonUniformDiscretization::get_grid() const {
  return get_non_uniform_grid(
    this->M
    , this->x_dom
    , this->disc_center
    , this->c
  );
}

#endif