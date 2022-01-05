#ifndef __DISCRETIZATION_CPP
#define __DISCRETIZATION_CPP

#include "../include/discretization.h"
#include <math.h>
#include <iostream>
#include <stddef.h>

// UNIFORM DISCRETIZATION

Discretization::~Discretization(){}

std::vector<double> get_uniform_x_grid(size_t M, double dx) {
  std::vector<double> x_values(M+1);

  for (size_t m=0; m <= M; m++) {
    double cur_spot = static_cast<double>(m) * dx;
    x_values.at(m) = cur_spot;
  }

  return x_values;
}

std::vector<double> UniformDiscretization::get_x_grid() const {
    return get_uniform_x_grid(
      this->M
      , this->x_dom / this->M
    );
}

std::vector<double> UniformDiscretization::get_t_grid() const {
    return get_uniform_x_grid(
      this->N
      , this->t_dom / this->N
    );
}

#endif