#ifndef __FDM_CPP
#define __FDM_CPP

#include <fstream>
#include <iostream>
#include <functional>
#include <math.h>
#include "../include/fdm.h"
#include "fdm_static.cpp"

FDMBase::FDMBase(double _x_dom, size_t _M,
                 double _t_dom, size_t _N)
  : x_dom(_x_dom), M(_M), t_dom(_t_dom), N(_N) {}


BSEuroImplicit::BSEuroImplicit(double _x_dom, size_t _M,
                                   double _t_dom, size_t _N,
                                   EuropeanOption * euro_option) 
  : FDMBase(_x_dom, _M, _t_dom, _N), pde(new BlackScholesPDE(euro_option)) {
  }



/*
  * Solve the PDE for all time steps and space steps
*/
std::vector<double> BSEuroImplicit::step_march(std::string output_file) {

  return step_march_uniform(
    output_file
    , this->x_dom
    , this->M
    , this->t_dom
    , this->N
    , this->pde
    , increment_european_price
    , this->results
  );
}
BSAmericanImplicitUniform::BSAmericanImplicitUniform(double _x_dom, size_t _M, double _t_dom, size_t _N, AmericanOption * american_option) 
: FDMBase(_x_dom, _M, _t_dom, _N), no_early_exercise_pde(new BlackScholesPDE(american_option)) {}

std::vector<double> BSAmericanImplicitUniform::step_march(std::string output_file) {
  return step_march_uniform(
    output_file
    , this->x_dom
    , this->M
    , this->t_dom
    , this->N
    , this->no_early_exercise_pde
    , increment_american_price
    , this->results
  );
}
BSAmericanImplicit::BSAmericanImplicit(double _x_dom, size_t _M, double _t_dom, size_t _N, AmericanOption * american_option) 
: FDMBase(_x_dom, _M, _t_dom, _N), no_early_exercise_pde(new BlackScholesPDE(american_option)) {}

std::vector<double> BSAmericanImplicit::step_march(std::string output_file) {
  return step_march_non_uniform(
    output_file
    , this->x_dom
    , this->M
    , this->t_dom
    , this->N
    , this->no_early_exercise_pde
    , increment_american_price
    , this->results
  );
}
#endif