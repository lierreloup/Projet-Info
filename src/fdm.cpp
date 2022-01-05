#ifndef __FDM_CPP
#define __FDM_CPP

#include <fstream>
#include <iostream>
#include <functional>
#include <math.h>
#include "../include/fdm.h"
#include "fdm_static.cpp" // This is where all important functions are

/*
  * Solve the PDE for all time steps and space steps
*/
std::vector<double> BSEuroImplicit::step_march(std::string output_file) {

  return step_march_function(
    output_file
    , this->pde
    , increment_european_price
    , this->results
    , this->disc
  );
}

std::vector<double> BSAmericanImplicit::step_march(std::string output_file) {
  return step_march_function(
    output_file
    , this->no_early_exercise_pde
    , increment_american_price
    , this->results
    , this->disc
  );
}
#endif