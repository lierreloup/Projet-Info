#ifndef __FDM_H
#define __FDM_H

#include "../include/pde.h"
#include <vector>
#include <string> 

// Finite Difference Method - Abstract Base Class
class FDMBase {
 protected:

  // Space discretisation
  double x_dom;      // Spatial extent [0.0, x_dom]
  size_t M;   // Number of spatial differencing points
  double dx;         // Spatial step size (calculated from above)

  // Time discretisation
  double t_dom;      // Temporal extent [0.0, t_dom]
  size_t N;   // Number of temporal differencing points

  // Constructor
  FDMBase(double _x_dom, size_t _M,
          double _t_dom, size_t _N);

 public:
  // Carry out the actual time-stepping
  virtual void step_march(std::string output_file) = 0;
  virtual ConvectionDiffusionPDE const * get_pde() = 0;
};


// TODO : rather than take a pde as argument, take a european option and create a BS for it
// TODO : remove most of methods

class BSEuroImplicit : public FDMBase {
 protected:
  BlackScholesPDE const * pde; // TODO : does this work with base constructor ?

 public:
  BSEuroImplicit(double _x_dom, size_t _M,
                   double _t_dom, size_t _N,
                   EuropeanOption * european_option);

  ConvectionDiffusionPDE const * get_pde() { return pde; }

  void step_march(std::string output_file);
  ~BSEuroImplicit() { delete pde; }
};

/*
struct AmericanOptionParameters {
  BlackScholesPDE const * no_early_exercise_pde;

};
*/

class BSAmericanImplicitUniform : public FDMBase {
  BlackScholesPDE const * no_early_exercise_pde;

 public:
  BSAmericanImplicitUniform(double _x_dom, size_t _M,
                   double _t_dom, size_t _N,
                   AmericanOption * american_option);
  
  ConvectionDiffusionPDE const * get_pde() { return no_early_exercise_pde; }
  void step_march(std::string output_file);
};
#endif