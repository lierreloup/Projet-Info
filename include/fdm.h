#ifndef __FDM_H
#define __FDM_H

#include "../include/pde.h"
#include <vector>
#include <string> 

/**
 * @brief this is a utility struct in case we have to store some info about the PDE solving
 * (e.g spot values)
 * 
 */
struct NecessaryResults {
  std::vector<double> x_values;
  std::vector<double> t_values;
};


// Finite Difference Method - Abstract Base Class
class FDMBase {
 protected:

  // Space discretisation
  double x_dom;      // Spatial extent [0.0, x_dom]
  size_t M;   // Number of spatial differencing points

  // Time discretisation
  double t_dom;      // Temporal extent [0.0, t_dom]
  size_t N;   // Number of temporal differencing points

  // Constructor
  FDMBase(double _x_dom, size_t _M,
          double _t_dom, size_t _N);

 public:
  /**
    * @brief necessary for pricing : store some info about the solving
    * 
    */
  NecessaryResults results;

  /**
   * @brief solve the pde
   * 
   * @param output_file 
   * @return std::vector<double> the function values for final time and all space values
   * @
   */
  virtual std::vector<double> step_march(std::string output_file) = 0;

  // TODO : maybe make this a shared pointer for safety, or delete completely
  virtual ConvectionDiffusionPDE const * get_pde() = 0;
};


// TODO : rather than take a pde as argument, take a european option and create a BS for it
// TODO : remove most of methods

class BSEuroImplicit : public FDMBase {
 protected:
  BlackScholesPDE const * pde;

// TODO : replace pointer to European option with reference
 public:
  BSEuroImplicit(double _x_dom, size_t _M,
                   double _t_dom, size_t _N,
                   EuropeanOption * european_option);

  ConvectionDiffusionPDE const * get_pde() { return pde; }

  // return price
  std::vector<double> step_march(std::string output_file) override;
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

  std::vector<double> step_march(std::string output_file) override;

  ~BSAmericanImplicitUniform() { delete no_early_exercise_pde; }

};
#endif