#ifndef __FDM_H
#define __FDM_H

#include "../include/pde.h"
#include "discretization.h"
#include <vector>
#include <string> 

/**
 * @brief this is a utility struct in case we have to store some info about the PDE solving
 * (e.g spot values)
 * 
 */
struct NecessaryResults {
  std::vector<double> x_values;
};


// Finite Difference Method - Abstract Base Class
class FDMBase {
 protected:

  // Space discretisation
  Discretization const & disc;

  // Time discretisation (better refactoring would require that this is part of the previous attribute "disc")
  double t_dom;      // Temporal extent [0.0, t_dom]
  size_t N;   // Number of temporal differencing points

  // Constructor
  FDMBase(double _t_dom, size_t _N,
          Discretization const & _disc) 
          : t_dom(_t_dom), N(_N), disc(_disc) {}

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

};


class BSEuroImplicit : public FDMBase {
 protected:
  BlackScholesPDE const * pde;

// TODO : replace pointer to European option with reference
 public:
  BSEuroImplicit(double _t_dom, size_t _N,
                   EuropeanOption * european_option,
          Discretization const & _disc)
          : FDMBase(_t_dom, _N, _disc), pde(new BlackScholesPDE(european_option)) {}

  // return price
  std::vector<double> step_march(std::string output_file) override;
  ~BSEuroImplicit() { delete pde; }
};


class BSAmericanImplicit : public FDMBase {
  BlackScholesPDE const * no_early_exercise_pde;

 public:
  BSAmericanImplicit(double _t_dom, size_t _N,
                   AmericanOption * american_option,
          Discretization const & _disc)
        : FDMBase(_t_dom, _N, _disc), no_early_exercise_pde(new BlackScholesPDE(american_option)) {}
  
  std::vector<double> step_march(std::string output_file) override;

  ~BSAmericanImplicit() { delete no_early_exercise_pde; }

};
#endif