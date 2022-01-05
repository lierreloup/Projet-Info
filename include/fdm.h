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

  // Space and time discretisation
  Discretization & disc;

  // Constructor
  FDMBase(Discretization & _disc) 
  : disc(_disc) {}

 public:
  /**
    * @brief necessary for pricing : store some info about the solving
    * 
    */
  NecessaryResults results;

  /**
   * @brief Solve the pde.
   * 
   * @param output_file A file where prices for the whole grid are stored
   * @return std::vector<double> the function values for final time and all space values (final line)
   */
  virtual std::vector<double> step_march(std::string output_file) = 0;

};


// IMPLEMENTATIONS FOR BLACK SCHOLES PRICING

class BSEuroImplicit : public FDMBase {
 protected:
  BlackScholesPDE const * pde;

 public:
  BSEuroImplicit(EuropeanOption * european_option, UniformDiscretization & _disc)
  : FDMBase(_disc), pde(new BlackScholesPDE(european_option)) {}

  // return price
  std::vector<double> step_march(std::string output_file) override;

  ~BSEuroImplicit() { delete pde; }
};





class BSAmericanImplicit : public FDMBase {
  BlackScholesPDE const * no_early_exercise_pde; // PDE in the no exercise region (see the Horng-Tien reference)

 public:
  BSAmericanImplicit(AmericanOption * american_option, UniformDiscretization & _disc)
  : FDMBase(_disc), no_early_exercise_pde(new BlackScholesPDE(american_option)) {}
  
  std::vector<double> step_march(std::string output_file) override;

  ~BSAmericanImplicit() { delete no_early_exercise_pde; }

};
#endif