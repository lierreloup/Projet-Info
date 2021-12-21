#ifndef __FDM_H
#define __FDM_H

#include "../include/pde.h"
#include <vector>
#include <string> 

// Finite Difference Method - Abstract Base Class
class FDMBase {
 protected:
  ConvectionDiffusionPDE* pde;

  // Space discretisation
  double x_dom;      // Spatial extent [0.0, x_dom]
  unsigned long M;   // Number of spatial differencing points
  double dx;         // Spatial step size (calculated from above)
  std::vector<double> x_values;  // Stores the coordinates of the x dimension

  // Time discretisation
  double t_dom;      // Temporal extent [0.0, t_dom]
  unsigned long N;   // Number of temporal differencing points
  double dt;         // Temporal step size (calculated from above)

  // Time-marching
  double prev_t, cur_t;   // Current and previous times

  // Differencing coefficients
  double alpha, beta, gamma;

  // Storage
  std::vector<double> new_result;   // New solution (becomes N+1)
  std::vector<double> old_result;   // Old solution (becomes N)

  // Constructor
  FDMBase(double _x_dom, unsigned long _M,
          double _t_dom, unsigned long _N,
          ConvectionDiffusionPDE* _pde);

  // Override these virtual methods in derived classes for 
  // specific FDM techniques, such as explicit Euler, Crank-Nicolson, etc.
  virtual void calculate_step_sizes() = 0;
  virtual void set_initial_conditions() = 0;
  virtual void calculate_boundary_conditions() = 0;
  virtual void calculate_inner_domain() = 0;

 public:
  // Carry out the actual time-stepping
  virtual void step_march(std::string output_file) = 0;
};

class FDMEulerExplicit : public FDMBase {
 protected:
  void calculate_step_sizes();
  void set_initial_conditions();
  void calculate_boundary_conditions();
  void calculate_inner_domain();

 public:
  FDMEulerExplicit(double _x_dom, unsigned long _M,
                   double _t_dom, unsigned long _N,
                   ConvectionDiffusionPDE* _pde);

  void step_march(std::string output_file);
};

class FDMEulerImplicit : public FDMBase {
 protected:
  void calculate_step_sizes();
  void set_initial_conditions();
  void calculate_boundary_conditions();
  void calculate_inner_domain();

 public:
  FDMEulerImplicit(double _x_dom, unsigned long _M,
                   double _t_dom, unsigned long _N,
                   ConvectionDiffusionPDE* _pde);

  void step_march(std::string output_file);
};

struct AmericanOptionParameters {
  BlackScholesPDE * no_early_exercise_pde;
};

class PriceAmericanOption : public FDMBase {
  void calculate_step_sizes();
  void set_initial_conditions();
  void calculate_boundary_conditions();
  void calculate_inner_domain();

 public:
  PriceAmericanOption(double _x_dom, unsigned long _M,
                   double _t_dom, unsigned long _N,
                   AmericanOptionParameters params);

  void step_march(std::string output_file);
};
#endif