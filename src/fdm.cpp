#ifndef __FDM_CPP
#define __FDM_CPP

#include <fstream>
#include <functional>
#include "../include/fdm.h"

FDMBase::FDMBase(double _x_dom, unsigned long _J,
                 double _t_dom, unsigned long _N,
                 ConvectionDiffusionPDE* _pde) 
  : x_dom(_x_dom), J(_J), t_dom(_t_dom), N(_N), pde(_pde) {}

FDMEulerExplicit::FDMEulerExplicit(double _x_dom, unsigned long _J,
                                   double _t_dom, unsigned long _N,
                                   ConvectionDiffusionPDE* _pde) 
  : FDMBase(_x_dom, _J, _t_dom, _N, _pde) {
  calculate_step_sizes();
  set_initial_conditions();
}

void FDMEulerExplicit::calculate_step_sizes() {
  dx = x_dom/static_cast<double>(J-1);
  dt = t_dom/static_cast<double>(N-1);
}

void FDMEulerExplicit::set_initial_conditions() {
  // Spatial settings
  double cur_spot = 0.0;

  old_result.resize(J, 0.0);
  new_result.resize(J, 0.0);
  x_values.resize(J, 0.0);

  for (unsigned long j=0; j<J; j++) {
    cur_spot = static_cast<double>(j)*dx;
    old_result[j] = pde->init_cond(cur_spot);
    x_values[j] = cur_spot;
  }

  // Temporal settings
  prev_t = 0.0;
  cur_t = 0.0;
}

void FDMEulerExplicit::calculate_boundary_conditions() {
  new_result[0] = pde->boundary_left(prev_t, x_values[0]);
  new_result[J-1] = pde->boundary_right(prev_t, x_values[J-1]);
}

void FDMEulerExplicit::calculate_inner_domain() {
  // Only use inner result indices (1 to J-2)
  for (unsigned long j=1; j<J-1; j++) {
    // Temporary variables used throughout
    double dt_sig = dt * (pde->diff_coeff(prev_t, x_values[j]));
    double dt_sig_2 = dt * dx * 0.5 * (pde->conv_coeff(prev_t, x_values[j]));

    // Differencing coefficients (see \alpha, \beta and \gamma in text)
    alpha = dt_sig - dt_sig_2;
    beta = dx * dx - (2.0 * dt_sig) + (dt * dx * dx * (pde->zero_coeff(prev_t, x_values[j])));
    gamma = dt_sig + dt_sig_2;

    // Update inner values of spatial discretisation grid (Explicit Euler)
    new_result[j] = ( (alpha * old_result[j-1]) + 
                      (beta * old_result[j]) + 
                      (gamma * old_result[j+1]) )/(dx*dx) - 
      (dt*(pde->source_coeff(prev_t, x_values[j])));
  }
}

void FDMEulerExplicit::step_march() { 
  std::ofstream fdm_out("fdm.csv");

  while(cur_t < t_dom) {
    cur_t = prev_t + dt;
    calculate_boundary_conditions();
    calculate_inner_domain();
    for (int j=0; j<J; j++) {
      fdm_out << x_values[j] << " " << prev_t << " " << new_result[j] << std::endl;
    }
    
    old_result = new_result;
    prev_t = cur_t;
  }

  fdm_out.close();
}


// La fête commence ici Castor





FDMEulerImplicit::FDMEulerImplicit(double _x_dom, unsigned long _J,
                                   double _t_dom, unsigned long _N,
                                   ConvectionDiffusionPDE* _pde) 
  : FDMBase(_x_dom, _J, _t_dom, _N, _pde) {
  calculate_step_sizes();
  set_initial_conditions();
}

void FDMEulerImplicit::calculate_step_sizes() {
  dx = x_dom/static_cast<double>(J-1);
  dt = t_dom/static_cast<double>(N-1);
}

// TODO : renommer J en N psk J c'est dégueu
void FDMEulerImplicit::set_initial_conditions() {
  // Spatial settings
  double cur_spot = 0.0;

  old_result.resize(J, 0.0);
  new_result.resize(J, 0.0);
  x_values.resize(J, 0.0);

  for (unsigned long j=0; j<J; j++) {
    cur_spot = static_cast<double>(j)*dx;
    old_result[j] = pde->init_cond(cur_spot);
    x_values[j] = cur_spot;
  }

  // Temporal settings
  prev_t = 0.0;
  cur_t = 0.0;
}

void FDMEulerImplicit::calculate_boundary_conditions() {
  // # balek absolu
}

void FDMEulerImplicit::calculate_inner_domain() {
  compute_next_line(
    x_values
    , dx
    , dt
    , prev_t
    , pde
    , old_result
    , new_result
  );
}

// Diagonal coefficients of the tridiagonal matrix A in implicit resolution
double compute_d(size_t n, double a, double b) {
  return 1 + b + a * n*n;
}

// Upper coefficients of the tridiagonal matrix A in implicit resolution
double compute_u(size_t n, double a, double b) {
  double n_min_1 = n-1;
  return -0.5 * (
    b * n_min_1 + a * n_min_1 * n_min_1
  );
}

// Lower coefficients of the tridiagonal matrix A in implicit resolution
double compute_l(size_t n, double a, double b) {
  double n_plus_1 = n+1;
  return 0.5 * (
    b * n_plus_1 - a * n_plus_1 * n_plus_1
  );
}
// TODO : maybe wrap the three functions in a Tridiagonal class for more clarity

// /!\ TODO : in the paper, indices of solution go from 1 to N-1...
std::vector<double> & Thomas_solver(const std::vector<double>& l, const std::vector<double>& d, const std::vector<double>& u, const std::vector<double> & b, std::vector<double>result, size_t begin_result, size_t end_result) {
  size_t n = end_result - begin_result;
  
  // cprime and dprime vectors are essential coefficients to Thomas' algorithm
  std::vector<double> cprime(n);
  std::vector<double> dprime(n);

  cprime[0] = u.at(1) / d.at(1);
  dprime[0] = b.at(0) / d.at(1);

  for (size_t i = 1; i < n - 1; i++) {
    // /!\ the indices for the current values of some variables start at 1 instead of 0
    double u_curr = u.at(i+1);
    double d_curr = d.at(i+1);
    double l_curr = l.at(i+1);
    double b_curr = b.at(i+1);

    cprime[i] = u_curr / (
      d_curr - l_curr * cprime[i-1]
    );

    dprime[i] = (
      b_curr - l_curr * dprime[i-1]
    ) / (
      d_curr - l_curr * cprime[i-1]
    );
  }

  // now solve the linear equation by finding x such that A x = b
  // TODO : we may replace at by [] for optimisation (at is safer because out of bounds is checked)
  result.at(end_result) = dprime[n-1];

  for (size_t i = end_result - 1; i >= begin_result; i--)  
  {
    double x_next = result.at(i+1);
    double x_curr = dprime[i] - cprime[i] * x_next;

    result.at(i) = x_curr;
  }

  return result;
}

std::vector<double> compute_b(std::vector<double> & old_result, double new_result_0, double new_result_N, double l0, double uN) {
  size_t n = old_result.size() - 2; // remove terms of index 0 and N

  std::vector<double> b(old_result); // calls the copy constructor

  b.erase(b.end());
  b.erase(b.begin());

  b.at(0) -= l0 * new_result_0;
  b.at(n-1) -= uN * new_result_N;

  return b;
}

std::vector<double> & calculate_boundary_conditions(double prev_t, std::vector<double> x_values, ConvectionDiffusionPDE * pde, std::vector<double> new_result) {
  size_t N = new_result.size() - 1;
  new_result[0] = pde->boundary_left(prev_t, x_values[0]);
  new_result[N-1] = pde->boundary_right(prev_t, x_values[N-1]);

  return new_result; // return the reference to the result array in case we need to affect it to a variable
}

// TODO : maybe rename "implicit Euler" as "BSImplicitEuler" because the scheme does'nt work for all convection PDEs
std::vector<double> & compute_next_line(std::vector<double> & x_values, double dx, double dt, double prev_t, ConvectionDiffusionPDE * pde, std::vector<double> & old_result, std::vector<double> & new_result) {
  size_t N = new_result.size() - 1; // results range from 0 to N

  calculate_boundary_conditions(prev_t, x_values, pde, new_result);

  // Define alpha, beta, d, u, l
  std::vector<double> alpha(N-2);
  std::vector<double> beta(N-2);
  std::vector<double> gamma(N-2);

  double a_coeff = 2.0 * pde->diff_coeff(prev_t, 1) * dt;
  double b_coeff = pde->conv_coeff(prev_t, 1) * dt;
  
  // Only use inner result indices (1 to J-2)
  for (size_t j=1; j<N-1; j++) {
    gamma.at(j-1) = compute_l(j, a_coeff, b_coeff);
    beta.at(j-1) = compute_d(j, a_coeff, b_coeff);
    alpha.at(j-1) = compute_u(j, a_coeff, b_coeff);
  }
    // I keep this just in case, TODO: delete later
    // // Temporary variables used throughout
    // double dt_sig = dt * (pde->diff_coeff(prev_t, j)); // = sigma² * j² * dt / 2 for Black Scholes
    // double dt_sig_2 = dt * 0.5 * (pde->conv_coeff(prev_t, j)); // = r * j * dt / 2 for Black Scholes

    // // Differencing coefficients (see \alpha, \beta and \gamma in text)
    // alpha.at(j) = dt_sig - dt_sig_2;
    // beta.at(j) = 1 - (2.0 * dt_sig) + (dt * (pde->zero_coeff(prev_t, j))); // = 1 - sigma² * j¹ * dt - r * dt for Black Scholes
    // gamma.at(j) = dt_sig + dt_sig_2;



  // Define vector b
  std::vector<double> b = compute_b(old_result, new_result[0], new_result[N], compute_l(0, a_coeff, b_coeff), compute_u(N, a_coeff, b_coeff));
 
  // Solve A v = b
  Thomas_solver(gamma, beta, alpha, b, new_result, 1, N-1);
  
  // Set result value equal to v
  return new_result; // return the reference to the result array in case we need to affect it to a variable

}

void FDMEulerImplicit::step_march() { 
  std::ofstream fdm_out("fdm.csv");

  while(cur_t < t_dom) {
    cur_t = prev_t + dt;
    calculate_boundary_conditions();
    calculate_inner_domain();
    for (int j=0; j<J; j++) {
      fdm_out << x_values[j] << " " << prev_t << " " << new_result[j] << std::endl;
    }
    
    old_result = new_result;
    prev_t = cur_t;
  }

  fdm_out.close();
}


#endif