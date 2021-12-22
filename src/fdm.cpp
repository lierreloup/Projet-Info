#ifndef __FDM_CPP
#define __FDM_CPP

#include <fstream>
#include <iostream>
#include <functional>
#include "../include/fdm.h"

FDMBase::FDMBase(double _x_dom, unsigned long _M,
                 double _t_dom, unsigned long _N,
                 ConvectionDiffusionPDE* _pde) 
  : x_dom(_x_dom), M(_M), t_dom(_t_dom), N(_N), pde(_pde) {}

FDMEulerExplicit::FDMEulerExplicit(double _x_dom, unsigned long _M,
                                   double _t_dom, unsigned long _N,
                                   ConvectionDiffusionPDE* _pde) 
  : FDMBase(_x_dom, _M, _t_dom, _N, _pde) {
  calculate_step_sizes();
  set_initial_conditions();
}

void FDMEulerExplicit::calculate_step_sizes() {
  dx = x_dom/static_cast<double>(M-1);
  dt = t_dom/static_cast<double>(N-1);
}

void FDMEulerExplicit::set_initial_conditions() {
  // Spatial settings
  double cur_spot = 0.0;

  old_result.resize(M, 0.0);
  new_result.resize(M, 0.0);
  x_values.resize(M, 0.0);

  for (unsigned long m=0; m<M; m++) {
    cur_spot = static_cast<double>(m)*dx;
    old_result[m] = pde->init_cond(cur_spot);
    x_values[m] = cur_spot;
  }

  // Temporal settings
  prev_t = 0.0;
  cur_t = 0.0;
}

void FDMEulerExplicit::calculate_boundary_conditions() {
  new_result[0] = pde->boundary_left(prev_t, x_values[0]);
  new_result[M-1] = pde->boundary_right(prev_t, x_values[M-1]);
}

void FDMEulerExplicit::calculate_inner_domain() {
  // Only use inner result indices (1 to M-2)
  for (unsigned long m=1; m<M-1; m++) {
    // Temporary variables used throughout
    double dt_sig = dt * (pde->diff_coeff(prev_t, x_values[m]));
    double dt_sig_2 = dt * dx * 0.5 * (pde->conv_coeff(prev_t, x_values[m]));

    // Differencing coefficients (see \alpha, \beta and \gamma in text)
    alpha = dt_sig - dt_sig_2;
    beta = dx * dx - (2.0 * dt_sig) + (dt * dx * dx * (pde->zero_coeff(prev_t, x_values[m])));
    gamma = dt_sig + dt_sig_2;

    // Update inner values of spatial discretisation grid (Explicit Euler)
    new_result[m] = ( (alpha * old_result[m-1]) + 
                      (beta * old_result[m]) + 
                      (gamma * old_result[m+1]) )/(dx*dx) - 
      (dt*(pde->source_coeff(prev_t, x_values[m])));
  }
}

void FDMEulerExplicit::step_march(std::string output_file) { 
  std::ofstream fdm_out("fdm.csv");

  while(cur_t < t_dom) {
    cur_t = prev_t + dt;
    calculate_boundary_conditions();
    calculate_inner_domain();
    for (int m=0; m<M; m++) {
      fdm_out << x_values[m] << " " << prev_t << " " << new_result[m] << std::endl;
    }
    
    old_result = new_result;
    prev_t = cur_t;
  }

  fdm_out.close();
}


// La fête commence ici Castor





FDMEulerImplicit::FDMEulerImplicit(double _x_dom, unsigned long _M,
                                   double _t_dom, unsigned long _N,
                                   ConvectionDiffusionPDE* _pde) 
  : FDMBase(_x_dom, _M, _t_dom, _N, _pde) {
  calculate_step_sizes();
  set_initial_conditions();
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
std::vector<double> & Thomas_solver(const std::vector<double>& lower, const std::vector<double>& diagonal, const std::vector<double>& upper, const std::vector<double> & right_side, std::vector<double> & result, size_t begin_result, size_t end_result) {
  try
  {  
  size_t diagonal_size = end_result - begin_result + 1;

  // cprime and dprime vectors are essential coefficients to Thomas' algorithm
  std::vector<double> cprime(diagonal_size);
  std::vector<double> dprime(diagonal_size);

  cprime.at(0) = upper.at(0) / diagonal.at(0);
  dprime.at(0) = right_side.at(0) / diagonal.at(0);

  for (size_t i = 1; i < diagonal_size - 1; i++) {
    std::cout << i << " diagonal " << diagonal_size << std::endl;
    // to avoid getting a headache, coefficients are named like in the wikipedia page https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    // even though Wikipedia indices range from 1 to n, ours range from 0 to diagonal_size - 1
    double ci = upper.at(i);
    double bi = diagonal.at(i);
    double ai = lower.at(i-1);
    double di = right_side.at(i);

    // I wrote fractions on two lines so that they can be read as fractions
    cprime.at(i) = ci / (
      bi - ai * cprime.at(i-1)
    );

    dprime.at(i) = (
      di - ai * dprime.at(i-1)
    ) / (
      bi - ai * cprime.at(i-1)
    );
  }

  double bi = diagonal.at(diagonal_size - 1);
  double ai = lower.at(diagonal_size-2);
  double di = right_side.at(diagonal_size - 1);

  dprime.at(diagonal_size - 1) = (
    di - ai * dprime.at(diagonal_size - 2)
  ) / (
    bi - ai * cprime.at(diagonal_size - 2)
  );

  std::cout << "after loop";
  // now solve the linear equation by finding x such that A x = b
  result.at(end_result) = dprime.at(diagonal_size-1);
  std::cout << "after end result";

  size_t i = diagonal_size - 2;
  for (size_t result_index = end_result - 1; result_index >= begin_result; result_index--)  
  {

    double x_next = result.at(result_index+1);
    double x_curr = dprime.at(i) - cprime.at(i) * x_next;

    result.at(result_index) = x_curr;
    i --;
  }

  return result;

  }
  catch(const std::exception& e)
  {
    std::cerr << "Thomas solver " << e.what() << '\n';
    throw e;
  }
}

std::vector<double> compute_b(std::vector<double> & old_result, double new_result_0, double new_result_M, double l0, double uM) {
  try
  {
  size_t M = old_result.size() - 1;
  size_t size_b = M - 1; // b indices in the paper range from 1 to M-1

  std::vector<double> b(size_b);

  for (size_t m = 1; m <= size_b; m++) {
    b.at(m-1) = old_result.at(m);
  }

  b.at(0) -= l0 * new_result_0;
  b.at(size_b-1) -= uM * new_result_M;

  return b;

  }

  catch(const std::exception& e)
  {
    std::cerr << "compute_b : " << e.what() << '\n';
    throw e;
  }
}

std::vector<double> & store_boundary_conditions(double prev_t, std::vector<double> & x_values, ConvectionDiffusionPDE * pde, std::vector<double> & new_result) {
  size_t M = new_result.size() - 1;
  new_result.at(0) = pde->boundary_left(prev_t, x_values[0]);
  new_result.at(M) = pde->boundary_right(prev_t, x_values[M]);

  return new_result; // return the reference to the result array in case we need to affect it to a variable
}

// TODO : maybe rename "implicit Euler" as "BSImplicitEuler" because the scheme does'nt work for all convection PDEs
std::vector<double> & compute_next_line(std::vector<double> & x_values, double dx, double dt, double prev_t, ConvectionDiffusionPDE * pde, std::vector<double> & old_result, std::vector<double> & new_result) {
  size_t M = new_result.size() - 1; // results range from 0 to N

  // store_boundary_conditions(prev_t, x_values, pde, new_result);

  // Define alpha, beta, d, u, l
  // in the paper, M is called N
  std::vector<double> upper(M-2); //in the paper, indices range from 2 to M
  std::vector<double> diagonal(M-1); //in the paper, indices range from 1 to M-1
  std::vector<double> lower(M-2); //in the paper, indices range from 0 to M-2

  double a_coeff = 2.0 * pde->diff_coeff(prev_t, 1) * dt; // this equals sigma² * dt
  double b_coeff = pde->conv_coeff(prev_t, 1) * dt; // this equals r * dt
  
  // Only use inner result indices (1 to M-1)
  for (size_t m = 1; m <= M-2; m++) {
    lower.at(m-1) = compute_l(m, a_coeff, b_coeff); // l indices go from 0 to M-2
  }

  for (size_t m=1; m <= M-1; m++) {
    diagonal.at(m-1) = compute_d(m, a_coeff, b_coeff); // d indices go from 1 to M-1
  }

  for (size_t m=2; m <= M-1; m++) {
    upper.at(m-2) = compute_u(m, a_coeff, b_coeff); // u indices go from 2 to M-1
  }


  // Define vector b
  std::vector<double> b = compute_b(old_result, new_result.at(0), new_result.at(M), compute_l(0, a_coeff, b_coeff), compute_u(M, a_coeff, b_coeff));
 
  // Solve A v = b
  Thomas_solver(lower, diagonal, upper, b, new_result, 1, M-1);
  
  // Set result value equal to v
  return new_result; // return the reference to the result array in case we need to affect it to a variable

}

void FDMEulerImplicit::calculate_step_sizes() {
  dx = x_dom/static_cast<double>(M);
  dt = t_dom/static_cast<double>(N);
}

// TO MAKE THIS GENERIC : use pointers/const refs/setters
void FDMEulerImplicit::set_initial_conditions() {
  // Spatial settings
  double cur_spot = 0.0;

  this->old_result.resize(M+1, 0.0);
  this->new_result.resize(M+1, 0.0);
  this->x_values.resize(M+1, 0.0);

  for (size_t m=0; m <= M; m++) {
    cur_spot = static_cast<double>(m)*dx;
    std::cout << "m " << m << " cur_spot " << cur_spot << " dx " << dx << std::endl;
    this->old_result.at(m) = pde->init_cond(cur_spot);
    this->x_values.at(m) = cur_spot;
  }

  // Temporal settings TODO : IS THIS TIME TO MAT OR TIME
  this->prev_t = 0;
  this->cur_t = 0;
}

void FDMEulerImplicit::calculate_boundary_conditions() {
  store_boundary_conditions(
    this->prev_t
    , this->x_values
    , this->pde
    , this->new_result
  );
}

void FDMEulerImplicit::calculate_inner_domain() {
  compute_next_line(
    this->x_values
    , this->dx
    , this->dt
    , this->prev_t
    , this->pde
    , this->old_result
    , this->new_result
  );
}


/*
  * Solve the PDE for all time steps and space steps
*/
void FDMEulerImplicit::step_march(std::string output_file) {
  try
  {
 
  std::ofstream fdm_out(output_file);

  for (int m=0; m <= this->M; m++) {
    fdm_out << this->x_values.at(m) << " " << this->cur_t << " " << this->old_result.at(m) << std::endl;
  }
  
  double upper_time = this->t_dom - dt;
  while(this->prev_t <= upper_time) {
    this->cur_t = this->prev_t + this->dt;
    this->calculate_boundary_conditions();
    this->calculate_inner_domain();
    for (int m=0; m <= this->M; m++) {
      //store x, t and price in a new line of the file
      // TODO : why prev_t ?
      fdm_out << this->x_values.at(m) << " " << this->cur_t << " " << this->new_result.at(m) << std::endl;
    }
    
    this->old_result = this->new_result;
    this->prev_t = this->cur_t;
  }

  fdm_out.close();
  }

  catch(const std::exception& e)
  {
    std::cerr << "Euler Implicit step_march : " << e.what() << '\n';
    throw e;
  }
 
}

  void PriceAmericanOption::calculate_step_sizes() {
    dx = x_dom/static_cast<double>(M-1);
    dt = t_dom/static_cast<double>(N-1);
  }

  void PriceAmericanOption::set_initial_conditions() {
    // Spatial settings
    double cur_spot = 0.0;

    this->old_result.resize(M, 0.0);
    this->new_result.resize(M, 0.0);
    this->x_values.resize(M, 0.0);

    for (unsigned long m=0; m<M; m++) {
      cur_spot = static_cast<double>(m)*dx;
      this->old_result.at(m) = pde->init_cond(cur_spot);
      this->x_values.at(m) = cur_spot;
    }

    // Temporal settings
    this->prev_t = 0.0;
    this->cur_t = 0.0;
  }

  //TODO : check
  void PriceAmericanOption::calculate_boundary_conditions() {
    store_boundary_conditions(
      this->prev_t
      , this->x_values
      , this->pde
      , this->new_result
    );
  }

  //calculate line exactly like european option, then check if Vm+1 - Vm strictly greater than 0 (in which case set to 0)
  std::vector<double> & increment_american_price(std::vector<double> & x_values, double dx, double dt, double prev_t, ConvectionDiffusionPDE * pde, std::vector<double> & old_result, std::vector<double> & new_result) {
    compute_next_line(
      x_values
      , dx
      , dt
      , prev_t
      , pde
      , old_result
      , new_result
    );

    size_t N = new_result.size();

    // We use the Horng-Tien method to compute American Option prices
    // That is, if the time derivative is strictly greater than 0, set it to 0
    for (int n = 0; n < N; n++) {
      try
      {
        double time_derivative_times_dt = new_result.at(n) - old_result.at(n);
        if (time_derivative_times_dt > 0) {
          new_result[n] = 0;
        }
      }
      catch(const std::exception& e)
      {
        std::cerr << "increment_american_price : " << e.what() << '\n';
      }
      
    }

    return new_result;
  }

  void PriceAmericanOption::calculate_inner_domain() {
    increment_american_price(
      this->x_values
      , this->dx
      , this->dt
      , this->prev_t
      , this->pde
      , this->old_result
      , this->new_result
    );
  }

  PriceAmericanOption::PriceAmericanOption(double _x_dom, unsigned long _M, double _t_dom, unsigned long _N, AmericanOptionParameters params) 
  : FDMBase(_x_dom, _M, _t_dom, _N, params.no_early_exercise_pde)
  {
    calculate_step_sizes();
    set_initial_conditions();
  }

  void PriceAmericanOption::step_march(std::string output_file) {
    std::ofstream fdm_out(output_file);

    double upper_time = this->t_dom - dt;
    while(this->prev_t <= upper_time) {
      this->cur_t = this->prev_t + this->dt;
      this->calculate_boundary_conditions();
      this->calculate_inner_domain();
      for (int m=0; m<this->M; m++) {
        //store x, t and price in a new line of the file
        fdm_out << this->x_values.at(m) << " " << this->prev_t << " " << this->new_result.at(m) << std::endl;
      }
      
      this->old_result = this->new_result;
      this->prev_t = this->cur_t;
    }

    fdm_out.close();
  }


#endif