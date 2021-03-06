#ifndef __FDM_CPP_STATIC
#define __FDM_CPP_STATIC

#include <fstream>
#include <iostream>
#include <functional>
#include <math.h>
#include "../include/pde.h"
#include "../include/fdm.h"



// IMPLICIT METHOD : UTILITY FUNCTIONS





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

// Right-hand side of the system Ax = b in implicit resolution
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

/**
 * @brief Solves a tridiagonal system A x = b.
 * 
 * @param lower part of A
 * @param diagonal part of A
 * @param upper part of A
 * @param right_side b
 * @param result x is stored in a result_vector, from index 'begin_result' to index 'end_result'
 * @param begin_result 
 * @param end_result 
 * @attention the reason begin_result and end_result are present is that, when using Thomas solver in implicit FDM,
 * the result vector can be written as : new_result[0->M] but we only need to compute : new_result[1->M-1] because
 * indices 0 and M are known as boundary conditions
 * @ref https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
 * @return std::vector<double>& a reference to the result vector (in case we want to store its reference elsewhere)
 */
std::vector<double> & Thomas_solver(const std::vector<double>& lower, const std::vector<double>& diagonal, const std::vector<double>& upper, const std::vector<double> & right_side, std::vector<double> & result, size_t begin_result, size_t end_result) {
  try
  { 
  size_t diagonal_size = end_result - begin_result + 1;

  // cprime and dprime vectors are essential coefficients to Thomas' algorithm
  std::vector<double> cprime(diagonal_size);
  std::vector<double> dprime(diagonal_size);

  // initial index
  cprime.at(0) = upper.at(0) / diagonal.at(0);
  dprime.at(0) = right_side.at(0) / diagonal.at(0);

  for (size_t i = 1; i < diagonal_size - 1; i++) {
    // to avoid getting a headache, coefficients are named like in the wikipedia page https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    // even though Wikipedia indices range from 1 to n, ours range from 0 to diagonal_size - 1
    double ci = upper.at(i);
    double bi = diagonal.at(i);
    double ai = lower.at(i-1);
    double di = right_side.at(i);

    // wrote fractions on two lines so that they can be read as fractions
    cprime.at(i) = ci / (
      bi - ai * cprime.at(i-1)
    );

    dprime.at(i) = (
      di - ai * dprime.at(i-1)
    ) / (
      bi - ai * cprime.at(i-1)
    );
  }
  // Last index is different, cprime is not defined here
  double bi = diagonal.at(diagonal_size - 1);
  double ai = lower.at(diagonal_size - 2);
  double di = right_side.at(diagonal_size - 1);

  dprime.at(diagonal_size - 1) = (
    di - ai * dprime.at(diagonal_size - 2)
  ) / (
    bi - ai * cprime.at(diagonal_size - 2)
  );

  // now use cprime and dprime to solve the linear equation by finding x such that A x = b
  result.at(end_result) = dprime.at(diagonal_size-1);

  size_t i = diagonal_size - 2;

  // The condition for the loop is weird because "size_t" handles "result_index--" pretty bad (in case it should go under 0)
  for (size_t result_index = end_result - 1; result_index >= begin_result && result_index <= end_result; result_index--)  
  {
    double x_next = result.at(result_index+1);

    double x_curr = dprime.at(i) - cprime.at(i) * x_next;

    result.at(result_index) = x_curr;

    i--;
  }

  return result;

  }
  catch(const std::exception& e)
  {
    std::cerr << "Thomas solver " << e.what() << '\n';
    throw e;
  }
}






// INCREMENT FUNCTIONS

/**
 * These are the most technical functions, which deduce the next price line from the previous one.
 * A price line consists of prices for all spots at one point in time.
**/






// Deduce prices at next time step from prices at previous step and linear system solving (European case)
std::vector<double> & increment_european_price(std::vector<double> const & x_values, double dt, double prev_t, ConvectionDiffusionPDE const * pde, std::vector<double> & old_result, std::vector<double> & new_result) {
  size_t M = new_result.size() - 1; // results range from 0 to N

  // Define alpha, beta, d, u, l
  // in the paper, M is called N
  std::vector<double> upper(M-2); //in the paper, indices range from 2 to M
  std::vector<double> diagonal(M-1); //in the paper, indices range from 1 to M-1
  std::vector<double> lower(M-2); //in the paper, indices range from 0 to M-2

  double a_coeff = 2.0 * pde->diff_coeff(prev_t, 1) * dt; // this equals sigma?? * dt
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


// This is an american pricing implementation ispired by the Horng Tien method (see references)
// Computes a price line exactly like for european option, then check if the time derivative is strictly smaller than 0, 
// in which case it is set to 0 because american prices can only increase backwards in time
std::vector<double> & increment_american_price(std::vector<double> const & x_values, double dt, double prev_t, ConvectionDiffusionPDE const * pde, std::vector<double> & old_result, std::vector<double> & new_result) {
  
  increment_european_price(
    x_values
    , dt
    , prev_t
    , pde
    , old_result
    , new_result
  );

  size_t M = new_result.size() - 1;

  // We use the Horng-Tien method to compute American Option prices
  // That is, if the backward time derivative is strictly smaller than 0, set it to 0
  for (int m = 0; m <= M; m++) {
    try
    {
      double time_derivative_times_dt = new_result.at(m) - old_result.at(m);
      if (time_derivative_times_dt < 0) {
        new_result[m] = old_result[m]; 
        //experience seems to show that this never happens when american and european exercise are
        //equivalent
      }
    }
    catch(const std::exception& e)
    {
      std::cerr << "increment_american_price : " << e.what() << '\n';
    }
    
  }

  return new_result;
}


// This typedef of a function pointer is such that we can pass american and european price increment functions as parameters
// This is useful in the step_march function
typedef std::vector<double> & (* increment_function)(std::vector<double> const & x_values, double dt, double prev_t, ConvectionDiffusionPDE const * pde, std::vector<double> & old_result, std::vector<double> & new_result);








// ACTUAL PDE SOLVING








// Values for the initial time step (0)
void BS_initial_conditions(size_t M, std::vector<double> & old_result, std::vector<double> & new_result, std::vector<double> const & x_values, ConvectionDiffusionPDE const * pde, double & prev_t, double & cur_t) {

  for (size_t m=0; m <= M; m++) {
    double cur_spot = x_values.at(m);
    old_result.at(m) = pde->init_cond(cur_spot);
  
  }
  prev_t = 0;
  cur_t = 0;
}



// At each time step, compute and store boundary conditions for the PDE
std::vector<double> & store_boundary_conditions(double prev_t, std::vector<double> const & x_values, ConvectionDiffusionPDE const * pde, std::vector<double> & new_result) {
  size_t M = new_result.size() - 1;
  new_result.at(0) = pde->boundary_left(prev_t, x_values[0]);
  new_result.at(M) = pde->boundary_right(prev_t, x_values[M]);

  return new_result; // return the reference to the result array in case we need to affect it to a variable
}




// Auxilliary function to step_march, which carries out most of the work
std::vector<double> step_march_aux(std::string output_file, BlackScholesPDE const * pde, size_t M, size_t N, std::vector<double> const & x_values, std::vector<double> const & time_to_maturity_values, increment_function increment_price) {
  
  // Initialize everything
  double prev_t = 0;
  double cur_t = 0;
  std::vector<double> old_result(M+1, 0), new_result(M+1, 0);
  
  BS_initial_conditions(
    M
    , old_result
    , new_result
    , x_values
    , pde
    , prev_t
    , cur_t
  );


  // Find prices for the whole spacetime grid
  try
  {
  // File to store the prices of all lines
  std::ofstream fdm_out(output_file);

  // Store the first time step
  for (int m=0; m <= M; m++) {
    fdm_out << x_values.at(m) << " " << cur_t << " " << old_result.at(m) << std::endl;
  }
  
  // All other time steps :
  for (size_t n = 1; n <= N; n++) {
    cur_t = time_to_maturity_values[n];
    double dt = cur_t - prev_t;

    store_boundary_conditions(
      cur_t
      , x_values
      , pde
      , new_result
    );

    increment_price(
      x_values
      , dt
      , prev_t
      , pde
      , old_result
      , new_result
    );

    for (int m=0; m <= M; m++) {
      //store x, t and price in a new line of the file
      fdm_out << x_values.at(m) << " " << cur_t << " " << new_result.at(m) << std::endl;
    }
    
    old_result = new_result;
    prev_t = cur_t;
  }

  return new_result;

  fdm_out.close();
  }

  catch(const std::exception& e)
  {
    std::cerr << "Euler Implicit step_march : " << e.what() << '\n';
    throw e;
  }
 
}



/**
 * @brief This function is meant to solve a pde (as in the step march method from FDM classes)
 * 
 * @param output_file where to store prices
 * @param pde 
 * @param increment_price the increment function (e.g american or european)
 * @param results in case we want to store some of the computations in RAM
 * @param disc the discretization grid
 * @return std::vector<double> the last price line
 */
std::vector<double> step_march_function(std::string output_file, BlackScholesPDE const * pde, increment_function increment_price, NecessaryResults & results, Discretization const & disc) {

  std::vector<double> time_to_maturity_values = disc.get_t_grid();

  std::vector<double> x_values = disc.get_x_grid();
  results.x_values = std::vector<double>(x_values);

  size_t M = disc.get_M()
  , N = disc.get_N()
  ;

  return step_march_aux(
    output_file
    , pde
    , M
    , N
    , x_values
    , time_to_maturity_values
    , increment_price
  );
}

#endif