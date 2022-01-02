#ifndef __FDM_STATIC_H
#define __FDM_STATIC_H

#include "../include/pde.h"
#include "../include/fdm.h"

// Diagonal coefficients of the tridiagonal matrix A in implicit resolution
double compute_d(size_t n, double a, double b) ;

// Upper coefficients of the tridiagonal matrix A in implicit resolution
double compute_u(size_t n, double a, double b) ;

// Lower coefficients of the tridiagonal matrix A in implicit resolution
double compute_l(size_t n, double a, double b) ;

std::vector<double> & Thomas_solver(const std::vector<double>& lower, const std::vector<double>& diagonal, const std::vector<double>& upper, const std::vector<double> & right_side, std::vector<double> & result, size_t begin_result, size_t end_result) ;
std::vector<double> compute_b(std::vector<double> & old_result, double new_result_0, double new_result_M, double l0, double uM) ;
std::vector<double> & store_boundary_conditions(double prev_t, std::vector<double> const & x_values, ConvectionDiffusionPDE const * pde, std::vector<double> & new_result) ;
std::vector<double> & increment_european_price(std::vector<double> const & x_values, double dt, double prev_t, ConvectionDiffusionPDE const * pde, std::vector<double> & old_result, std::vector<double> & new_result) ;
std::vector<double> get_uniform_x_grid(size_t M, double dx) ;
void BS_initial_conditions(size_t M, std::vector<double> & old_result, std::vector<double> & new_result, std::vector<double> const & x_values, ConvectionDiffusionPDE const * pde, double & prev_t, double & cur_t) ;
  typedef std::vector<double> & (* increment_function)(std::vector<double> const & x_values, double dt, double prev_t, ConvectionDiffusionPDE const * pde, std::vector<double> & old_result, std::vector<double> & new_result);

std::vector<double> step_march_aux(std::string output_file, BlackScholesPDE const * pde, size_t M, size_t N, std::vector<double> const & x_values, std::vector<double> const & time_to_maturity_values, increment_function increment_price) ;
 
/*
  * Solve the PDE for all time steps and space steps
*/
std::vector<double> step_march_uniform(std::string output_file, double x_dom, size_t M,
                   double t_dom, size_t N, BlackScholesPDE const * pde, increment_function increment_price, NecessaryResults & results) ;
 

// NON UNIFORM DISCRETIZATION

double compute_delta_eta(double M, double x_dom, double K, double c) ;
 double compute_eta(double i, double K, double c, double delta_eta) ;

std::vector<double> get_non_uniform_grid(size_t M, double x_dom, double K, double c) ;

std::vector<double> step_march_non_uniform(std::string output_file, double x_dom, size_t M,

                   double t_dom, size_t N, BlackScholesPDE const * pde, increment_function increment_price, NecessaryResults & results);

#endif