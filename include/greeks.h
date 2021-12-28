#ifndef __GREEK_H
#define __GREEK_H

#include <math.h>
#include <fstream>
#include <string>
#include <iostream>
#include <string.h>

double delta_option(double spot, double time_to_maturity, double strike, double rate, double volatility, std::string output_pde, const char* type_of_option);

double gamma_option(double spot, double time_to_maturity, double strike, double rate, double volatility, std::string output_pde, const char* type_of_option);

double thet_option(double spot, double time_to_maturity, double strike, double rate, double volatility, std::string output_pde, const char* type_of_option);

double rho_option(double spot, double time_to_maturity, double strike, double rate, double volatility, std::string output_pde, const char* type_of_option);

double vega_option(double spot, double time_to_maturity, double strike, double rate, double volatility, std::string output_pde, const char* type_of_option);


#endif
