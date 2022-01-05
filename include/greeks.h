#ifndef __GREEK_H
#define __GREEK_H

#include <math.h>
#include <fstream>
#include <string>
#include <iostream>
#include <string.h>
#include "../include/pricers.h"
#include "../include/pde.h"
#include "../include/discretization.h"
#include "../include/vba_interface.h"

double delta_option(price_inputs in, std::string output_pde, const char* type_of_option);

double gamma_option(price_inputs in, std::string output_pde, const char* type_of_option);

double theta_option(price_inputs in, std::string output_pde, const char* type_of_option);

double rho_option(price_inputs in, std::string output_pde, const char* type_of_option);

double vega_option(price_inputs in, std::string output_pde, const char* type_of_option);


#endif
