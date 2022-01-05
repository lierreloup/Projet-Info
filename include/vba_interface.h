#ifndef __VBA_INTERFACE_H
#define __VBA_INTERFACE_H

#include <string>
#include "../include/pricers.h"
#include "../include/discretization.h"
#include "../include/option.h"
#include "../include/pde.h"
#include "../include/fdm.h"
#include "../include/greeks.h"

/**
 * @brief structure which contains all inputs sent by the VBA app
 * 
 */
struct Input {
    std::string option_type;
    double spot, time_to_maturity, strike, rate, volatility, M, N; 
};
/**
 * @brief Reads the vba input (vba must write its input in a file) and assigns it to the proper parameters
 * 
 * @param filename input file ; the file should contain each parameter separated by a new line (or a space)
 */
Input get_params_from_file(std::string filename);

/**
 * @brief structure which contains all outputs needed by the VBA app
 * 
 */
struct Output {
    double price;
    double delta;
    double gamma;
    double theta;
    double rho;
    double vega;

    //etc
};

/**
 * @brief Fill an output with greeks so that VBA app gets it back
 * 
 * @param output 
 */
void store_greeks(Output &output,Input input);

/**
 * @brief Create an output file so that VBA app gets it back
 * 
 * @param filename 
 * @param output 
 */
void create_output_file(std::string filename, Output output);

double price_option(Input input);


#endif