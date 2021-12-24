#ifndef __VBA_INTERFACE_H
#define __VBA_INTERFACE_H

#include <string>

/**
 * @brief structure which contains all inputs sent by the VBA app
 * 
 */
struct Input {
    std::string option_type;
    double spot, time_to_maturity, strike, rate, volatility;
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
    //etc
};

/**
 * @brief Create an output file so that VBA app gets it back
 * 
 * @param filename 
 * @param output 
 */
void create_output_file(std::string filename, Output output);

#endif