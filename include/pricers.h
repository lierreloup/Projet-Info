#ifndef __PRICERS_H
#define __PRICERS_H

#include <string>
#include "discretization.h"

struct price_inputs {
    double spot
        , time_to_maturity
        , strike
        , rate
        , volatility
    ;

    Discretization * disc; // here we use a pointer so we don't have to make a constructor for price_inputs (wouldn't work for a reference)
    size_t N;
};

/**
 * @brief get price for european call
 * 
 * @param spot 
 * @param time_to_maturity 
 * @param strike 
 * @param rate 
 * @param volatility 
 * @param output_pde 
 * @return double 
 */
double price_european_call(price_inputs in, std::string output_pde = "european_call.csv");

/**
 * @brief get price for european put
 * 
 * @param spot 
 * @param time_to_maturity 
 * @param strike 
 * @param rate 
 * @param volatility 
 * @param output_pde 
 * @return double 
 */
double price_european_put(price_inputs in, std::string output_pde = "european_put.csv");

/**
 * @brief get price for american call
 * 
 * @param spot 
 * @param time_to_maturity 
 * @param strike 
 * @param rate 
 * @param volatility 
 * @param output_pde 
 * @return double 
 */
double price_american_call(price_inputs in, std::string output_pde = "american_call.csv");

/**
 * @brief get price for american put
 * 
 * @param spot 
 * @param time_to_maturity 
 * @param strike 
 * @param rate 
 * @param volatility 
 * @param output_pde 
 * @return double 
 */
double price_american_put(price_inputs in, std::string output_pde = "american_put.csv");

UniformDiscretization default_UniformDiscretization(price_inputs in);

NonUniformDiscretization default_NonUniformDiscretization(price_inputs in);

size_t default_N(price_inputs in);

// "BONUS" functions

// Generates put prices from call prices by put call parity ;
// gives pretty good results, though values for big spots may get negative
// (but very close to 0)
// Allowed to make a quick "check" of our pde put prices
void getPutPricesFromCallFile(std::string call_prices_file, double rate, double strike, std::string output_file);

#endif