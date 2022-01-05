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

    UniformDiscretization * disc; // here we use a pointer so we don't have to make a constructor for price_inputs (wouldn't work for a reference)
};

/**
 * @brief get price for european call
 */
double price_european_call(price_inputs in, std::string output_pde = "european_call.csv");

/**
 * @brief get price for european put
 */
double price_european_put(price_inputs in, std::string output_pde = "european_put.csv");

/**
 * @brief get price for american call
 */
double price_american_call(price_inputs in, std::string output_pde = "american_call.csv");

/**
 * @brief get price for american put
 */
double price_american_put(price_inputs in, std::string output_pde = "american_put.csv");

/**
 * @brief set in.disc to a default discretization scheme
 * 
 * @param in the input struct to pricing functions
 * @return UniformDiscretization an discretization for pde solving
 */
UniformDiscretization default_UniformDiscretization(price_inputs in,size_t M, size_t N);

/**
 * @brief Get the exercise boundary from the grid prices of an american option (Horng Tien part IV)
 * 
 * @param in the pricing parameters
 * @param pde_output_file the file where the finite differences for the pde (with given parameters) were stored
 * @param exercise_output_file the file where we wish to store the exercise boundaries
 * @attention if an exercise boundary could not be found, it is set to -1
 */
void get_exercise_boundaries_from_pde(price_inputs in, std::string pde_file, std::string exercise_output_file);


// "BONUS" functions

// Generates put prices from call prices by put call parity ;
// gives pretty good results, though values for big spots may get negative
// (but very close to 0)
// Allowed to make a quick "check" of our pde put prices
void getPutPricesFromCallFile(std::string call_prices_file, double rate, double strike, std::string output_file);

#endif