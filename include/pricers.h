#ifndef __PRICERS_H
#define __PRICERS_H

#include <string>


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
double price_european_call(double spot, double time_to_maturity, double strike, double rate, double volatility, std::string output_pde = "european_call.csv");

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
double price_european_put(double spot, double time_to_maturity, double strike, double rate, double volatility, std::string output_pde = "european_put.csv");

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
double price_american_call(double spot, double time_to_maturity, double strike, double rate, double volatility, std::string output_pde = "american_call.csv");

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
double price_american_put(double spot, double time_to_maturity, double strike, double rate, double volatility, std::string output_pde = "american_put.csv");



// "BONUS" functions

// Generates put prices from call prices by put call parity ;
// gives pretty good results, though values for big spots may get negative
// (but very close to 0)
// Allowed to make a quick "check" of our pde put prices
void getPutPricesFromCallFile(std::string call_prices_file, double rate, double strike, std::string output_file);

#endif