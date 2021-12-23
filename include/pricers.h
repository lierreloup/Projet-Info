// TODO : rename this module as "utils" in .h, .cpp, Makefile

#include <string>

void getPutPricesFromCallFile(std::string call_prices_file, double maturity, double rate, double strike, std::string output_file);

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
 * @brief Reads the vba input (vba must write its input in a file) and assigns it to the proper parameters
 * 
 * @param filename input file
 * @param spot 
 * @param time_to_maturity 
 * @param strike 
 * @param rate 
 * @param volatility 
 */
void set_params_from_file(std::string filename, double & spot, double & time_to_maturity, double & strike, double & rate, double & volatility);

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
