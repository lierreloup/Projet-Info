// TODO : rename this module as "utils" in .h, .cpp, Makefile

#include <string>

void getPutPricesFromCallFile(std::string call_prices_file, double maturity, double rate, double strike, std::string output_file);

double price_european_call(double spot, double time_to_maturity, double strike, double rate, double volatility, std::string output_pde = "european_call.csv");

double price_american_call(double spot, double time_to_maturity, double strike, double rate, double volatility, std::string output_pde = "american_call.csv");

void set_params_from_file(std::string filename, double & spot, double & time_to_maturity, double & strike, double & rate, double & volatility);

struct Output {
    double price;
    double delta;
    //etc
};

void create_output_file(std::string filename, Output output);
