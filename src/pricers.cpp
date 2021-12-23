#include "../include/option.h"
#include "../include/pde.h"
#include "../include/fdm.h"
#include "../include/pricers.h"
#include <math.h>
#include <fstream>
#include <string>
#include <iostream>
#include <stdexcept>      // std::invalid_argument
/**
 * @brief Get the Put Price From Call Price using Put Call parity
 * 
 * @param callPrice 
 * @param spot 
 * @param time 
 * @param rate 
 * @param strike 
 */
double getPutPriceFromCallPrice(double call_price, double spot, double time, double maturity, double rate, double strike) {
    return call_price - spot + exp(-rate * (maturity - time)) * strike;
}

void getPutPricesFromCallFile(std::string call_prices_file, double maturity, double rate, double strike, std::string output_file) {
    std::fstream call_prices_stream(call_prices_file, std::ios::in);
    std::ofstream put_out(output_file);

    while(call_prices_stream) {
        double spot, time, call_price;           

        call_prices_stream >> spot >> time >> call_price;    

        double put_price = getPutPriceFromCallPrice(call_price, spot, time, maturity, rate, strike);
        put_out << spot << ' ' << time << ' ' << put_price << std::endl;
    }

    put_out.close();
}

/**
 * @brief finds indices in x_values such that x_values[left_index] <= target <= x_values[right_index] 
 *  * 
 * @param x_values in increasing order
 * @param target 
 * @param left_index 
 * @param right_index
 * 
 * @attention if target = last value of x_values, right_index = -1
 * @throws std::invalid_argument
 */
void closest(std::vector<double> x_values, double target, size_t & left_index, size_t & right_index) {
    size_t I = x_values.size();
        std::cout << "my size is " << I << std::endl;

    for (size_t i = 0; i < I-1; i++) {
        std::cout << "watcha doin" << x_values[i] << std::endl;
        if (x_values[i] <= target && target <= x_values[i+1]) {
            left_index = i, right_index = i+1;
            return;
        }
    }

    if (target == x_values[I-1]) {
        left_index = I-1, right_index = -1;
    }
    else {
        throw std::invalid_argument("closest : no closest indices were found ; is x_values properly ordered ? is target too big ?");
    }

    return;
}

double interpolate(double x_target, double x0, double y0, double x1, double y1) {
    return y0 + (x_target - x0) * (y1 - y0) / (x1 - x0);
}

double price_aux(double spot, FDMBase & solve_pde, std::string output_pde) {
    std::vector<double> prices_at_time_to_maturity = solve_pde.step_march(output_pde);
    std::vector<double> spot_values = solve_pde.results.x_values;

    // find prices for 2 closest spots
    size_t left_index, right_index;
    closest(spot_values, spot, left_index, right_index);

    // rare case where the last value is the target
    if (right_index == -1) return prices_at_time_to_maturity[left_index];

    // interpolate to find price
    double left_price = prices_at_time_to_maturity[left_index]
    , right_price = prices_at_time_to_maturity[right_index];

    return interpolate(spot, spot_values[left_index], left_price, spot_values[right_index], right_price);
}


double price_european_call(double spot, double time_to_maturity, double strike, double rate, double volatility, std::string output_pde) {
    // create option
    EuropeanCallOption option(strike, rate, volatility);

    // determine discretization precision
    double t_dom = time_to_maturity, x_dom = 4 * strike; // 4 * strike is usually enough for the boundary conditions to hold approximately
    size_t M = 100, N = 100;

    // solve pde and get the last line (which corresponds to prices for all spots and given time_to_maturity)
    BSEuroImplicit solve_pde(x_dom, M, t_dom, N, &option);

    return price_aux(spot, solve_pde, output_pde);
}

double price_american_call(double spot, double time_to_maturity, double strike, double rate, double volatility, std::string output_pde) {
    // create option
    AmericanCallOption option(strike, rate, volatility);

    // determine discretization precision
    double t_dom = time_to_maturity, x_dom = 4 * strike; // 4 * strike is usually enough for the boundary conditions to hold approximately
    size_t M = 100, N = 100;

    // solve pde and get the last line (which corresponds to prices for all spots and given time_to_maturity)
    BSAmericanImplicitUniform solve_pde(x_dom, M, t_dom, N, &option);

    return price_aux(spot, solve_pde, output_pde);
}

/**
 * @brief Set the params from file
 * 
 * @param filename path to file ; the file should contain each parameter separated by a new line (or a space)
 * @param spot 
 * @param time_to_maturity 
 * @param strike 
 * @param rate 
 * @param volatility 
 */
void set_params_from_file(std::string filename, double & spot, double & time_to_maturity, double & strike, double & rate, double & volatility) {
    std::fstream params_stream(filename, std::ios::in);

    while(params_stream) {
        params_stream >> spot >> time_to_maturity >> strike >> rate >> volatility;    
    }

    params_stream.close();
}


void create_output_file(std::string filename, Output output) {
    std::ofstream put_out(filename);
    put_out << output.price << '\n';
    put_out << output.delta << '\n';

    put_out.close();
}