#ifndef __PRICERS_CPP
#define __PRICERS_CPP

#include "../include/option.h"
#include "../include/pde.h"
#include "../include/fdm.h"
#include "../include/pricers.h"
#include <math.h>
#include <fstream>
#include <string>
#include <iostream>
#include <stdexcept>      // std::invalid_argument
#include <limits>

/**
 * @brief finds indices in x_values such that x_values[left_index] <= target <= x_values[right_index] 
 *  * 
 * @param x_values in increasing order
 * @param target 
 * @param left_index 
 * @param right_index
 * 
 * @attention if target = last value of x_values, right_index is set to max size_t
 * @throws std::invalid_argument
 */
void closest(std::vector<double> x_values, double target, size_t & left_index, size_t & right_index) {
    size_t I = x_values.size();

    // look for the proper indices
    for (size_t i = 0; i < I-1; i++) {
        if (x_values[i] <= target && target <= x_values[i+1]) {
            left_index = i, right_index = i+1;
            return;
        }
    }

    // rare case : the last x value is actually the target value
    if (target == x_values[I-1]) {
        left_index = I-1, right_index = std::numeric_limits<size_t>::max();;
    }

    // exception case : no interval was found
    else {
        throw std::invalid_argument("closest : no closest indices were found ; is x_values properly ordered ? is target too big ?");
    }

    return;
}

/**
 * @brief performs linear interpolation
 * 
 * @param x_target 
 * @param x0 
 * @param y0 
 * @param x1 
 * @param y1 
 * @return double the estimated y for x_target
 */
double interpolate(double x_target, double x0, double y0, double x1, double y1) {
    double res = y0 + (x_target - x0) * (y1 - y0) / (x1 - x0);
    std::cout << "x " << x_target << " x0 " << x0 << " x1 "<<x1<<" y0 "<<y0<<" y1 "<<y1<<'\n';
    std::cout << res << '\n';
    return res;
}

/**
 * @brief Auxilliary function for pricers ; solves the PDE and interpolates the price
 * 
 * @param spot 
 * @param solve_pde 
 * @param output_pde 
 * @return double 
 */
double price_aux(double spot, FDMBase & solve_pde, std::string output_pde) {
    std::vector<double> prices_at_time_to_maturity = solve_pde.step_march(output_pde);
    std::vector<double> spot_values(solve_pde.results.x_values);

    // find prices for 2 closest spots
    size_t left_index, right_index;
    closest(spot_values, spot, left_index, right_index);

    // rare case where the last value is the target
    if (right_index == std::numeric_limits<size_t>::max()) return prices_at_time_to_maturity[left_index];

    // interpolate to find price
    double left_price = prices_at_time_to_maturity[left_index]
    , right_price = prices_at_time_to_maturity[right_index];

    return interpolate(spot, spot_values[left_index], left_price, spot_values[right_index], right_price);
}

double price_european_call(price_inputs in, std::string output_pde) {
    // create option
    EuropeanCallOption option(in.strike, in.rate, in.volatility);
    // determine discretization precision
    double t_dom = in.time_to_maturity; //, x_dom = 4 * std::max(in.strike * exp(-rate * time_to_maturity), spot); // is usually enough for the boundary conditions to hold approximately

    std::cout << " N in price " << in.N << '\n';
    // instantiate pde solver
    BSEuroImplicit solve_pde(t_dom, in.N, &option, *in.disc);

    //solve and interpolate
    return price_aux(in.spot, solve_pde, output_pde);
}

double price_european_put(price_inputs in, std::string output_pde) {
    EuropeanPutOption option(in.strike, in.rate, in.volatility);
    double t_dom = in.time_to_maturity; 


    // instantiate pde solver
    BSEuroImplicit solve_pde(t_dom, in.N, &option, *in.disc);

    return price_aux(in.spot, solve_pde, output_pde);
}

double price_american_call(price_inputs in, std::string output_pde) {
    // create option
    AmericanCallOption option(in.strike, in.rate, in.volatility);

    double t_dom = in.time_to_maturity;

    // solve pde and get the last line (which corresponds to prices for all spots and given time_to_maturity)
    BSAmericanImplicit solve_pde(t_dom, in.N, &option, *in.disc);

    return price_aux(in.spot, solve_pde, output_pde);
}


double price_american_put(price_inputs in, std::string output_pde) {
    // create option
    AmericanPutOption option(in.strike, in.rate, in.volatility);

    double t_dom = in.time_to_maturity;

    // solve pde and get the last line (which corresponds to prices for all spots and given time_to_maturity)
    BSAmericanImplicit solve_pde(t_dom, in.N, &option, *in.disc);

    return price_aux(in.spot, solve_pde, output_pde);
}

UniformDiscretization default_UniformDiscretization(price_inputs in) {
    double x_dom = 3 * std::max(in.spot, exp(-in.rate * in.time_to_maturity));
    size_t M = std::min(x_dom * 10, 1e3);

    UniformDiscretization disc(x_dom, M);
    return disc;
}

NonUniformDiscretization default_NonUniformDiscretization(price_inputs in) {
    double x_dom = 3 * std::max(in.spot, exp(-in.rate * in.time_to_maturity));
    size_t M = std::min(static_cast<size_t>(x_dom * 10.), static_cast<size_t>(1e4));

    double disc_center = in.strike
        , c = in.strike / 3.
    ;

    NonUniformDiscretization disc(x_dom, M, disc_center, c);
    return disc;
}

size_t default_N(price_inputs in) {
    double result = std::min(in.time_to_maturity * 10., 1e3);
    std::cout << " in default " << result << '\n';

    return static_cast<size_t>(result);
}








// "BONUS" FUNCTIONS


/**
 * @brief Get the Put Price From Call Price using Put Call parity
 * 
 * @param callPrice 
 * @param spot 
 * @param time 
 * @param rate 
 * @param strike 
 */
double getPutPriceFromCallPrice(double call_price, double spot, double time_to_maturity, double rate, double strike) {
    return call_price - spot + exp(-rate * time_to_maturity) * strike;
}

void getPutPricesFromCallFile(std::string call_prices_file, double rate, double strike, std::string output_file) {
    std::fstream call_prices_stream(call_prices_file, std::ios::in);
    std::ofstream put_out(output_file);

    while(call_prices_stream) {
        double spot, time_to_maturity, call_price;           

        call_prices_stream >> spot >> time_to_maturity >> call_price;    

        double put_price = getPutPriceFromCallPrice(call_price, spot, time_to_maturity, rate, strike);
        put_out << spot << ' ' << time_to_maturity << ' ' << put_price << std::endl;
    }

    put_out.close();
    call_prices_stream.close();
}

#endif