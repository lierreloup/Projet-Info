#include "../include/payoff.h"
#include "../include/option.h"
#include "../include/pde.h"
#include "../include/fdm.h"
#include <math.h>
#include <fstream>
#include <string>
#include <iostream>

/**
 * @brief Get the Put Price From Call Price using Put Call parity
 * 
 * @param callPrice 
 * @param spot 
 * @param time 
 * @param rate 
 * @param strike 
 */
double getPutPriceFromCallPrice(double call_price, double spot, double time, double rate, double strike) {
    return call_price - spot + exp(-rate * time) * strike;
}

void getPutPricesFromCallFile(std::string call_prices_file, double spot, double time, double rate, double strike, std::string output_file) {
    std::fstream call_prices_stream(call_prices_file, std::ios::in);
    std::ofstream put_out(output_file);

    while(call_prices_stream) {
        double spot, time, call_price;           

        for (size_t i = 0; i < 3; i++) {
            call_prices_stream >> spot >> time >> call_price;    
        }

        double put_price = getPutPriceFromCallPrice(call_price, spot, time, rate, strike);
        put_out << spot << ' ' << time << ' ' << put_price << std::endl;
    }
}