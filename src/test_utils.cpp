#ifndef __TEST_UTILS_CPP
#define __TEST_UTILS_CPP

#include <bits/stdc++.h>
//#include <math.h>


double random_float(double upper_bound) {
    return static_cast <double> (rand()) / (static_cast <double> (RAND_MAX)) * upper_bound - upper_bound / 2;
}

double randmax() {
    return random_float(1000);
}
/**
 * @brief random vector of size n
 * 
 * @param n 
 * @return std::vector<double> 
 */
std::vector<double> rand_vals(size_t n) {
    std::vector<double> v(n);
    generate(v.begin(), v.end(), randmax);

    return v;
}

double normalCDF(double value)
{
   return 0.5 * erfc(-value * M_SQRT1_2);
}

double d(double plus_or_minus, double underlying, double strike, double rate, double volatility, double time_to_maturity) {
    double top_left = log(underlying/strike),
    top_right = (rate + plus_or_minus * volatility * volatility / 2) * time_to_maturity,
    bottom = volatility * sqrt(time_to_maturity);

    return (top_left + top_right) / bottom;
}

double BS_call(double underlying, double strike, double rate, double volatility, double time_to_maturity) {
    double d1 = d(+1, underlying, strike, rate, volatility, time_to_maturity),
    d2 = d(-1, underlying, strike, rate, volatility, time_to_maturity),
    left = underlying * normalCDF(d1),
    right = strike * exp(-rate * time_to_maturity) * normalCDF(d2);

    return left - right;
}

double BS_put(double underlying, double strike, double rate, double volatility, double time_to_maturity) {
    double d1 = d(+1, underlying, strike, rate, volatility, time_to_maturity),
    d2 = d(-1, underlying, strike, rate, volatility, time_to_maturity),
    left = strike * exp(-rate * time_to_maturity) * normalCDF(-d2),
    right = underlying * normalCDF(-d1);

    return left - right;
}


#endif