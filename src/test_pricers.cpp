#include "../include/pricers.h"
#include "test_utils.cpp"
#include <math.h>
#include <criterion/criterion.h>

//
// EUROPEAN PRICING

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

struct prices_test {
    double underlying,
    strike,
    rate,
    volatility,
    time_to_maturity,
    price,
    theoretical_price
    ;
};



prices_test generate_pricing_test(std::string option_type, double upper_bound) {
    prices_test tst;
    tst.underlying = random_float(upper_bound) + upper_bound / 2, // we add half the uuper bound to get a positive number
    tst.strike = random_float(upper_bound) + upper_bound / 2,
    tst.rate = random_float(upper_bound),
    tst.volatility = random_float(upper_bound) + upper_bound / 2,
    tst.time_to_maturity = random_float(upper_bound)  + upper_bound / 2;

    tst.price = 
        option_type == "call" ?
            price_european_call(
                tst.underlying
                , tst.time_to_maturity
                , tst.strike
                , tst.rate
                , tst.volatility
            )

            : price_european_put(
                tst.underlying
                , tst.time_to_maturity
                , tst.strike
                , tst.rate
                , tst.volatility
            )
    ;
    
    tst.theoretical_price =
        option_type == "call" ?
            BS_call(
                tst.underlying
                , tst.strike
                , tst.rate
                , tst.volatility
                , tst.time_to_maturity
            )

            : BS_put(
                tst.underlying
                , tst.strike
                , tst.rate
                , tst.volatility
                , tst.time_to_maturity
            )
    ;

    return tst;
}

Test(fdm, european_prices) {
    double upper_bound = 100;
    double tol = 0.01;
    size_t count = 10;

    for (size_t i = 0; i < count; i++)
    {
        prices_test tst_call = generate_pricing_test("call", upper_bound);
        cr_assert(abs(tst_call.price - tst_call.theoretical_price) < tol);

        prices_test tst_put = generate_pricing_test("put", upper_bound);
        cr_assert(abs(tst_put.price - tst_put.theoretical_price) < tol);
    }
}