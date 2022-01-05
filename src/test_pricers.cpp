#include "../include/pricers.h"
#include "test_utils.cpp"
#include <math.h>
#include <criterion/criterion.h>

//
// EUROPEAN PRICING



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

void print_prices_test(prices_test tst) {
    std::cout << "Underlying : " << tst.underlying;
    std::cout << "\nTime to maturity : " << tst.time_to_maturity;
    std::cout << "\nStrike : " << tst.strike;
    std::cout << "\nRate : " << tst.rate;
    std::cout << "\nVolatility : " << tst.volatility;

    std::cout << "\nPrice : " << tst.price;
    std::cout << "\nTheoretical price : " << tst.theoretical_price;
    std::cout << '\n';
}



prices_test generate_pricing_test(std::string option_type, double upper_bound) {
    prices_test tst;
    tst.underlying = random_float(upper_bound) + upper_bound / 2, // we add half the upper bound to get a positive number
    tst.strike = random_float(upper_bound) + upper_bound / 2,
    tst.rate = random_float(upper_bound) / 100,
    tst.volatility = (random_float(upper_bound) + upper_bound / 2) / 100,
    tst.time_to_maturity = random_float(upper_bound)  + upper_bound / 2;

    //print_prices_test(tst);

    price_inputs in;
    in.spot = tst.underlying;
    in.strike = tst.strike;
    in.rate = tst.rate;
    in.volatility = tst.volatility;
    in.time_to_maturity = tst.time_to_maturity;
    in.N = default_N(in);
    NonUniformDiscretization disc = default_NonUniformDiscretization(in);
    in.disc = &disc;

    //print_vector(disc.get_grid());
    tst.price = 
        option_type == "call" ?
            price_european_call(in)

            : price_european_put(in)
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

    //print_prices_test(tst);
    return tst;
}

Test(fdm, european_prices) {
    double upper_bound = 100;
    double tol = 0.01;
    size_t count = 10;

    for (size_t i = 0; i < count; i++)
    {
        prices_test tst_call = generate_pricing_test("call", upper_bound);
        std::cout << "call \n";
        print_prices_test(tst_call);
        cr_assert(abs(tst_call.price - tst_call.theoretical_price) < tol);


        prices_test tst_put = generate_pricing_test("put", upper_bound);
        std::cout << "put \n";
        print_prices_test(tst_put);
        cr_assert(abs(tst_put.price - tst_put.theoretical_price) < tol);
    }
}