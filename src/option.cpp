#ifndef __VANILLA_OPTION_CPP
#define __VANILLA_OPTION_CPP

#include "../include/option.h"
#include <math.h>
#include <iostream>

VanillaOption::VanillaOption(double _K, double _r,
                             double _sigma) : 
  K(_K), r(_r), sigma(_sigma) {}

/**
 *
 * EUROPEAN OPTIONS
 *
*/

// CALL

EuropeanCallOption::EuropeanCallOption(double _K, double _r,
                             double _sigma) 
: EuropeanOption(_K, _r, _sigma) {}


double european_call_price_for_0_spot() {
  return 0;
}

double EuropeanCallOption::option_price_for_0_spot(double time_to_maturity) const {
  return european_call_price_for_0_spot();
}

double european_call_price_for_big_spot(double time_to_maturity, double spot, double strike, double rate) {
  return spot - strike * exp(-rate * time_to_maturity);
}

double EuropeanCallOption::option_price_for_big_spot(double time_to_maturity, double spot) const {
  return european_call_price_for_big_spot(
    time_to_maturity
    , spot
    , this->K
    , this->r
  );
}

double european_call_payoff(double spot, double strike) {
  double diff = spot - strike;
  return diff > 0 ? diff : 0;
}

double EuropeanCallOption::option_price_at_maturity(double spot) const {
  return european_call_payoff(spot, this->K);
}





//____________________________________________________________________________
// PUT

EuropeanPutOption::EuropeanPutOption(double _K, double _r,
                             double _sigma) 
: EuropeanOption(_K, _r, _sigma) {}


double european_put_price_for_0_spot(double time_to_maturity, double strike, double rate) {
  return strike * exp(-rate * time_to_maturity);
}

double EuropeanPutOption::option_price_for_0_spot(double time_to_maturity) const {
  return european_put_price_for_0_spot(
    time_to_maturity
    , this->K
    , this->r
    );
}

double european_put_price_for_big_spot() {
  return 0;
}

double EuropeanPutOption::option_price_for_big_spot(double time_to_maturity, double spot) const {
  return european_put_price_for_big_spot();
}

double european_put_payoff(double spot, double strike) {
  double diff = strike - spot;
  return diff > 0 ? diff : 0;
}

double EuropeanPutOption::option_price_at_maturity(double spot) const {
  return european_put_payoff(spot, this->K);
}








//____________________________________________________________________________
//____________________________________________________________________________

/**
 *
 * AMERICAN OPTIONS
 *
*/

// CALL

AmericanCallOption::AmericanCallOption(double _K, double _r, double _sigma)
: AmericanOption(_K, _r, _sigma) {}

double AmericanCallOption::option_price_for_0_spot(double time_to_maturity) const {
  // American and European calls have the same boundary condition for 0 spot
  return european_call_price_for_0_spot();
}

double american_call_price_for_big_spot(double time_to_maturity, double spot, double strike) {
  return spot - strike;
}

double AmericanCallOption::option_price_for_big_spot(double time_to_maturity, double spot) const {
  return american_call_price_for_big_spot(time_to_maturity, spot, this->K);
}

double AmericanCallOption::option_price_at_maturity(double spot) const {
  return european_call_payoff(spot, this->K); //same as European
}








// PUT
AmericanPutOption::AmericanPutOption(double _K, double _r, double _sigma)
: AmericanOption(_K, _r, _sigma) {}

double american_put_price_for_0_spot(double strike) {
  return strike;
}

double AmericanPutOption::option_price_for_0_spot(double time_to_maturity) const {
  return american_put_price_for_0_spot(this->K);
}

double AmericanPutOption::option_price_for_big_spot(double time_to_maturity, double spot) const {
  return european_put_price_for_big_spot(); //same as European
}

double AmericanPutOption::option_price_at_maturity(double spot) const {
  return european_put_payoff(spot, this->K); //same as European
}




#endif
