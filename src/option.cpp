#ifndef __VANILLA_OPTION_CPP
#define __VANILLA_OPTION_CPP

#include "../include/option.h"
#include <math.h>
#include "../include/payoff.h"
#include <iostream>

VanillaOption::VanillaOption() {}

VanillaOption::VanillaOption(double _K, double _r, double _T, 
                             double _sigma, PayOff* _pay_off) : 
  K(_K), r(_r), T(_T), sigma(_sigma), pay_off(_pay_off) {
  
  std::cout << "K init in vanilla" << this->K << std::endl;

}



EuropeanCallOption::EuropeanCallOption(double _K, double _r, double _T, 
                             double _sigma) 
: VanillaOption(_K, _r, _T, _sigma, new PayOffCall(_K)) {

  std::cout << "K init in european" << this->K << std::endl;
}


double european_call_price_for_0_spot(double time_to_maturity) {
  return 0;
}

double EuropeanCallOption::option_price_for_0_spot(double time_to_maturity) const {
  return european_call_price_for_0_spot(time_to_maturity);
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
  //std::cout << "spot " << spot << " strike " << strike << std::endl;
  double diff = spot - strike;
  return diff > 0 ? diff : 0;
}

double EuropeanCallOption::option_price_at_maturity(double spot) const {
  return european_call_payoff(spot, this->K);
}

AmericanCallOption::AmericanCallOption(double _K, double _r, double _T, double _sigma)
: VanillaOption(_K, _r, _T, _sigma, new PayOffCall(_K)) {}

double AmericanCallOption::option_price_for_0_spot(double time_to_maturity) const {
  // american and european options have the same boundary condition for 0 spot
  return european_call_price_for_0_spot(time_to_maturity);
}

double american_call_price_for_big_spot(double time_to_maturity, double spot, double strike) {
  return spot - strike;
}

double AmericanCallOption::option_price_for_big_spot(double time_to_maturity, double spot) const {
  return american_call_price_for_big_spot(time_to_maturity, spot, this->K);
}

double AmericanCallOption::option_price_at_maturity(double spot) const {
  return european_call_payoff(spot, this->K);
}


#endif