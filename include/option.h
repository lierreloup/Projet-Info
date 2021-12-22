#ifndef __VANILLA_OPTION_H
#define __VANILLA_OPTION_H

#include "../include/payoff.h"

class VanillaOption {
 public:

  /* TODO : completely remove payoff objects
  */
  PayOff* pay_off;

  double K;
  double r;
  double T;
  double sigma;

  //double (*option_value_for_0_spot)(double time);
  //double (*option_value_for_infinite_spot) (double time);
  
  virtual double option_price_for_0_spot(double time) const = 0;
  virtual double option_price_for_big_spot(double time, double spot) const = 0;
  virtual double option_price_at_maturity(double spot) const = 0;

  VanillaOption();
  VanillaOption(double _K, double _r, double _T, 
                double _sigma, PayOff* _pay_off);
  
};

class EuropeanCallOption : public VanillaOption {
 public:

  double option_price_for_0_spot(double time) const;
  double option_price_for_big_spot(double time, double spot) const;
  double option_price_at_maturity(double spot) const;

  EuropeanCallOption();
  EuropeanCallOption(double _K, double _r, double _T, 
                double _sigma);
  
};

class AmericanCallOption : public VanillaOption {
 public:

  double option_price_for_0_spot(double time) const;
  double option_price_for_big_spot(double time, double spot) const;
  double option_price_at_maturity(double spot) const;

  AmericanCallOption();
  AmericanCallOption(double _K, double _r, double _T, 
                double _sigma);
  
};

#endif