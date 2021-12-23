#ifndef __VANILLA_OPTION_H
#define __VANILLA_OPTION_H

class VanillaOption {
 public:

  double K;
  double r;
  double sigma;

  //double (*option_value_for_0_spot)(double time);
  //double (*option_value_for_infinite_spot) (double time);
  
  virtual double option_price_for_0_spot(double time) const = 0;
  virtual double option_price_for_big_spot(double time, double spot) const = 0;
  virtual double option_price_at_maturity(double spot) const = 0;

  // TODO : why default constructor ?
  //VanillaOption();
  VanillaOption(double _K, double _r, 
                double _sigma);
  
};


class EuropeanOption : public VanillaOption {
    protected:
    EuropeanOption(double _K, double _r, 
                double _sigma)
    : VanillaOption(_K, _r, _sigma) {}
};

class AmericanOption : public VanillaOption {
    protected:
    AmericanOption(double _K, double _r, 
                double _sigma)
    : VanillaOption(_K, _r, _sigma) {}
};

class EuropeanCallOption : public EuropeanOption {
 public:

  double option_price_for_0_spot(double time) const;
  double option_price_for_big_spot(double time, double spot) const;
  double option_price_at_maturity(double spot) const;

  //EuropeanCallOption();
  EuropeanCallOption(double _K, double _r,
                double _sigma);
  
};

class EuropeanPutOption : public EuropeanOption {
 public:

  double option_price_for_0_spot(double time) const;
  double option_price_for_big_spot(double time, double spot) const;
  double option_price_at_maturity(double spot) const;

  //EuropeanPutOption();
  EuropeanPutOption(double _K, double _r,
                double _sigma);
  
};

/*
union EuroOption {
    EuropeanCallOption call;
    EuropeanPutOption put;
};
*/

class AmericanCallOption : public AmericanOption {
 public:

  double option_price_for_0_spot(double time) const;
  double option_price_for_big_spot(double time, double spot) const;
  double option_price_at_maturity(double spot) const;

  //AmericanCallOption();
  AmericanCallOption(double _K, double _r,
                double _sigma);
  
};

#endif