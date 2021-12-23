#ifndef __VANILLA_OPTION_H
#define __VANILLA_OPTION_H

/**
 * @brief  abstract class for vanilla options, with parameters required for PDE pricing
 * 
 */
class VanillaOption {
 public:

  double K;
  double r;
  double sigma;

  //double (*option_value_for_0_spot)(double time);
  //double (*option_value_for_infinite_spot) (double time);
  
  /**
   * @brief returns P(time, 0) where P(t, S) is the price of the option, S the spot, t the time to maturity
   * 
   * @param time 
   * @return double 
   */
  virtual double option_price_for_0_spot(double time) const = 0;

  /**
   * @brief returns approximate P(time, spot) where P(t, S) is the price of the option, S the spot (close to "infinity"), t the time to maturity
   * 
   * @param time time to maturity
   * @param spot close to "infinity" (e.g 4 * strike)
   * @return double 
   */
  virtual double option_price_for_big_spot(double time, double spot) const = 0;

  /**
   * @brief returns P(0, spot) where P(t, S) is the price of the option, S the spot, t the time to maturity

   * 
   * @param spot 
   * @return double (e.g the payoff)
   */
  virtual double option_price_at_maturity(double spot) const = 0;

  VanillaOption(double _K, double _r, 
                double _sigma);
  
};

// Abstract classes for european and american options
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

//Concrete classes
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