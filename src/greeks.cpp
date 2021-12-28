#ifndef __GREEK_CPP
#define __GREEK_CPP

#include "../include/greeks.h"
#include "../include/pricers.h"
#include <math.h>
#include <fstream>
#include <string>
#include <iostream>
#include <string.h>

/**
 * @brief calculate the delta of an option
 *  * 
 * @param spot value
 * @param time_to_maturity
 * @param strike
 * @param rate
 * @param volatility
 * @param output_pde
 * @param type_of_option 
 *
 */


double delta_option(double spot, double time_to_maturity, double strike, double rate, double volatility, std::string output_pde, const char* type_of_option) {

    // precision
  double h = spot/10000;
  
    // determine derivative
  if (strcmp(type_of_option,"european_call")==0){
    double Price_S1 = price_european_call(spot,time_to_maturity, strike, rate, volatility, output_pde);
    double Price_S2 = price_european_call(spot+h,time_to_maturity, strike, rate, volatility, output_pde);
    return((Price_S2-Price_S1)/h);
      };


  if (strcmp(type_of_option,"european_put")==0){
    double Price_S1 = price_european_put(spot,time_to_maturity, strike, rate, volatility, output_pde);

    
    double Price_S2 = price_european_put(spot+h,time_to_maturity, strike, rate, volatility, output_pde);
    return((Price_S2-Price_S1)/h);
      };

  
  if (strcmp(type_of_option,"american_call")==0){
    double Price_S1 = price_american_call(spot,time_to_maturity, strike, rate, volatility, output_pde);
    double Price_S2 = price_american_call(spot+h,time_to_maturity, strike, rate, volatility, output_pde);
    return((Price_S2-Price_S1)/h);
      };
  
  if (strcmp(type_of_option,"american_put")==0){
    double Price_S1 = price_american_put(spot,time_to_maturity, strike, rate, volatility, output_pde);
    double Price_S2 = price_american_put(spot+h,time_to_maturity, strike, rate, volatility, output_pde);
    return((Price_S2-Price_S1)/h);
      };

  return 0.0;

}


double gamma_option(double spot, double time_to_maturity, double strike, double rate, double volatility, std::string output_pde, const char* type_of_option) {

    // precision
  double h = spot/100000;
  
    // determine derivative
  double delta_d1 = delta_option(spot, time_to_maturity, strike, rate, volatility, output_pde, type_of_option);

  double delta_d2 = delta_option(spot+h, time_to_maturity, strike, rate, volatility, output_pde, type_of_option);

  return((delta_d1-delta_d2)/h);

};

double theta_option(double spot, double time_to_maturity, double strike, double rate, double volatility, std::string output_pde, const char* type_of_option) {

    // precision
  double h = time_to_maturity/10000;
  
    // determine derivative
  if (strcmp(type_of_option,"european_call")==0){
    double Price_S1 = price_european_call(spot,time_to_maturity, strike, rate, volatility, output_pde);
    double Price_S2 = price_european_call(spot,time_to_maturity+h, strike, rate, volatility, output_pde);
    return((Price_S2-Price_S1)/h);
      };


  if (strcmp(type_of_option,"european_put")==0){
    double Price_S1 = price_european_put(spot,time_to_maturity, strike, rate, volatility, output_pde);

    
    double Price_S2 = price_european_put(spot,time_to_maturity+h, strike, rate, volatility, output_pde);
    return((Price_S2-Price_S1)/h);
      };

  
  if (strcmp(type_of_option,"american_call")==0){
    double Price_S1 = price_american_call(spot,time_to_maturity, strike, rate, volatility, output_pde);
    double Price_S2 = price_american_call(spot,time_to_maturity+h, strike, rate, volatility, output_pde);
    return((Price_S2-Price_S1)/h);
      };
  
  if (strcmp(type_of_option,"american_put")==0){
    double Price_S1 = price_american_put(spot,time_to_maturity, strike, rate, volatility, output_pde);
    double Price_S2 = price_american_put(spot,time_to_maturity+h, strike, rate, volatility, output_pde);
    return((Price_S2-Price_S1)/h);
  };
}



double rho_option(double spot, double time_to_maturity, double strike, double rate, double volatility, std::string output_pde, const char* type_of_option) {

    // precision
  double h = rate/10000;
  
    // determine derivative
  if (strcmp(type_of_option,"european_call")==0){
    double Price_S1 = price_european_call(spot,time_to_maturity, strike, rate, volatility, output_pde);
    double Price_S2 = price_european_call(spot,time_to_maturity, strike, rate+h, volatility, output_pde);
    return((Price_S2-Price_S1)/h);
      };


  if (strcmp(type_of_option,"european_put")==0){
    double Price_S1 = price_european_put(spot,time_to_maturity, strike, rate, volatility, output_pde);

    
    double Price_S2 = price_european_put(spot,time_to_maturity, strike, rate+h, volatility, output_pde);
    return((Price_S2-Price_S1)/h);
      };

  
  if (strcmp(type_of_option,"american_call")==0){
    double Price_S1 = price_american_call(spot,time_to_maturity, strike, rate, volatility, output_pde);
    double Price_S2 = price_american_call(spot,time_to_maturity, strike, rate+h, volatility, output_pde);
    return((Price_S2-Price_S1)/h);
      };
  
  if (strcmp(type_of_option,"american_put")==0){
    double Price_S1 = price_american_put(spot,time_to_maturity, strike, rate, volatility, output_pde);
    double Price_S2 = price_american_put(spot,time_to_maturity, strike, rate+h, volatility, output_pde);
    return((Price_S2-Price_S1)/h);
  };
}



double vega_option(double spot, double time_to_maturity, double strike, double rate, double volatility, std::string output_pde, const char* type_of_option) {

    // precision
  double h = volatility/10000;
  
    // determine derivative
  if (strcmp(type_of_option,"european_call")==0){
    double Price_S1 = price_european_call(spot,time_to_maturity, strike, rate, volatility, output_pde);
    double Price_S2 = price_european_call(spot,time_to_maturity, strike, rate, volatility+h, output_pde);
    return((Price_S2-Price_S1)/h);
      };


  if (strcmp(type_of_option,"european_put")==0){
    double Price_S1 = price_european_put(spot,time_to_maturity, strike, rate, volatility, output_pde);

    
    double Price_S2 = price_european_put(spot,time_to_maturity, strike, rate, volatility+h, output_pde);
    return((Price_S2-Price_S1)/h);
      };

  
  if (strcmp(type_of_option,"american_call")==0){
    double Price_S1 = price_american_call(spot,time_to_maturity, strike, rate, volatility, output_pde);
    double Price_S2 = price_american_call(spot,time_to_maturity, strike, rate, volatility+h, output_pde);
    return((Price_S2-Price_S1)/h);
      };
  
  if (strcmp(type_of_option,"american_put")==0){
    double Price_S1 = price_american_put(spot,time_to_maturity, strike, rate, volatility, output_pde);
    double Price_S2 = price_american_put(spot,time_to_maturity, strike, rate, volatility+h, output_pde);
    return((Price_S2-Price_S1)/h);
  };
}
#endif
