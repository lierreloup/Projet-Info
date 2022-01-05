#ifndef __GREEK_CPP
#define __GREEK_CPP

#include "../include/greeks.h"
#include "../include/pricers.h"
#include "../include/pde.h"
#include "../include/discretization.h"
#include "../include/vba_interface.h"
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

double delta_option(price_inputs in, std::string output_pde,const char *type_of_option) {

    // precision
  double h = in.spot/10000;
  
    // determine derivative
  if (strcmp(type_of_option,"european_call")==0){
    double Price_S1 = price_european_call(in, output_pde);
    in.spot += h;
    double Price_S2 = price_european_call(in, output_pde);
    return((Price_S2-Price_S1)/h);
      };


  if (strcmp(type_of_option,"european_put")==0){
    double Price_S1 = price_european_put(in, output_pde);

    in.spot+=h;
    double Price_S2 = price_european_put(in, output_pde);
    return((Price_S2-Price_S1)/h);
      };

  
  if (strcmp(type_of_option,"american_call")==0){
    double Price_S1 = price_american_call(in, output_pde);
    in.spot += h;
    double Price_S2 = price_american_call(in, output_pde);
    return((Price_S2-Price_S1)/h);
      };
  
  if (strcmp(type_of_option,"american_put")==0){
    double Price_S1 = price_american_put(in, output_pde);
    in.spot += h;
    double Price_S2 = price_american_put(in, output_pde);
    return((Price_S2-Price_S1)/h);
      };

  return 0.0;

}


double gamma_option(price_inputs in, std::string output_pde,const char *type_of_option) {

    // precision
  double h = in.spot/1000;
  
    // determine derivative
  double delta_d1 = delta_option(in, output_pde, type_of_option);

  in.spot += h;

  double delta_d2 = delta_option(in, output_pde, type_of_option);

  return((delta_d1-delta_d2)/h);

};

double theta_option(price_inputs in, std::string output_pde,const char *type_of_option) {

  // precision
  double h = in.time_to_maturity/100;

  // determine derivative
  if (strcmp(type_of_option,"european_call")==0){
    double Price_S1 = price_european_call(in, output_pde);
    in.time_to_maturity += h;
    double Price_S2 = price_european_call(in, output_pde);
    return((Price_S2-Price_S1)/h);
      };


  if (strcmp(type_of_option,"european_put")==0){
    double Price_S1 = price_european_put(in, output_pde);

    in.time_to_maturity += h;
    double Price_S2 = price_european_put(in, output_pde);
    return((Price_S2-Price_S1)/h);
      };

  
  if (strcmp(type_of_option,"american_call")==0){
    double Price_S1 = price_american_call(in, output_pde);
    in.time_to_maturity += h;
    double Price_S2 = price_american_call(in, output_pde);
    return((Price_S2-Price_S1)/h);
      };


  if (strcmp(type_of_option,"american_put")==0){
    double Price_S1 = price_american_put(in, output_pde);
    in.time_to_maturity += h;
    double Price_S2 = price_american_put(in, output_pde);
    return((Price_S2-Price_S1)/h);
  };
  return 0.0;
}



double rho_option(price_inputs in, std::string output_pde,const char *type_of_option) {

    // precision
  double h = in.rate/100;
    
    // determine derivative
  if (strcmp(type_of_option,"european_call")==0){
    double Price_S1 = price_european_call(in, output_pde);
    in.rate += h;
    double Price_S2 = price_european_call(in, output_pde);
    return((Price_S2-Price_S1)/h);
      };


  if (strcmp(type_of_option,"european_put")==0){
    double Price_S1 = price_european_put(in, output_pde);

    in.rate += h;
    double Price_S2 = price_european_put(in, output_pde);
    return((Price_S2-Price_S1)/h);
      };

  
  if (strcmp(type_of_option,"american_call")==0){
    double Price_S1 = price_american_call(in, output_pde);

    in.rate += h;
    double Price_S2 = price_american_call(in, output_pde);
    return((Price_S2-Price_S1)/h);
      };
  
  if (strcmp(type_of_option,"american_put")==0){
    double Price_S1 = price_american_put(in, output_pde);

    in.rate += h;
    double Price_S2 = price_american_put(in, output_pde);
    return((Price_S2-Price_S1)/h);
  };
  return 0.0;
}



double vega_option(price_inputs in, std::string output_pde, const char *type_of_option) {

    // precision
  double h = in.volatility/100;
    
    // determine derivative
  if (strcmp(type_of_option,"european_call")==0){
    double Price_S1 = price_european_call(in, output_pde);
    in.volatility += h;
    double Price_S2 = price_european_call(in, output_pde);
    return((Price_S2-Price_S1)/h);
      };


  if (strcmp(type_of_option,"european_put")==0){
    double Price_S1 = price_european_put(in, output_pde);
    in.volatility += h;    
    double Price_S2 = price_european_put(in, output_pde);
    return((Price_S2-Price_S1)/h);
      };

  
  if (strcmp(type_of_option,"american_call")==0){
    double Price_S1 = price_american_call(in, output_pde);
    in.volatility += h;
    double Price_S2 = price_american_call(in, output_pde);
    return((Price_S2-Price_S1)/h);
      };
  
  if (strcmp(type_of_option,"american_put")==0){
    double Price_S1 = price_american_put(in, output_pde);
    in.volatility += h;
    double Price_S2 = price_american_put(in, output_pde);
    return((Price_S2-Price_S1)/h);
  };
  return 0.0;
}

#endif