#include "../include/payoff.h"
#include "../include/option.h"
#include "../include/pde.h"
#include "../include/fdm.h"
#include "../include/pricers.h"

#include <iostream>

int main(int argc, char **argv) {
  // Create the option parameters
  double K = 0.5;  // Strike price
  double r = 0.05;   // Risk-free rate (5%)
  double v = 0.2;    // Volatility of the underlying (20%)
  double T = 1.00;    // One year until expiry

  // FDM discretisation parameters
  double x_dom = 1;       // Spot goes from [0.0, 1.0]
  unsigned long M = 20;
  double t_dom = T;         // Time period as for the option
  unsigned long N = 20;     

  // Create the PayOff and Option objects
  //PayOff* pay_off_call = new PayOffCall(K);
  std::cout << "in main K is " << K << std::endl;
  EuropeanCallOption* call_option = new EuropeanCallOption(K, r, T, v);

  // Create the PDE and FDM objects
  BlackScholesPDE* bs_pde = new BlackScholesPDE(call_option);
  FDMEulerImplicit fdm_euler(x_dom, M, t_dom, N, bs_pde);

  //AmericanOptionParameters call_params;
  //call_params.no_early_exercise_pde = bs_pde;

  //PriceAmericanOption price(x_dom, M, t_dom, N, call_params);

  // Run the FDM solver
  fdm_euler.step_march("call_fdm_implicit_backwards.csv");
  //price.step_march("america.csv");

  // Compute Put Price
  //getPutPricesFromCallFile("call_fdm_implicit.csv", T, r, K, "put_fdm_implicit");

  // Delete the PDE, PayOff and Option objects
  delete bs_pde;
  delete call_option;
  //delete pay_off_call;

  return 0;
}