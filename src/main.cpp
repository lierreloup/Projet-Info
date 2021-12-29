#include "../include/option.h"
#include "../include/pde.h"
#include "../include/fdm.h"
#include "../include/pricers.h"
#include "../include/greeks.h"
#include "../include/vba_interface.h"
#include "../include/graph_builder.h"
#include <iostream>

int main(int argc, char **argv) {

  //argv[1] is the name of a file which contains arguments
  Input input = get_params_from_file(argv[1]);
  
  std::cout << "\nThe gamma is :\n" << gamma_option(input.strike, input.time_to_maturity, input.strike, input.rate, input.volatility, "european_call.csv", "european_call") << "\n";
  
  //TODO : in file vba_interface, define a function which selects the right pricer according to vba input
 
   Price_graph(input.spot, input.time_to_maturity, input.strike, input.rate, input.volatility, "european_call.csv", "european_call");
  
  double price = price_european_call(
    input.strike
    , input.time_to_maturity
    , input.strike
    , input.rate
    , input.volatility
  );

  Output output;
  output.price = price;
  output.delta = 0xbadc0de;
  create_output_file("out_eu_call", output);

  return 0;
}
