#include "../include/option.h"
#include "../include/pde.h"
#include "../include/fdm.h"
#include "../include/pricers.h"
#include "../include/greeks.h"
#include "../include/vba_interface.h"
#include "../include/graph_builder.h"
#include "../include/discretization.h"
#include <iostream>

int main(int argc, char **argv) {

  //argv[1] is the name of a file which contains arguments
  Input input = get_params_from_file(argv[1]);

  
  //std::cout << "\nThe gamma is :\n" << gamma_option(input.strike, input.time_to_maturity, input.strike, input.rate, input.volatility, "european_call.csv", "european_call") << "\n";
  
  //TODO : in file vba_interface, define a function which selects the right pricer according to vba input
 
  //Price_graph(input.spot, input.time_to_maturity, input.strike, input.rate, input.volatility, "european_call.csv", "european_call");
  

  double price = price_option(input);

  const char* type_of_option = input.option_type.c_str();

  std::remove("../output");
  std::remove("../output.csv");
  Output output;
  output.price = price;
  store_greeks( output, input);
  create_output_file(argv[2], output);

Price_graph(input, input.option_type, type_of_option);
  delta_graph(input,  input.option_type,  type_of_option);

  return 0;
}
