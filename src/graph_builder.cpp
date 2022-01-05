#ifndef __GRAPH_BUILDER_CPP
#define __GRAPH_BUILDER_CPP

#include "../include/pricers.h"
#include "../include/greeks.h"
#include <math.h>
#include <fstream>
#include <string>
#include <iostream>


void Price_graph(double spot, double time_to_maturity, double strike, double rate, double volatility, std::string output_pde, const char* type_of_option){
  std::ofstream myfile;
  myfile.open ("price_graph.csv");
  myfile << "option_price,underlying_price\n";

  price_inputs in;
  in.spot = spot, in.time_to_maturity = time_to_maturity, in.strike = strike, in.rate = rate, in.volatility = volatility;
  UniformDiscretization disc = default_UniformDiscretization(in);
  in.disc = &disc;
  double k=strike/1000;
  for (int i=-50; i<50; i++){
    in.strike += strike + k*i;
    myfile << price_european_call(in, output_pde);
    myfile << ",";
    myfile << strike+i*k;
    myfile << "\n";
  }
  myfile.close();
}

#endif
