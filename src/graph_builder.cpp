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
  double k=strike/1000;
  for (int i=-50; i<50; i++){
    myfile << price_european_call(spot, time_to_maturity, strike+k*i, rate, volatility, output_pde);
    myfile << ",";
    myfile << strike+i*k;
    myfile << "\n";
  }
  myfile.close();
}

#endif
