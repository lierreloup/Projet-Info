#ifndef __GRAPH_BUILDER_CPP
#define __GRAPH_BUILDER_CPP

#include "../include/pricers.h"
#include "../include/greeks.h"
#include <math.h>
#include <fstream>
#include <string>
#include <iostream>


void Price_graph(price_inputs in, std::string output_pde, const char* type_of_option){
  std::ofstream myfile;
  myfile.open ("price_graph.csv");
  myfile << "option_price,underlying_price\n";


  double k=in.strike/1000;
  for (int i=-50; i<50; i++){
    in.strike += in.strike + k*i;
    myfile << price_european_call(in, output_pde);
    myfile << ",";
    myfile << in.strike+i*k;
    myfile << "\n";
  }
  myfile.close();
}

#endif