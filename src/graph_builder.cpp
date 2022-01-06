#ifndef __GRAPH_BUILDER_CPP
#define __GRAPH_BUILDER_CPP

#include "../include/pricers.h"
#include "../include/greeks.h"
#include <math.h>
#include <fstream>
#include <string>
#include <iostream>

void Price_graph(Input in, std::string output_pde, const char* type_of_option){
  std::remove("price_graph.csv");
  std::ofstream myfile;
  myfile.open ("price_graph.csv");
  myfile << "option_price,underlying_price\n";

  double k=in.spot/100;
  for (int i=-50; i<50; i++)
  {  
    in.spot += k;
    myfile << std::to_string(price_option(in)) + "," + std::to_string(in.spot+i*k);
    myfile << "\n";
  }
  myfile.close();
}

void delta_graph(Input input, std::string output_pde,  const char* type_of_option){
  
  price_inputs in;
    in.spot = input.spot, in.time_to_maturity = input.time_to_maturity, in.strike = input.strike, in.rate = input.rate, in.volatility = input.volatility;
    UniformDiscretization disc = default_UniformDiscretization(in, input.M, input.N);
    in.disc = &disc;
  
  
  std::remove("delta_graph.csv");
  std::ofstream myfile;
  myfile.open ("delta_graph.csv");
  myfile << "option_price,underlying_price\n";

  double k=in.spot/100;
  for (int i=-30; i<30; i++){
    in.spot += k;
    myfile << delta_option(in, output_pde,type_of_option)<< "," << in.spot+i*k;
    myfile << "\n";
  }
  myfile.close();
}
#endif