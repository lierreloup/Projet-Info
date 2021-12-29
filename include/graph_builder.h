#ifndef __GRAPH_BUILDER_H
#define __GRAPH_BUILDER_H

#include <math.h>
#include <fstream>
#include <string>
#include <iostream>
#include <string.h>

void Price_graph(double spot, double time_to_maturity, double strike, double rate, double volatility, std::string output_pde, const char* type_of_option);

#endif
