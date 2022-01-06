#ifndef __GRAPH_BUILDER_H
#define __GRAPH_BUILDER_H

#include <math.h>
#include <fstream>
#include <string>
#include <iostream>
#include <string.h>
#include "../include/pricers.h"
#include "../include/greeks.h"

void Price_graph(Input in, std::string output_pde,  const char* type_of_option);

void delta_graph(Input in, std::string output_pde, const char* type_of_option);
#endif