#ifndef __GRAPH_BUILDER_H
#define __GRAPH_BUILDER_H

#include <math.h>
#include <fstream>
#include <string>
#include <iostream>
#include <string.h>
#include "../include/pricers.h"
#include "../include/greeks.h"

void Price_graph(price_inputs in, std::string output_pde, const char* type_of_option);

#endif
