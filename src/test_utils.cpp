#ifndef __TEST_UTILS_CPP
#define __TEST_UTILS_CPP

#include <bits/stdc++.h>
//#include <math.h>


double random_float(double upper_bound) {
    return static_cast <double> (rand()) / (static_cast <double> (RAND_MAX)) * upper_bound - upper_bound / 2;
}

double randmax() {
    return random_float(1000);
}
/**
 * @brief random vector of size n
 * 
 * @param n 
 * @return std::vector<double> 
 */
std::vector<double> rand_vals(size_t n) {
    std::vector<double> v(n);
    generate(v.begin(), v.end(), randmax);

    return v;
}

#endif