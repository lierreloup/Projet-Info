#include "../include/criterion/criterion.h"
#include "fdm_static.cpp"
#include <bits/stdc++.h>
#include <math.h>

#include "test_utils.cpp"

//

// LINEAR SYSTEM SOLVING



/**
 * @brief matrix product between a tridiagonal matrix and a vector
 * 
 * @param lower 
 * @param diagonal of length n strictly greater than 1
 * @param upper 
 * @param multiplied_vector 
 * @return std::vector<double> 
 */
std::vector<double> tridiagonal_product(std::vector<double> lower, std::vector<double> diagonal, std::vector<double> upper, std::vector<double> multiplied_vector) {
    try
    {
    size_t n = diagonal.size();
    std::vector<double> result(n);

    result.at(0) = diagonal.at(0) * multiplied_vector.at(0) + upper.at(0) * multiplied_vector.at(1);

    for (size_t i = 1; i < n - 1; i++) {
        result.at(i) = lower.at(i-1) * multiplied_vector.at(i-1) + diagonal.at(i) * multiplied_vector.at(i) + upper.at(i) * multiplied_vector.at(i+1);
    }

    result.at(n-1) = lower.at(n-2) * multiplied_vector.at(n-2) + diagonal.at(n-1) * multiplied_vector.at(n-1);

    return result;
    }


    catch(const std::exception& e)
    {
        std::cerr << e.what() << " in tridiagonal_product" << '\n';
    }
}

/**
 * @brief determinant of a tridiagonal matrix
 * 
 * @param lower 
 * @param diagonal 
 * @param upper 
 * @return double 
 */
double tridiagonal_determinant(std::vector<double> lower, std::vector<double> diagonal, std::vector<double> upper) {

    size_t size = diagonal.size();
    
    double fn, fn_min_1 = abs(diagonal.at(0)), fn_min_2 = 1;
    for (size_t n = 2; n < size; n++)
    {
        double an = diagonal.at(n), bn_min_1 = upper.at(n-1), cn_min_1 = lower.at(n-1);
        fn = an * fn_min_1 - cn_min_1 * bn_min_1 * fn_min_2;
        fn_min_2 = fn_min_1, fn_min_1 = fn; 
    }
    
    return fn_min_1;
}

/**
 * @brief all necessary info for linear system solving test
 * 
 */
struct tri_test {
    std::vector<double> lower,
        diagonal,
        upper,
        x,
        z;
};

tri_test generate_tridiagonal_test(size_t square_matrix_size, double eps = 0.0001, size_t max_iter = 100) {
    try {
    size_t n = square_matrix_size, count = 0;
    std::vector<double> lower, diagonal, upper;

    bool invertible = false;
    while(!invertible && count < max_iter) {
        lower = rand_vals(n-1)
            , diagonal = rand_vals(n)
            , upper = rand_vals(n-1)
        ;

        invertible = abs(tridiagonal_determinant(lower, diagonal, upper)) > eps;
        count++;
    }

    if (count == max_iter) throw "generate_tridiagonal_test : could not generate an invertible matrix";

    std::vector<double> x = rand_vals(n);

    std::vector<double> z = tridiagonal_product(lower, diagonal, upper, x);

    tri_test result;

    result.diagonal = diagonal, result.lower = lower, result.upper = upper;
    result.x = x, result.z = z;

    return result;
    }

    catch(const std::exception& e)
    {
        std::cerr << e.what() << "in generate_tridiagonal_matrix" << '\n';
        throw e;
    }
    
}

double eucdist(std::vector<double> u, std::vector<double> v) {
    double dist = 0;
    size_t n = u.size();

    for (size_t i = 0; i < n; i++)
    {
        double diff = u.at(i) - v.at(i);
        dist += diff * diff;
    }
    
    return sqrt(dist);
}

Test(fdm, tridiagonal_system_solver) {
    size_t max_matrix_size = 10;
    double tol = 0.01; // The tolerated error between actual and found solution

    // Generate random tests for varying matrix size
    for (size_t size = 2; size < max_matrix_size; size++)
    {
        tri_test tst = generate_tridiagonal_test(size);
        std::vector<double> found_x(size);

        Thomas_solver(
            tst.lower
            , tst.diagonal
            , tst.upper
            , tst.z
            , found_x
            , 0
            , size-1
        );

        cr_assert(eucdist(found_x, tst.x) < tol);
    }
}

