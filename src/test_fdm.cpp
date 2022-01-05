#include <criterion/criterion.h>
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
        throw e;
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

    srand(time(0));
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

// PRICE INCREMENT TEST

typedef std::vector<double> vec;

vec BS_call_line(vec x_values, double time_to_maturity, double strike, double rate, double volatility) {
    size_t n = x_values.size();
    vec bs_line(n);

    for (size_t i = 0; i < n; i++)
    {
        bs_line.at(i) = BS_call(
            x_values.at(i)
            , strike
            , rate
            , volatility
            , time_to_maturity
        );
    }
    
    return bs_line;
}

struct increment_test {
    double underlying,
    strike,
    rate,
    volatility;
    double dt; double prev_t; vec old_result; vec new_result;
    vec new_result_theoretical;
};

// From a line of true BS prices, find prices at next time step and compare to Black Scholes
increment_test increment_old_result(vec x_values, double dt, double upper_bound) {
    increment_test tst;
    
    tst.prev_t = random_float(upper_bound) + upper_bound / 2,
    tst.strike = random_float(upper_bound) + upper_bound / 2,
    tst.rate = random_float(upper_bound) / 100, 
    tst.volatility = (random_float(upper_bound) + upper_bound / 2) / 100;

    EuropeanCallOption option(tst.strike, tst.rate, tst.volatility);
    BlackScholesPDE pde = BlackScholesPDE(&option);
    
    tst.dt = dt;
    tst.old_result = BS_call_line(
        x_values
        , tst.prev_t
        , tst.strike
        , tst.rate
        , tst.volatility
    );

    tst.new_result = vec(x_values.size());
    store_boundary_conditions(
        tst.prev_t
        , x_values
        , &pde
        , tst.new_result
    );
    tst.new_result = increment_european_price(
        x_values
        , tst.dt
        , tst.prev_t
        , &pde
        , tst.old_result
        , tst.new_result
    );

    tst.new_result_theoretical = BS_call_line(
        x_values
        , tst.prev_t + tst.dt
        , tst.strike
        , tst.rate
        , tst.volatility
    );

    return tst;
}

void print_increment_test(increment_test tst) {
    std::cout << "\nStrike : " << tst.strike;
    std::cout << "\nRate : " << tst.rate;
    std::cout << "\nVolatility : " << tst.volatility;
    std::cout << "\nprev_t : " << tst.prev_t;
    std::cout << "\nnext_t : " << tst.prev_t + tst.dt;

    std::cout << "\nold result : ";
    print_vector(tst.old_result);
    std::cout << "\nnew result : ";
    print_vector(tst.new_result);
    std::cout << "\nnew result theoretical : ";
    print_vector(tst.new_result_theoretical);
    std::cout << '\n';
}

/**
 * @brief checks whether prediciton is close enough to truth
 * 
 * @param predicted the predictions vector
 * @param ground_truth the truth vector
 * @param tol the error tolerance
 * @param max_index for BS PDE solving, for instance, predictions become unusable when the spot gets too big
 * due to limit conditions which are taken at infinity ; hence the error gets too big for big spots and we have to ignore it
 * @return true if for all coordinates (up until max_index), prediction is close enough to ground truth
 * @return false else
 */
bool all_abs_dist_under_tol(vec predicted, vec ground_truth, double tol, size_t max_index) {
    for (size_t i = 0; i < max_index; i++)
    {
        if (abs(predicted.at(i) - ground_truth.at(i)) > tol) {
            return false;
        }
    }
    return true;
}

Test(fdm, increment_european) {

    srand(time(0));
    size_t M = 100, max_index = M / 2;
    double dx = 1, dt = 0.05, upper_bound = 10, tol = 1;
    vec x_values = get_uniform_x_grid(M, dx);

    size_t nb_tests = 10;
    for (size_t i = 0; i < nb_tests; i++)
    {
        increment_test tst = increment_old_result(x_values, dt, upper_bound);
        print_increment_test(tst);
        cr_assert(
            all_abs_dist_under_tol(tst.new_result, tst.new_result_theoretical, tol, max_index)
        );
    }

}