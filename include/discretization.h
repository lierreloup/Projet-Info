#ifndef __DISCRETIZATION_H
#define __DISCRETIZATION_H

#include <vector>

// Function to get a uniform grid
std::vector<double> get_uniform_x_grid(size_t M, double dx);

/**
 * @brief 
 * Interface which wraps a particular discretization method and its parameters
 */
class Discretization {
    protected :
    double x_dom;      // Spatial extent [0.0, x_dom]
    size_t M;   // Number of spatial differencing points

    double t_dom; // Temporal extent [0.0, t_dom]
    size_t N; // Number of temporal differencing points

    Discretization(double _x_dom, size_t _M, double _t_dom, size_t _N)
    : x_dom(_x_dom), M(_M), t_dom(_t_dom), N(_N) {}

    public :
    virtual std::vector<double> get_x_grid() const = 0;
    virtual std::vector<double> get_t_grid() const = 0;

    double get_x_dom() const { return this->x_dom; }
    double get_M() const{ return this->M; }
    double get_t_dom() const { return this->t_dom; }
    double get_N() const{ return this->N; }
};

/**
 * @brief 
 * Uniform implementation of the Discretization interface
 */
class UniformDiscretization : public Discretization {
    public :
    UniformDiscretization(double _x_dom, size_t _M, double _t_dom, size_t _N)
    : Discretization(_x_dom, _M, _t_dom, _N) {} 

    std::vector<double> get_x_grid() const  override;
    std::vector<double> get_t_grid() const  override;
};


#endif