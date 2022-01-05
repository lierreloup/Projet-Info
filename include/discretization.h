#ifndef __DISCRETIZATION_H
#define __DISCRETIZATION_H

#include <vector>

std::vector<double> get_uniform_x_grid(size_t M, double dx);

class Discretization {
    protected :
    double x_dom;      // Spatial extent [0.0, x_dom]
    size_t M;   // Number of spatial differencing points

    Discretization(double _x_dom, double _M) : x_dom(_x_dom), M(_M) {}

    public :
    virtual std::vector<double> get_grid() const = 0;

    double get_x_dom() const { return this->x_dom; }
    double get_M() const{ return this->M; }
};

class UniformDiscretization : public Discretization {
    public :
    UniformDiscretization(double _x_dom, double _M)
        : Discretization(_x_dom, _M) {} 

    std::vector<double> get_grid() const  override;
};

class NonUniformDiscretization : public Discretization {
    double disc_center; // where the non-uniform discretization is most concentrated
    double c; // the smaller c is, the most precise the discretization is around the center (and imprecise elsewhere)

    public :
    NonUniformDiscretization(double _x_dom, double _M, double _disc_center, double _c)
    : Discretization(_x_dom, _M), disc_center(_disc_center), c(_c) {}

    std::vector<double> get_grid() const override;
};



#endif