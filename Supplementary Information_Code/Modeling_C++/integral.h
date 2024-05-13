#pragma once

#include <math.h>

// boost integration
#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/constants/constants.hpp>

using namespace std;

// method: "gauss", "gauss_kronrod"
// extern const int GK_depth;

// extern const char* integration_method;

const int GK_depth = 15; // for boost integral
const char* integration_method = "gauss_kronrod"; // "gauss_kronrod"; // "gauss";

double integral1d(double l, double r, function<double(double)> func)
{    
    double value = 0;
    {
        using namespace boost::math::quadrature;

        // value = gauss<double, 15>::integrate(func, l, r);
        
        value = gauss_kronrod<double, 15>::integrate(func, l, r, GK_depth, 1e-9);
    }
    return value;

}

double integral2d(double l1, double r1, double l2, double r2, const char* method, function<double(double, double)> func)
{
    double value = 0;

    if (strcmp(method, "gauss") == 0)
    {
        using namespace boost::math::quadrature;

        auto f = [&](double t) {
            auto g = [&](double s) {
                return func(t, s);
            };
            return gauss<double, 7>::integrate(g, l1, r1);
        };
        value = gauss<double, 7>::integrate(f, l2, r2);
    }

    if (strcmp(method, "gauss_kronrod") == 0)
    {
        using namespace boost::math::quadrature;

        auto f = [&](double t) {
            auto g = [&](double s) {
                return func(t, s);
            };
            return gauss_kronrod<double, 15>::integrate(g, l1, r1, GK_depth);
            //return gauss<double, 7>::integrate(g, 0, 2 * pi);
        };
        value = gauss_kronrod<double, 15>::integrate(f, l2, r2, GK_depth, 1e-9);

        //double error = 0;
        //value = gauss_kronrod<double, 15>::integrate(f, l2, r2, GK_depth, 1e-9, &error);
        //cout << "err: " << error << endl;
    }
    return value;
}
