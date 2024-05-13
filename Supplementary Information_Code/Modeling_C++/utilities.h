#pragma once

#include <vector>
#include <algorithm>
#include <cmath>
#include <functional>
#include <integral.h>

#include <random>
#include <chrono>

#include <project_constants.h>

using namespace std;

double interp_fun(const std::vector<double>& a, const std::vector<double>& b, double z_loc) {
    auto it = std::lower_bound(a.begin(), a.end(), z_loc);
    if (it == a.begin()) {
        return b.front();
    }
    if (it == a.end()) {
        return b.back();
    }
    size_t idx = std::distance(a.begin(), it) - 1;
    double t = (z_loc - a[idx]) / (a[idx+1] - a[idx]);
    return b[idx] + t * (b[idx+1] - b[idx]);
}
// example:
// std::vector<double> a = {a1, a2, a3, a4, a5};
// std::vector<double> b = {b1, b2, b3, b4, b5};
// double z_loc = 2.5;

// double result = interp_fun(a, b, z_loc);


std::vector<double> interp_fun(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& z_loc) {
    std::vector<double> result;
    result.reserve(z_loc.size());

    for (double z : z_loc) {
        auto it = std::lower_bound(a.begin(), a.end(), z);
        if (it == a.begin()) {
            result.push_back(b.front());
        } else if (it == a.end()) {
            result.push_back(b.back());
        } else {
            size_t idx = std::distance(a.begin(), it) - 1;
            double t = (z - a[idx]) / (a[idx+1] - a[idx]);
            result.push_back(b[idx] + t * (b[idx+1] - b[idx]));
        }
    }

    return result;
}
// example:
// std::vector<double> a = {a1, a2, a3, a4, a5};
// std::vector<double> b = {b1, b2, b3, b4, b5};
// std::vector<double> z_loc = {z1, z2, z3, z4, z5};

// std::vector<double> interpolated_values = interp_fun(a, b, z_loc);

double width_fun(double z, double H, double theta) {
    double value = 0;
    switch(melting_region_shape_global)
    {
        case Curved:
            value =  z*z / H / tan(theta) * 1.5;
            break;

        case Triangle:
            value =  z/ tan(theta);
            break;

        default:
            value =  z*z / H / tan(theta) * 1.5;
            break;
    }
    return value;
}
// example:
// double z = 3.0;
// double H = 10.0;
// double theta = 0.5;

// double width = width_fun(z, H, theta);

double dFdz_fun(double z, double z_loc, double Fmax, double H, double FAD, std::function<double(double)> X_fun) {
    double threshold = H * (1 - FAD / Fmax);

    if (z >= threshold) {
        double value =X_fun(z_loc + z);
        return Fmax / H * (1 - X_fun(z_loc + z));
    } else {
        return Fmax / H;
    }
}
//example:
// double z = 3.0;
// double z_loc = 1.5;
// double Fmax = 10.0;
// double H = 5.0;
// double FAD = 2.0;

// auto X_fun = [](double z){ return /* define your X_fun calculation here */; };

// double result = dFdz_fun(z, z_loc, Fmax, H, FAD, X_fun);


double rand_01()
{
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<> distribution(0, 1.0);

	//double val = dis(gen);
	//std::cout << val << std::endl;
	//return val;

	return distribution(gen);
}
double rand_normal()
{
	random_device rd;
	mt19937 gen(rd());
	normal_distribution<> distribution{ 0, 1 };

	//double val = dis(gen);
	//std::cout << val << std::endl;
	//return val;

	return distribution(gen);
}
int rand_int_between_ab(int a, int b)
{
	random_device rd;
	mt19937 gen(rd());
	uniform_int_distribution<> distribution(a, b);

	//double val = dis(gen);
	//std::cout << val << std::endl;
	//return val;

	return distribution(gen);
}

const int COST_CUT_OFF = 100;
double post_p(double cost_value)
{
	if (cost_value > COST_CUT_OFF)
		cost_value = COST_CUT_OFF;
	return exp(-cost_value * cost_value);
}

double post_p(double* cost_value_list, double* p_list, int list_length)
{
	double sum_p = 0;
	for (int ii = 0; ii < list_length; ii++)
	{
		p_list[ii] = post_p(cost_value_list[ii]);
		sum_p += p_list[ii];
	}
	return sum_p;
}

double cost_w_fun(double* y, double* pred, double* y_sigma, int N)
{
    double value = 0;
    for (int ii = 0; ii < N; ii++)
    {
        value += (y[ii] - pred[ii]) * (y[ii] - pred[ii]) / (y_sigma[ii] * y_sigma[ii]);
    }
    return value/N;
}

class Clock
{
    public:
    void start_clock()
    {
        start = chrono::high_resolution_clock::now();
    };
    chrono::high_resolution_clock::time_point start;
    chrono::high_resolution_clock::time_point stop;
    void stop_clock()
    {
        stop = chrono::high_resolution_clock::now();
    };
    void print_time()
    {
        stop_clock();
        auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
        cout << "Execution time: " << duration.count() << " ms" << endl;
    };
};