#pragma once

#include <string>
#include <iostream>
#include <algorithm>

#include <boost/math/constants/constants.hpp>

#include <MCMC.h>
#include <project_constants.h>

using namespace std;

double crustal_thickness(double z_loc, double* var)
{
    double FAD = var[index_FAD];
    double XRP_max = var[index_XRP_MAX];
    double F_max = var[index_F_MAX];
    
    std::vector<double> segments(var+index_segments_start, var + index_segments_end);
    double sum = 0;
    for (int i = 0; i < segments.size(); ++i) {
        sum += segments[i];
    }
    for (int i = 0; i < segments.size(); ++i) {
        segments[i] = segments[i]/sum*segments_length_total;
    }
    for (int i = 1; i < segments.size(); ++i) {
        segments[i] = segments[i] + segments[i-1];
    }
    segments.insert(segments.begin(), 0);
    segments.pop_back();
    for (int i = 0; i < segments.size(); ++i) {
        segments[i] = segments[i];
    }

    std::vector<double> xrp_drops(var+index_xrp_start, var + index_xrp_end);
    sum = 0;
    for (int i = 0; i < xrp_drops.size(); ++i) {
        sum += xrp_drops[i];
    }
    for (int i = 0; i < xrp_drops.size(); ++i) {
        xrp_drops[i] = xrp_drops[i]/sum*XRP_max;
    }
    std::vector<double> xrps = {0, 0, 0, 0};
    xrps[0] = XRP_max;
    for (int i = 1; i < xrps.size(); ++i) {
        xrps[i] = xrps[i-1] - xrp_drops[i-1];
    }

    // std::cout<< interp_fun(segments, xrps, 0)<<std::endl;
    // std::cout<< interp_fun(segments, xrps, 60)<<std::endl;
    // std::cout<< interp_fun(segments, xrps, 100)<<std::endl;
    
    double value = 0;
    auto int_fun = [&](double z)
    { 
        auto X_fun =[&] (double z)
        {
            // std::cout<<z<<" "<<interp_fun(segments, xrps, z_loc+z-H)<<std::endl;
            return interp_fun(segments, xrps, z-H);
        };
        // std::cout<<dFdz_fun(z, z_loc, F_max, H, FAD, X_fun)<<" "<<width_fun(z, H, theta)
        // << " "<<z<<" "<<X_fun(z)
        // <<std::endl;
        return dFdz_fun(z, z_loc, F_max, H, FAD, X_fun)*width_fun(z, H, theta);
    };
    value = integral1d(h_min, H, int_fun);
    return value;
}

void kane_crustal_thickness(double* var, double* result)
{
    for (int ii=0; ii<y_size_global; ii++)
    {
        double zloc = x_data[0] - x_data[ii];
        result[ii] = crustal_thickness(zloc*tan(theta), var);
    }
}
class Project_Kane
{
public:
    Project_Kane(
            ):
        var_size(var_size_global),
        y_size(y_size_global),
        output_file(output_file_format)
    {
        // var_size, var_init, var_range, 
        set_var_and_fun();

        // prior_fun
        set_prior_fun();

        // cost_fun
        set_cost_fun();

        // pred_fun
        set_pred_fun();

        set_output_file_name();
    };

    const int var_size;
    const int y_size;

    double* var_init;
    double* var_min;
    double* var_max;
    double* step_size;

    string output_file;

    // 
    void (*prior_fun) (double* var, double** var_range) = nullptr;
    double (*cost_w) (double* y, double* yred, double* y_sigma, int y_size) = nullptr;
    void (*pred_fun) (double* var, double* result) = nullptr;

    void (*calculate_properties_fun) (double* var, double* result);

    void set_var_and_fun()
    {
        var_init = new double[var_size]; 
        std::copy(var_init_input, var_init_input + var_size, var_init);
        
        var_min = new double [var_size]; 
        var_max = new double [var_size]; 

        var_min[0] = 0;
        var_max[0] = F_max_global;
        
        var_min[0] = 0;
        var_max[0] = 0.16;

        var_min[1] = 0;
        var_max[1] = 1;

        for (int i = index_xrp_start; i < index_segments_end; ++i) {
            var_min[i] = 0;
            var_max[i] = 1;
        }
        var_min[index_F_MAX] = 0;
        var_max[index_F_MAX] = F_max_global;

        // var_min[index_F_MAX] = F_max_global;
        // var_max[index_F_MAX] = F_max_global;

        step_size = new double [var_size]; 
        
        std::copy(step_size_input, step_size_input + var_size, step_size);

    }
    void set_prior_fun()
    {
    }
    void set_cost_fun()
    {
        cost_w=&cost_w_fun;
    }
    void set_pred_fun()
    {
        pred_fun = &kane_crustal_thickness;
    }
    void set_output_file_name()
    {
        string file_head("mcmc_x_Kane_mantle");
        string prefix("mcmc_x_Kane_mantle_cv");
        switch(melting_region_shape_global)
        {
            case Curved:
                output_file.replace(output_file.find(file_head), file_head.size(), prefix);
                
                break;

            case Triangle:
                {
                    string alt_prefix("mcmc_x_Kane_mantle_tg");
                    output_file.replace(output_file.find(file_head), file_head.size(), alt_prefix);
                }
                break;

            default:
                output_file.replace(output_file.find(file_head), file_head.size(), prefix);
                break;
        }
        string alt_prefix("mcmc_x_Kane_mantle_hmin_");
        if (h_min>0)
        {
            output_file += "_hmin_" + to_string(h_min);
        }
    }

    void run();

    ~Project_Kane()
    {
        delete[] var_init;
        delete[] var_min;
        delete[] var_max;
        delete[] step_size;
    }
};
void Project_Kane::run()
{
        MCMC mcmc_problem(
        MCMC_length_this_global,
        var_size,
        y_size,

        var_init,
        step_size,
        var_min,
        var_max,

        const_cast<double*> (data),
        const_cast<double*> (sigma),

        pred_fun,
        cost_w,
        output_file);

        mcmc_problem.display_detail_level = 0;
        mcmc_problem.run();
        mcmc_problem.write();
}
