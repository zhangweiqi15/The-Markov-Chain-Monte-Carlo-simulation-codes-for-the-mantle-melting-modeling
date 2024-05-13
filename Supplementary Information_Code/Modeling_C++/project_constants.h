#pragma once

#include <string>
#include <iostream>

#include <boost/math/constants/constants.hpp>

using namespace std;

enum MELTING_REGION_SHAPE {Triangle = 0, Curved = 1};
const MELTING_REGION_SHAPE melting_region_shape_global = Triangle;

const int MCMC_length_this_global = 10000;

string output_file_format{ "/mnt/d/local_project/Kane_mantle_array/output/mcmc_x_Kane_mantle" };
    
const double pi = boost::math::constants::pi<double>();
const double h_min = 30;
const double H = 60;
const double theta = pi/4.0;
const double F_max_global = 0.32;

const int var_size_global = 10;
const int y_size_global = 4;

const double x_data[] = {66.7656, 38.9466, 13.3531, 0};
const double data[] = {1.2, 1.5, 2.0, 5.5};
const double sigma[] = {0.05, 0.1, 0.25, 0.1};

const int index_FAD = 0;
const int index_XRP_MAX = 1;
const int index_xrp_start = 2;
const int index_xrp_end = index_xrp_start + 3;
const int index_segments_start = 5;
const int index_segments_end = index_segments_start + 4;
const int index_F_MAX = 9;
// const double var_init_input[] = {
//     0.225, // ancient depletion
//     0.85, // XRP
//     0.13, 0.20,  0.52, // drops of xrp
//     0.5000, 0.3833, 0.0017, 0.1150, // four segments, locations where xrp change, the last segment has 0 RP
//     0.25, //F_max
// };
const double var_init_input[] = {
    0.16, // ancient depletion
    0.85, // XRP
    0.13, 0.20,  0.52, // drops of xrp
    0.5000, 0.3833, 0.0017, 0.1150, // four segments, locations where xrp change, the last segment has 0 RP
    0.32, //F_max
};
// const double var_init_input[] = {
//     0.20, // ancient depletion
//     0.8, // XRP
//     0.33, 0.33,  0.33, // drops of xrp
//     0.2500, 0.2500, 0.2500, 0.2500, // four segments, locations where xrp change, the last segment has 0 RP
//     0.25, // F_max
// };
const double segments_length_total = 60.0;

// const double step_size_input[] = {
//     0.05, // ancient depletion
//     0.05, // XRP
//     0.1, 0.1,  0.1, //relative// drops of xrp
//     0.1, 0.1, 0.1, 0.1, //relative// four segments, locations where xrp change, the last segment has 0 RP
//     0.05, // F_max
// };
const double step_size_input[] = {
    0.05, // ancient depletion
    0.05, // XRP
    0.1, 0.1,  0.1, //relative// drops of xrp
    0.1, 0.1, 0.1, 0.1, //relative// four segments, locations where xrp change, the last segment has 0 RP
    0.02, // F_max
};