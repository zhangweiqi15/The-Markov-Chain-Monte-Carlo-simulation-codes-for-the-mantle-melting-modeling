#include <iostream>
#include <vector>
#include <string>

#include <omp.h>

#include <utilities.h>
#include <project.h>
using namespace std;


int main() {
    // #pragma omp parallel
    // {
    //     int thread_id = omp_get_thread_num();
    //     // #pragma omp critical
    //     // std::cout << "Hello, World! from Thread " << thread_id << std::endl;
    // }

    Project_Kane project_kane;
    Clock clock;
    clock.start_clock();
    project_kane.run();
    clock.print_time();
    
    double zloc[]= {0,   27.8190,   53.4125,   66.7656};

    // for (int ii = 0; ii < 4; ++ii) 
    // {
    //     double thickness = crustal_thickness(zloc[ii]*tan(theta), const_cast<double*>(var_init_input));
    //     std::cout << "thickness: " << thickness << std::endl;
    // }

    // double thickness = crustal_thickness(0, const_cast<double*>(var_init_input));
    // std::cout << "thickness: " << thickness << std::endl;

    return 0;
}