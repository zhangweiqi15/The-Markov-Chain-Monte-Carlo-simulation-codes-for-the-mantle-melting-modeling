Modelling crustal production at Kane fracture zone
Author: Boda Liu
Date: 2023/12/25
Email: bodaliu@mail.iggcas.ac.cn

1.	Get ready
1) 	Make sure you have a C++ compiler. We use g++ in Windows Subsystem for Linux. There are many tutorials, e.g. https://code.visualstudio.com/docs/cpp/config-wsl
2) 	We need to install boost library (https://www.boost.org/). Please follow the steps in the official documentation. The version utilized in this study is 1.71.0.0ubuntu2.
3) 	We need to enable openmp parallel computation.  You can add flag to the compiler. In visual studio, we add one line “-fopenmp” in task.json
 

2.	Input data, set up the model
1) 	To change input data, modify lines 26-27 in project_constants.h
 
where x_data is the distance away from the axis, data is the crustal thickness in km.

2) 	To change the starting and final melting depth change lines 18-19 in project_constants.h
 

3.	Specify parameters in Monte Carlo sampling
1) 	To change the number of Gibbs steps, modify line 13 in project_constants.h
 

2) 	To change the number of models in each Gibbs step, modify line 11 in MCMC.h
 

3) 	To change the number of parallel threads, modify line 13 in MCMC.h
 

4.	Set up output file
Specify the output folder. Here, I save the output in disk D of the Windows system. Please keep the default filename, “mcmc_x_Kane_mantle”.
 

5.	Run
Compile the main cpp file (hellowworld.cpp) and run. You can use command line. We just hit the run button in visual studio. Outputs are model parameters and cost in each line. The execution time for a 64-core workstation is typically a couple of hours. Of cause the time depends on the computation power, Gibbs steps, and the number of models to search in each step.
 
