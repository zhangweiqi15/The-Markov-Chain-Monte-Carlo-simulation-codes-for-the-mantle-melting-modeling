#pragma once

#include <string>
#include <iostream>

#include <utilities.h>
#include <data_handling.h>

using namespace std;

const int Gibbs_neighbor_number = 120;

const int use_thread_number = 60;

enum MCMC_type {Gibbs = 0, MH = 1 };

class MCMC
{
public:	
	MCMC(
		const int length,
		const int var_size,
		const int y_size,

		double* var_init,
		double* var_step_size,
		double* var_min,
		double* var_max,

		double* y_data,
		double* y_sigma,

		void (*ypred) (double* var, double* y),
		double (*cost_w) (double* yred, double* y, double* y_sigma, int y_size),
		string output_file_pass
	) :
		length(length),
		var_size(var_size),
		y_size(y_size),

		var_init(var_init),
		var_step_size(var_step_size),
		var_min(var_min),
		var_max(var_max),

		y_data(y_data),
		y_sigma(y_sigma),

		ypred_fun(ypred),
		cost_w(cost_w),
		output_file(output_file_pass)
	{
        mcmc_type = Gibbs;
		var_list = new double* [length];
		pred_list = new double* [length];
		cost_list = new double[length];

		for (int ii = 0; ii < length; ii++)
		{
			var_list[ii] = new double[var_size];
			pred_list[ii] = new double[y_size];
		}

		// string file_head("kane_model");
		// string prefix("mcmc_kane_model");
		// output_file.replace(output_file.find(file_head), file_head.size(), prefix);
	};

    enum MCMC_type mcmc_type;

	const int length;
	const int var_size;
	const int y_size;

	const double* var_init;
	const double* var_step_size;
	const double* var_min;
	const double* var_max;
	double* y_data;
	double* y_sigma;

	void (*ypred_fun) (double* var, double* result);

	double (*cost_w) (double* y, double* yred, double* y_sigma, int y_size);

	double** var_list;
	double** pred_list;
	double* cost_list;

	string output_file;

	int flag_constrain_parameters = 0;
	void Constrain_parameters(int flag, double** neighbor_list);

	void Gibbs_neighbor_generator(double** neighbor_list, double* var_current, int ivar);

	void Neighbor_generator(double** neighbor_list, double* var_current);

	void Batch_Calculator(
		double** var_batch,
		double** pred_batch, // to be calculated
		double* cost_batch // to be calculated
	);

	void Gibbs_sampler(
		double* current_var, // to be updated
		double* current_pred, // to be updated
		double& current_cost, // to be updated
		double** var_neighbor_list,
		double** pred_neighbor_list,
		double* cost_neighbor_list);

	void Gibbs_run();
	void MH_run();
    void run()
    {
        switch (mcmc_type) 
        {
            case Gibbs:
                Gibbs_run();
                break;
            case MH:
                MH_run();
                break;
            default:
                Gibbs_run();
        }
    };

	int display_detail_level = 0;
	void display_state(double* current_var, double* current_pred, double cost)
	{
		cout << " var: ";
		for (int iv = 0; iv < var_size; iv++)
			cout << current_var[iv] << " ";

		cout << " pred: ";
		for (int ip = 0; ip < y_size; ip++)
			cout << current_pred[ip] << " ";

		cout << " cost = " << cost;
		cout << endl;
	};

	void write();

	~MCMC()
	{
		for (int ii = 0; ii < length; ii++)
		{
			delete[] var_list[ii];
			delete[] pred_list[ii];
		}
		delete[] var_list;
		delete[] pred_list;
		delete[] cost_list;
	};
};
void MCMC::Constrain_parameters(int flag, double** neighbor_list)
{
	if (flag == 1)
	{
		#pragma omp parallel for num_threads(use_thread_number) 
        for (int ii = 0; ii < Gibbs_neighbor_number; ii++)
        {
            double tmp = neighbor_list[ii][0];
            if (tmp > neighbor_list[ii][var_size-1])
            {
                neighbor_list[ii][0] = neighbor_list[ii][var_size-1];
                neighbor_list[ii][var_size-1] = tmp;
            }
        }
    }
	if (flag == 2) //  nz-->cc-->aa
	{
	}
}

void MCMC::Gibbs_neighbor_generator(double** neighbor_list, double* var_current, int ivar)
{
	// could be paralleled
	#pragma omp parallel for num_threads(use_thread_number) 
	for (int ii = 0; ii < Gibbs_neighbor_number; ii++)
	{
		for (int iv = 0; iv < var_size; iv++)
			neighbor_list[ii][iv] = var_current[iv];

		if (var_min[ivar] < var_max[ivar] && ii>0)
		{
			double move = (rand_01() - 0.5) * var_step_size[ivar];

			neighbor_list[ii][ivar] += move;

			if (neighbor_list[ii][ivar] < var_min[ivar] || neighbor_list[ii][ivar] > var_max[ivar])
				neighbor_list[ii][ivar] = neighbor_list[ii][ivar] - 2 * move;
		}
	}
}
void MCMC::Neighbor_generator(double** neighbor_list, double* var_current)
{
	// could be paralleled
	#pragma omp parallel for num_threads(use_thread_number) 
	for (int ii = 0; ii < Gibbs_neighbor_number; ii++)
	{
		for (int iv = 0; iv < var_size; iv++)
		{
			neighbor_list[ii][iv] = var_current[iv];

			if (var_min[iv] < var_max[iv] && ii>=0)
			{
				double move = (rand_01() - 0.5) * var_step_size[iv];

				neighbor_list[ii][iv] += move;

				if (neighbor_list[ii][iv] < var_min[iv] || neighbor_list[ii][iv] > var_max[iv])
					neighbor_list[ii][iv] = neighbor_list[ii][iv] - 2 * move;
			}
		}
	}
}
void MCMC::Batch_Calculator(
	double** var_batch,
	double** pred_batch, 
	double* cost_batch // to be calculated
	)
{
	// could be paralleled

	#pragma omp parallel for num_threads(use_thread_number)
	for (int ii = 0; ii < Gibbs_neighbor_number; ii++)
	{
		ypred_fun(var_batch[ii], pred_batch[ii]);
		cost_batch[ii] = cost_w(y_data, pred_batch[ii], y_sigma, y_size);
	}
}

void MCMC::Gibbs_sampler(
	double* current_var,
	double* current_pred,
	double& current_cost,
	double** var_neighbor_list,
	double** pred_neighbor_list,
	double* cost_neighbor_list)
{
	double* p_list = new double[Gibbs_neighbor_number];
	double sum_p = post_p(cost_neighbor_list, p_list, Gibbs_neighbor_number); // update p_list and calculate the sum

	double rand_p = rand_01() * sum_p;
	for (int ii = 0; ii < Gibbs_neighbor_number; ii++)
	{
		rand_p -= p_list[ii];
		if (rand_p <= 0)
		{
			for (int iv = 0; iv < var_size; iv++)
				current_var[iv] = var_neighbor_list[ii][iv];

			for (int ip = 0; ip < y_size; ip++)
				current_pred[ip] = pred_neighbor_list[ii][ip];

			current_cost = cost_neighbor_list[ii];

			break;
		}
	}

	delete[] p_list;
}

void MCMC::Gibbs_run()
{
	double** var_neighbor_list = new double* [Gibbs_neighbor_number];
	double** pred_neighbor_list = new double* [Gibbs_neighbor_number];
	double* cost_neighbor_list = new double[Gibbs_neighbor_number];

	for (int ii = 0; ii < Gibbs_neighbor_number; ii++)
	{
		var_neighbor_list[ii] = new double[var_size];
		pred_neighbor_list[ii] = new double[y_size];
	}

	double* current_var = new double[var_size];
	double* current_pred = new double[y_size];
	double current_cost = 0;

	for (int ii = 0; ii < var_size; ii++)
		current_var[ii] = var_init[ii];

	for (int ii = 0; ii < length; ii++)
	{
		for (int iv = 0; iv < var_size; iv++)
		{
			if (display_detail_level > 0)
				cout << "Gnerate neighbors..." << endl;
			Gibbs_neighbor_generator(
				var_neighbor_list, // to be generated
				current_var, iv);
			if (flag_constrain_parameters>0)
				Constrain_parameters(flag_constrain_parameters, var_neighbor_list);

			if (display_detail_level > 0)
				cout << "Calculate costs of neighbors..." << endl;
			Batch_Calculator(
				var_neighbor_list,
				pred_neighbor_list, // to be calculated
				cost_neighbor_list // to be calculated
			);

			if (display_detail_level > 0)
				cout << "Sampling neighbors..." << endl;
			Gibbs_sampler(
				current_var, // to be updated
				current_pred, // to be updated
				current_cost, // to be updated
				var_neighbor_list,
				pred_neighbor_list,
				cost_neighbor_list);

			if (display_detail_level > 0)
				cout << "Gibbs sampling done for step_" << ii << "-" << iv << endl;

			if (iv == var_size - 1)
			{
				for (int iiv = 0; iiv < var_size; iiv++)
					var_list[ii][iiv] = current_var[iiv];

				for (int iip = 0; iip < y_size; iip++)
					pred_list[ii][iip] = current_pred[iip];

				cost_list[ii] = current_cost;

				printf("%5d ", int(ii));
				for (int iiv = 0; iiv < var_size; iiv++)
				{
					printf("%10.3e ", var_list[ii][iiv]);
				}
				printf("cost = %7.3f\n", cost_list[ii]);

				if (display_detail_level > 0)
				{
					display_state(current_var, current_pred, current_cost);
				}
			}

			if (display_detail_level > 1)
			{
				display_state(current_var, current_pred, current_cost);
			}
		}
	}

	for (int ii = 0; ii < Gibbs_neighbor_number; ii++)
	{
		delete[] var_neighbor_list[ii];
		delete[] pred_neighbor_list[ii];
	}
	delete[] var_neighbor_list;
	delete[] pred_neighbor_list;
	delete[] cost_neighbor_list;

	delete[] current_var;
	delete[] current_pred;
}

void MCMC::MH_run()
{
	double** var_neighbor_list = new double* [Gibbs_neighbor_number];
	double** pred_neighbor_list = new double* [Gibbs_neighbor_number];
	double* cost_neighbor_list = new double[Gibbs_neighbor_number];

	for (int ii = 0; ii < Gibbs_neighbor_number; ii++)
	{
		var_neighbor_list[ii] = new double[var_size];
		pred_neighbor_list[ii] = new double[y_size];
	}

	double* current_var = new double[var_size];
	double* current_pred = new double[y_size];
	double current_cost = 0;

	for (int ii = 0; ii < var_size; ii++)
		current_var[ii] = var_init[ii];

	for (int ii = 0; ii < length; ii++)
	{
		//for (int iv = 0; iv < var_size; iv++)
		{
			if (display_detail_level > 0)
				cout << "Gnerate neighbors..." << endl;
			Neighbor_generator(
				var_neighbor_list, // to be generated
				current_var);

			if (display_detail_level > 0)
				cout << "Calculate costs of neighbors..." << endl;
			Batch_Calculator(
				var_neighbor_list,
				pred_neighbor_list, // to be calculated
				cost_neighbor_list // to be calculated
				);

			if (display_detail_level > 0)
				cout << "Sampling neighbors..." << endl;
			Gibbs_sampler(
				current_var, // to be updated
				current_pred, // to be updated
				current_cost, // to be updated
				var_neighbor_list,
				pred_neighbor_list,
				cost_neighbor_list);

			if (display_detail_level > 0)
				cout << "Gibbs sampling done for step_" << ii << endl;

			//if (iv == var_size - 1)
			{
				for (int iiv = 0; iiv < var_size; iiv++)
					var_list[ii][iiv] = current_var[iiv];

				for (int iip = 0; iip < y_size; iip++)
					pred_list[ii][iip] = current_pred[iip];

				cost_list[ii] = current_cost;

				printf("%5d ", int(ii));
				for (int iiv = 0; iiv < var_size; iiv++)
				{
					printf("%10.3e ", var_list[ii][iiv]);
				}
				printf("cost = %7.3f\n", cost_list[ii]);

				if (display_detail_level > 0)
				{					
					display_state(current_var, current_pred, current_cost);
				}
			}

			if (display_detail_level > 1)
			{
				display_state(current_var, current_pred, current_cost);
			}
		}
	}

	for (int ii = 0; ii < Gibbs_neighbor_number; ii++)
	{
		delete[] var_neighbor_list[ii];
		delete[] pred_neighbor_list[ii];
	}
	delete[] var_neighbor_list;
	delete[] pred_neighbor_list;
	delete[] cost_neighbor_list;

	delete[] current_var;
	delete[] current_pred;
}

void MCMC::write()
{
    string file_head("mcmc_x_");
    string prefix("mcmc_");

    switch (mcmc_type) 
    {
        case Gibbs:
            output_file.replace(output_file.find(file_head), file_head.size(), prefix);
            break;

         case MH:
            {
                string  MH_prefix("mcmc_MH_"); // Declare MH_prefix inside a block
                output_file.replace(output_file.find(file_head), file_head.size(), MH_prefix);
            }
            break;

        default:
            output_file.replace(output_file.find(file_head), file_head.size(), prefix);
            break;
    }

	write_table(var_list, pred_list, cost_list,
		length,
		var_size, y_size, 1,
		output_file);
}