
#pragma once

#include <iostream>
#include <math.h>

#include <iomanip>  // setprecision
#include <fstream>

#include <chrono>

#include <algorithm>
#include <numeric>
#include <iostream>
#include <vector>

using namespace std;

void write_table(double** data1, double** data2, double* data3, int length,
	int y_size1, int y_size2, int y_size3, string output_file)
{
	std::cout << std::endl << "Write to the file: ";
	std::cout << output_file << std::endl;

	std::ofstream output_stream;

	output_stream.open(output_file);
	output_stream << length << std::endl;
	output_stream << y_size1 << std::endl;
	output_stream << y_size2 << std::endl;
	output_stream << y_size3 << std::endl;

	for (int ii = 0; ii < length; ii++)
	{
		for (int ip = 0; ip < y_size1; ip++)
			output_stream << data1[ii][ip] << " ";

		for (int ip = 0; ip < y_size2; ip++)
			output_stream << data2[ii][ip] << " ";

		for (int ip = 0; ip < y_size3; ip++)
		{
			//cout << ii << " " << endl;
			output_stream << data3[ii] << " ";
		}

		output_stream << std::endl;
	}
	output_stream.close();
}
