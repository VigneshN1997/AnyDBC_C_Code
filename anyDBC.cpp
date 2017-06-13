#include <iostream>
#include <string>
#include <set>
#include <map>
#include <fstream>
#include <cstdlib>
#include <limits>
using namespace std;


typedef std::numeric_limits< double > db1;

double** readData(char* file_name,long long int* num_records, int* dimension)
{
	ifstream fp;
	fp.open(file_name);
	
	fp >> *num_records;
	fp >> *dimension;
	double** dataSet = (double**)malloc((*num_records)*sizeof(double*));
	long long int record_num,field_num;
	for(record_num = 0; record_num < *num_records; record_num++)
	{
		dataSet[record_num] = (double*)malloc((*dimension)*sizeof(double));
	}

	for(record_num = 0; record_num < *num_records; record_num++)
	{
		for(field_num = 0; field_num < *dimension; field_num++)
		{
			fp >> dataSet[record_num][field_num];
		}
	}
	// fclose(fp);
	return dataSet;
}

int main(int argc, char* argv[])
{
	char* file_name = argv[1];
	long long int num_records = 0; // number of records in file
	int dimension = 0; // number of fields in each record;
	

	double** dataSet = readData(file_name,&num_records,&dimension);
	// int i,j;
	// for(i = 0; i < num_records; i++)
	// {
	// 	for (int j = 0; j < dimension; j++)
	// 	{
	// 		cout.precision(7);
	// 		cout << fixed << dataSet[i][j] << " ";
	// 	}
	// 	cout << "\n";
	// }

	return 0;
}