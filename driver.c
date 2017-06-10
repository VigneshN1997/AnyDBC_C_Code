#include "header.h"

int main(int argc, char *argv[])
{
	char* file_name = argv[1];
	long long int num_records = 0; // number of records in file
	int dimension = 0; // number of fields in each record;
	

	double** dataSet = readData(file_name,&num_records,&dimension);


	return 0;
}