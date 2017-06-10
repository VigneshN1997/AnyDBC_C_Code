#include "header.h"

double** readData(char* file_name,long long int* num_records, int* dimension)
{
	FILE* fp = fopen(file_name,"r");
	fscanf(fp,"%lld\n",num_records);
	fscanf(fp,"%d\n",dimension);

	double** dataSet = (double**)malloc((*num_records)*sizeof(double*));
	int record_num,field_num;
	for(record_num = 0; record_num < *num_records; record_num++)
	{
		dataSet[record_num] = (double*)malloc((*dimension)*sizeof(double));
	}

	for(record_num = 0; record_num < *num_records; record_num++)
	{
		for(field_num = 0; field_num < *dimension; field_num++)
		{
			fscanf(fp,"%lf",&dataSet[record_num][field_num]);
		}
	}
	fclose(fp);
	return dataSet;
}

void anyDBC(double** dataSet, int minPts, int epsilon, int alpha, int beta)
{
	
}