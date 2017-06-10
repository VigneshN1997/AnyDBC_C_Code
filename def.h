#include <stdio.h>
#include <stdlib.h>

enum Status
{
	CORE,
	BORDER,
	NOISE,
	NONE
};

enum Process
{
	PROCESSED,
	UNPROCESSED
};

typedef enum Status Status;

struct object
{
	double* coordinate;
	long long int count_of_neighbours;
	int processed;
	int touched;

};