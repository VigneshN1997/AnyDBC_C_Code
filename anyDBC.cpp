#include <iostream>
#include <string>
#include <set>
#include <map>
#include <fstream>
#include <cstdlib>
#include <limits>

#include <vector>
#include <random>
#include <unordered_set>
#include <algorithm>
using namespace std;

map<int, set<int> > neighbourMap; // this will map a point to a set(the set will be the epsilon neighbours of the point)
double** dataSet; // this is a list of all points
vector<long long int> rangeQueryPerformed;
map<int, string> coreList;
int alpha = 100;
int minDist = 2;
int minPts = 5;

double** readData(char* file_name,long long int* num_records, int* dimension);
vector<long long int> pick(long long int N, int k);
std::unordered_set<int> pickSet(long long int N, int k, std::mt19937& gen);
void anyDBC(int num_records, int dimension);
vector<long long int> getRandomPoints(vector<long long int> untouchedList);
double distance(double* a,double* b, int dimension);
set<long long int> performRangeQuery(int point, long long int num_records,int dimension);

//typedef std::numeric_limits< double > db1;

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


vector<long long int> pick(long long int N, int k) {
    std::random_device rd;
    std::mt19937 gen(rd());

    std::unordered_set<int> elems = pickSet(N, k, gen);

    // ok, now we have a set of k elements. but now
    // it's in a [unknown] deterministic order.
    // so we have to shuffle it:

    vector<long long int> result(elems.begin(), elems.end());
    std::shuffle(result.begin(), result.end(), gen);
    return result;
}

std::unordered_set<int> pickSet(long long int N, int k, std::mt19937& gen)
{
    std::unordered_set<int> elems;
    for (long long int r = N - k; r < N; ++r) {
        int v = std::uniform_int_distribution<>(1, r)(gen);

        // there are two cases.
        // v is not in candidates ==> add it
        // v is in candidates ==> well, r is definitely not, because
        // this is the first iteration in the loop that we could've
        // picked something that big.

        if (!elems.insert(v).second) {
            elems.insert(r);
        }   
    }
    return elems;
}


vector<long long int> getRandomPoints(vector<long long int> untouchedList)
{
	vector<long long int> randomPoints;
	int i;
	if(untouchedList.size() <= alpha)
	{
		
		for(i = 0; i < untouchedList.size(); i++)
		{
			randomPoints.push_back(untouchedList.at(i));
		}
	}
	else
	{
		vector<long long int> ranList = pick(untouchedList.size(),alpha);
		for(i = 0; i < ranList.size(); i++)
		{
			randomPoints.push_back(untouchedList.at(ranList.at(i)));
		}
	}
	return randomPoints;
}

double distance(double* a,double* b, int dimension)
{
	double distance = 0;
	int i;
	for(i = 0; i < dimension; i++)
	{
		distance = distance + pow((a[i]-b[i]),2);
		cout << distance << "\n";
	}
	return sqrt(distance);
}

set<long long int> performRangeQuery(int point, long long int num_records,int dimension)
{
	double* point_a_coordinates = dataSet[point];
	double* point_b_coordinates;
	rangeQueryPerformed.push_back(point);
	set<long long int> neighbourPoints;
	long long int i;
	for(i = 0; i < num_records; i++)
	{
		point_b_coordinates = dataSet[i];
		if(point_b_coordinates == point_a_coordinates)
		{
			continue;
		}
		double dist = distance(point_a_coordinates,point_b_coordinates,dimension);
		cout.precision(7);
		cout << fixed <<point_a_coordinates[0]<<" "<< point_b_coordinates[0] << " " <<dist << "\n";
		if(dist <= minDist)
		{
			neighbourPoints.insert(i);
		}
	}
	return neighbourPoints;
}

vector<long long int> assignStateNei(long long int point, vector<long long int> listOfNeighbors)
{
	vector<long long int> touchList;
	neighbourMap[point] = listOfNeighbors;
	if(listOfNeighbors.size() >= minPts)
	{

	}
}

void anyDBC(int num_records, int dimension)
{	//STARTING STEP:1
	vector<long long int> untouchedList;
	vector<long long int>::iterator itr;
	long long int i;
	for(i = 0; i < num_records; i++)
	{
		untouchedList.push_back(i);		
	}
	while(untouchedList.size() > 0)
	{
		vector<long long int> randomPoints = getRandomPoints(untouchedList);
		for(i = 0; i < randomPoints.size(); i++)
		{
			long long int point = randomPoints.at(i);
			itr = find(untouchedList.begin(),untouchedList.end(),point);
			if(itr != untouchedList.end()) // test this function
			{
				untouchedList.erase(itr);
				set<long long int> neighbours_of_point = performRangeQuery(point,num_records,dimension);
				touch = assignStateNei(point,neighbours_of_point);

			}
		}

	}

}

int main(int argc, char* argv[])
{
	char* file_name = argv[1];
	long long int num_records = 0; // number of records in file
	int dimension = 0; // number of fields in each record;
	

	dataSet = readData(file_name,&num_records,&dimension);
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
	anyDBC(num_records,dimension);
	return 0;
}