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
set<int> rangeQueryPerformed;
map<int, string> coreList;
map<int, string> borderList;
set<int> noiseList;
map<int, set<int> > coreForPointMap;
set<int> neiNoise;
int alpha = 100;
int minDist = 2;
int minPts = 5;

double** readData(char* file_name,int* num_records, int* dimension);
vector<int> pick(int N, int k);
std::unordered_set<int> pickSet(int N, int k, std::mt19937& gen);
void anyDBC(int num_records, int dimension);
vector<int> getRandomPoints(vector<int> untouchedList);
double distance(double* a,double* b, int dimension);
set<int> performRangeQuery(int point, int num_records,int dimension);

//typedef std::numeric_limits< double > db1;

double** readData(char* file_name,int* num_records, int* dimension)
{
	ifstream fp;
	fp.open(file_name);
	
	fp >> *num_records;
	fp >> *dimension;
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
			fp >> dataSet[record_num][field_num];
		}
	}
	// fclose(fp);
	return dataSet;
}


vector<int> pick(int N, int k) {
    std::random_device rd;
    std::mt19937 gen(rd());

    std::unordered_set<int> elems = pickSet(N, k, gen);

    // ok, now we have a set of k elements. but now
    // it's in a [unknown] deterministic order.
    // so we have to shuffle it:

    vector<int> result(elems.begin(), elems.end());
    std::shuffle(result.begin(), result.end(), gen);
    return result;
}

std::unordered_set<int> pickSet(int N, int k, std::mt19937& gen)
{
    std::unordered_set<int> elems;
    for (int r = N - k; r < N; ++r) {
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


vector<int> getRandomPoints(vector<int> untouchedList)
{
	vector<int> randomPoints;
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
		vector<int> ranList = pick(untouchedList.size(),alpha);
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

set<int> performRangeQuery(int point, int num_records,int dimension)
{
	double* point_a_coordinates = dataSet[point];
	double* point_b_coordinates;
	rangeQueryPerformed.insert(point);
	set<int> neighbourPoints;
	int i;
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

vector<int> assignStateNei(int point, set<int> listOfNeighbors)
{
	vector<int> touchList;
	map<int, string>::iterator itr;			//iterator for borderList
	map<int, string>::iterator itr_core;	//iterator for coreList
	vector<int>::iterator itr_vec;			//iterator for touchList		
	map<int, set<int> >::iterator it_cfpm; 	//iterator for coreForPointMap
	map<int, set<int> >::iterator it_nm;	//iterator for neighbourMap
	set<int>::iterator it_rqp;			//iterator for rangeQueryPerformed
	set<int>::iterator it_nn;			//iterator for neiNoise
	set<int>::iterator itr_nei;				//iterator for listOfNeighbours
	neighbourMap.insert(pair <int, set<int> >(point, listOfNeighbors));
	int i;
	if(listOfNeighbors.size() >= minPts)
	{
		coreList.insert(pair <int, string>(point,"PROCESSED"));
		itr = borderList.find(point);
		if(itr != borderList.end())
		{
			borderList.erase(itr);
		}
		
		for(itr_nei = listOfNeighbors.begin(); itr_nei != listOfNeighbors.end(); itr++)
		{
			int nei = *itr_nei;
			itr_vec = find(touchList.begin(),touchList.end(),nei);
			if(itr_vec == touchList.end())
			{
				touchList.push_back(nei);
			}
			it_cfpm = coreForPointMap.find(nei);
			if(it_cfpm == coreForPointMap.end())
			{
				set<int> s;
				coreForPointMap.insert(pair<int,set<int> >(nei,s));
			}
			coreForPointMap[nei].insert(point);
			it_rqp = find(rangeQueryPerformed.begin(),rangeQueryPerformed.end(),nei);
			if(it_rqp == rangeQueryPerformed.end())
			{
				it_nm = neighbourMap.find(nei);
				if(it_nm == neighbourMap.end())
				{
					set<int> s;
					neighbourMap.insert(pair<int,set<int> >(nei,s));
				}
				neighbourMap[nei].insert(point);
				if(neighbourMap[nei].size() < minPts)
				{
					itr = borderList.find(nei);
					if(itr == borderList.end())
					{
						it_rqp = find(rangeQueryPerformed.begin(),rangeQueryPerformed.end(),nei);
						if(it_rqp != rangeQueryPerformed.end())
						{
							borderList.insert(pair<int,string>(nei,"PROCESSED"));
						}
						else
						{
							borderList.insert(pair<int,string>(nei,"UNPROCESSED"));
						}
					}
				}
				else
				{
					itr_core = coreList.find(nei);
					if(itr_core == coreList.end())
					{
						coreList.insert(pair<int,string>(nei,"UNPROCESSED"));
						itr = borderList.find(nei);
						if(itr != borderList.end())
						{
							borderList.erase(itr);
						}
					}
				}
			}
		}
	}
	else
	{
		noiseList.insert(point);
		it_nn = find(neiNoise.begin(),neiNoise.end(),point);
		if(it_nn != neiNoise.end())
		{
			neiNoise.erase(point);
		}
		for(itr_nei = listOfNeighbors.begin(); itr_nei != listOfNeighbors.end(); itr++)
		{
			int nei = *itr_nei;
			itr_vec = find(touchList.begin(),touchList.end(),nei);
			if(itr_vec == touchList.end())
			{
				touchList.push_back(nei);
			}
			itr_core = coreList.find(nei);
			itr = borderList.find(nei);
			if(itr_core == coreList.end() && itr == borderList.end())
			{
				neiNoise.insert(nei);
			}
			it_nm = neighbourMap.find(nei);
			if(it_nm == neighbourMap.end())
			{
				set<int> s;
				neighbourMap.insert(pair<int,set<int> >(nei,s));
			}
			neighbourMap[nei].insert(point);
		}
	}
	return touchList;
}

void anyDBC(int num_records, int dimension)
{	//STARTING STEP:1
	vector<int> untouchedList;
	vector<int>::iterator itr;
	int i;
	for(i = 0; i < num_records; i++)
	{
		untouchedList.push_back(i);		
	}
	while(untouchedList.size() > 0)
	{
		vector<int> randomPoints = getRandomPoints(untouchedList);
		for(i = 0; i < randomPoints.size(); i++)
		{
			int point = randomPoints.at(i);
			itr = find(untouchedList.begin(),untouchedList.end(),point);
			if(itr != untouchedList.end()) // test this function
			{
				untouchedList.erase(itr);
				set<int> neighbours_of_point = performRangeQuery(point,num_records,dimension);
				vector<int> touch = assignStateNei(point,neighbours_of_point);

			}
		}

	}

}

int main(int argc, char* argv[])
{
	char* file_name = argv[1];
	int num_records = 0; // number of records in file
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