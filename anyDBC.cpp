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
map<int, set<int> > clusters;
map<int, set<int> > edgeNo;
map<int, set<int> > edgeYes;
map<int, set<int> > edgeWeak;
map<int, set<int> > edgeUnknown;
map<int, int> visitedNode;


map<int, int> usizeList;
map<int, double> statList;
map<int, double> degList;


int alpha = 100;
int beta = 50;
int minDist = 2;
int minPts = 5;

double** readData(char* file_name,int* num_records, int* dimension);
vector<int> pick(int N, int k);
std::unordered_set<int> pickSet(int N, int k, std::mt19937& gen);
void anyDBC(int num_records, int dimension);
vector<int> getRandomPoints(vector<int> untouchedList);
double distance(double* a,double* b, int dimension);
set<int> performRangeQuery(int point, int num_records,int dimension);
vector<int> assignStateNei(int point, set<int> listOfNeighbors);
void createPCIR(int point);
int ddcBetPCIR(int core1,int core2,int dimension);
void connComp();
void DFS(int u, int rep);
void calculateStatDegree(int num_records);
vector<int> calculateScore(int num_records);
bool stoppingCondition();
void processOutliers();
//typedef std::numeric_limits< double > db1;


template<typename A, typename B>
std::pair<B,A> flip_pair(const std::pair<A,B> &p)
{
    return std::pair<B,A>(p.second, p.first);
}

template<typename A, typename B>
std::multimap<B,A> flip_map(const std::map<A,B> &src)
{
    std::multimap<B,A> dst;
    std::transform(src.begin(), src.end(), std::inserter(dst, dst.begin()), 
                   flip_pair<A,B>);
    return dst;
}





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
		
		for(itr_nei = listOfNeighbors.begin(); itr_nei != listOfNeighbors.end(); itr_nei++)
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
		for(itr_nei = listOfNeighbors.begin(); itr_nei != listOfNeighbors.end(); itr_nei++)
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


void createPCIR(int point)
{
	set<int> s;
	clusters.insert(pair<int,set<int> >(point,s));
	clusters[point].insert(point);
}

int ddcBetPCIR(int core1,int core2,int dimension)
{
	double dist = distance(dataSet[core1],dataSet[core2],dimension);
	if(dist > 3*minDist)
	{
		return 1;
	}
	set<int> pointPCIR = neighbourMap[core1];
	pointPCIR.insert(core1);
	set<int> neiPointPCIR = neighbourMap[core2];
	neiPointPCIR.insert(core2);
	set<int> intersectionPoints;
	set_intersection(pointPCIR.begin(),pointPCIR.end(),neiPointPCIR.begin(),neiPointPCIR.end(),std::inserter(intersectionPoints,intersectionPoints.begin()));
	set<int> coreListKeys;
	map<int, string>::iterator map_itr;
	for(map_itr = coreList.begin(); map_itr != coreList.end(); map_itr++)
	{
		coreListKeys.insert(map_itr->first);
	}
	set<int> coreIntersectionPoints;
	set_intersection(intersectionPoints.begin(),intersectionPoints.end(),coreListKeys.begin(),coreListKeys.end(),std::inserter(coreIntersectionPoints,coreIntersectionPoints.begin()));
	int iCore = coreIntersectionPoints.size();
	if((iCore > 0) || ((dist > 1.732*minDist)&&(intersectionPoints.size() >= minPts)))
	{
		return 0;
	}
	else if(intersectionPoints.size() > 0)
	{
		return 2;
	}
	else
	{
		return 3;
	}
}

void connComp()
{
	set<int> repComp;
	visitedNode.clear();
	set<int> edgeYesKeys;
	map<int,set<int> >::iterator edgeYesItr;
	map<int,int>::iterator visited_itr;
	for(edgeYesItr = edgeYes.begin(); edgeYesItr != edgeYes.end(); edgeYesItr++)
	{
		edgeYesKeys.insert(edgeYesItr->first);
	}
	set<int>::iterator keys_itr;
	for(keys_itr = edgeYesKeys.begin(); keys_itr != edgeYesKeys.end(); keys_itr++)
	{
		int u = *keys_itr;
		visited_itr = visitedNode.find(u);
		if(visited_itr == visitedNode.end())
		{
			DFS(u,u);
			repComp.insert(u);
		}
	}
	set<int> deleteClust;
	set<int>::iterator del_itr;
	map<int, set<int> >::iterator clus_itr;
	set_difference(edgeYesKeys.begin(),edgeYesKeys.end(),repComp.begin(),repComp.end(),std::inserter(deleteClust,deleteClust.begin()));
	for(del_itr = deleteClust.begin(); del_itr != deleteClust.end(); del_itr++)
	{
		int u = *del_itr;
		clus_itr = clusters.find(u);
		if(clus_itr != clusters.end())
		{
			clusters.erase(clus_itr);
		}
	}
	edgeYes.clear();
	set<int>::iterator rep_itr;
	set<int>::iterator no_itr;
	set<int>::iterator weak_itr;
	map<int,set<int> >::iterator edgeNoItr;
	map<int,set<int> >::iterator edgeWeakItr;

	for(visited_itr = visitedNode.begin(); visited_itr != visitedNode.end(); visited_itr++)
	{
		int u = visited_itr->first;
		int rep = visited_itr->second;
		for(rep_itr = clusters[rep].begin(); rep_itr != clusters[rep].end(); rep_itr++)
		{
			int point = *rep_itr;
			if(point == u)
			{
				continue;
			}
			edgeNoItr = edgeNo.find(u);
			if(edgeNoItr != edgeNo.end())
			{
				no_itr = find(edgeNo[u].begin(),edgeNo[u].end(),point);
				if(no_itr != edgeNo[u].end())
				{
					edgeNo[u].erase(no_itr);
				}
			}
			edgeWeakItr = edgeWeak.find(u);
			if(edgeWeakItr != edgeWeak.end())
			{
				weak_itr = find(edgeWeak[u].begin(),edgeWeak[u].end(),point);
				if(weak_itr != edgeWeak[u].end())
				{
					edgeWeak[u].erase(weak_itr);
				}
			}
		}
		edgeNoItr = edgeNo.find(u);
		if(edgeNoItr != edgeNo.end())
		{
			if(edgeNo[u].size() == 0)
			{
				edgeNo.erase(edgeNoItr);
			}
		}
		edgeWeakItr = edgeWeak.find(u);
		if(edgeWeakItr != edgeWeak.end())
		{
			if(edgeWeak[u].size() == 0)
			{
				edgeWeak.erase(edgeWeakItr);
			}
		}
	}
}

void DFS(int u, int rep)
{
	visitedNode.insert(pair<int, int>(u,rep));
	set<int>::iterator clus_itr;
	for(clus_itr = clusters[u].begin(); clus_itr != clusters[u].end(); clus_itr++)
	{
		int v = *clus_itr;
		clusters[rep].insert(v);
	}
	set<int>::iterator yes_itr;
	map<int,int>::iterator visited_itr;
	for(yes_itr = edgeYes[u].begin(); yes_itr != edgeYes[u].end(); yes_itr++)
	{
		int v = *yes_itr;
		visited_itr = visitedNode.find(v);
		if(visited_itr == visitedNode.end())
		{
			DFS(v,rep);
		}
	}
}

void calculateStatDegree(int num_records)
{
	int w = clusters.size();
	map<int, int> numBorderPoints;
	map<int, set<int> >::iterator clus_itr;
	vector<int> alreadyCountedPoint;
	set<int>::iterator v_itr;		// v iterator
	set<int>::iterator nei_itr;		//neighbourMap[p] iterator
	vector<int>::iterator al_itr; 	// alreadyCountedPoint iterator
	set<int>::iterator range_itr;	//rangeQueryPerformed iterator
	map<int, string>::iterator border_itr;//borderList iterator
	for(clus_itr = clusters.begin(); clus_itr != clusters.end(); clus_itr++)
	{
		int k = clus_itr->first;
		set<int> v = clus_itr->second;
		usizeList.insert(pair<int,int>(k,0));
		statList.insert(pair<int,double>(k,0));
		numBorderPoints.insert(pair<int,int>(k,0));
		int noOfpointsInCluster = 0;
		alreadyCountedPoint.clear();
		for(v_itr = v.begin(); v_itr != v.end(); v_itr++)
		{
			int p = *v_itr;
			for(nei_itr = neighbourMap[p].begin(); nei_itr!= neighbourMap[p].end(); nei_itr++)
			{
				int x = *nei_itr;
				al_itr = find(alreadyCountedPoint.begin(),alreadyCountedPoint.end(),x);
				if(al_itr == alreadyCountedPoint.end())
				{
					noOfpointsInCluster++;
					range_itr = find(rangeQueryPerformed.begin(),rangeQueryPerformed.end(),x);
					if(range_itr == rangeQueryPerformed.end())
					{
						usizeList[k]++;
					}
					border_itr = borderList.find(x);
					if(border_itr != borderList.end())
					{
						numBorderPoints[k]++;
					}
					alreadyCountedPoint.push_back(x);
				}
			}
		}
		statList[k] = ((double)usizeList[k]/(double)noOfpointsInCluster) + ((double)noOfpointsInCluster/(double)num_records);
	}
	map<int, set<int> >::iterator edgeWeakItr;
	map<int, set<int> >::iterator edgeWeakItr2;
	map<int, set<int> >::iterator edgeUnknownItr;
	map<int, set<int> >::iterator edgeUnknownItr2;
	set<int>::iterator weak_itr;
	set<int>::iterator un_itr;
	map<int, double>::iterator stat_itr; //statList iterator
	for(clus_itr = clusters.begin(); clus_itr != clusters.end(); clus_itr++)
	{
		int u = clus_itr->first;
		int siValue = 0;
		degList.insert(pair<int,double>(u,0));
		edgeWeakItr = edgeWeak.find(u);
		if(edgeWeakItr != edgeWeak.end())
		{
			for(weak_itr = edgeWeak[u].begin(); weak_itr != edgeWeak[u].end(); weak_itr++)
			{
				int v = *weak_itr;
				stat_itr = statList.find(v);
				if(stat_itr != statList.end())
				{
					degList[u]+=statList[v];
					siValue++;
				}
				else
				{
					edgeWeakItr2 = edgeWeak.find(v);
					if(edgeWeakItr2 != edgeWeak.end())
					{
						edgeWeak.erase(edgeWeakItr2);
					}
				}
			}
			degList[u] *= w;
		}
		edgeUnknownItr = edgeUnknown.find(u);
		if(edgeUnknownItr != edgeUnknown.end())
		{
			for(un_itr = edgeUnknown[u].begin(); un_itr != edgeUnknown[u].end(); un_itr++)
			{
				int v = *un_itr;
				stat_itr = statList.find(v);
				if(stat_itr != statList.end())
				{
					degList[u]+=statList[v];
					siValue++;
				}
				else
				{
					edgeUnknownItr2 = edgeUnknown.find(v);
					if(edgeUnknownItr2 != edgeUnknown.end())
					{
						edgeUnknown.erase(edgeUnknownItr2);
					}
				}
			}
		}
		if(numBorderPoints[u] == 0)
		{
			siValue = 0;
		}
		degList[u] -= siValue;
	}
}

vector<int> calculateScore(int num_records)
{
	map<int, double> scoreSet;
	set<int> unprocessedPoints;
	set<int> range;
	set<int>::iterator unp_itr;
	map<int, string>::iterator core_itr;
	map<int, string>::iterator border_itr;
	map<int, set<int> >::iterator clus_itr;
	for(int i = 0; i < num_records; i++)
	{
		range.insert(i);
	}
	set_difference(range.begin(),range.end(),rangeQueryPerformed.begin(),rangeQueryPerformed.end(),std::inserter(unprocessedPoints,unprocessedPoints.begin()));
	set<int> unprocessedPoints1;
	set_difference(unprocessedPoints.begin(),unprocessedPoints.end(),neiNoise.begin(),neiNoise.end(),std::inserter(unprocessedPoints1,unprocessedPoints1.begin()));
	unprocessedPoints = unprocessedPoints1;
	for(unp_itr = unprocessedPoints.begin(); unp_itr != unprocessedPoints.end(); unp_itr++)
	{
		int i = *unp_itr;
		core_itr = coreList.find(i);
		border_itr = borderList.find(i);
		if(core_itr != coreList.end() || border_itr != borderList.end())
		{
			double score = 0;
			for(clus_itr = clusters.begin(); clus_itr != clusters.end(); clus_itr++)
			{
				int rep = clus_itr->first;
				set<int> coreListRep = clus_itr->second;
				set<int> intersectCore;
				set_intersection(coreForPointMap[i].begin(),coreForPointMap[i].end(),coreListRep.begin(),coreListRep.end(),std::inserter(intersectCore,intersectCore.begin()));
				if(intersectCore.size() > 0)
				{
					score += degList[rep];
				}
			}
			score += ((double)1)/((double)neighbourMap[i].size());
			scoreSet.insert(pair<int,double>(i,score));
		}
	}

   	multimap<double, int> sorted_x = flip_map(scoreSet);
   	multimap<double, int>::reverse_iterator sorted_itr;
   	vector<int> returnList;
   	int betaF = beta;
   	if(sorted_x.size() > 0)
   	{
   	    if(sorted_x.size() < betaF)
   	    {
   	        betaF = sorted_x.size();
   	    }
   	    sorted_itr=sorted_x.rbegin();
   	    for(int i = 0; i < betaF; i++)
   	    {
   	        returnList.push_back(sorted_itr->second);
   	        sorted_itr++;
   	    }
   	}
   	return returnList;

}

bool stoppingCondition()
{
	return false;
}

void processOutliers()
{
	
}


void anyDBC(int num_records, int dimension)
{	//STARTING STEP:1
	vector<int> untouchedList;
	vector<int>::iterator itr; 				//iterator for untouchedList
	map<int,string>::iterator itr_core; 	//iterator for coreList
	map<int,set<int> >::iterator itr_clus1;	//iterator for clusters
	map<int,set<int> >::iterator itr_clus2;	//iterator for clusters
	set<int>::iterator itr_done;			//iterator for donePoints

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
				itr_core = coreList.find(point);
				if(itr_core != coreList.end())
				{
					createPCIR(point);
				}
				int j;
				for(j = 0; j < touch.size(); j++)
				{
					int p = touch[i];
					itr = find(untouchedList.begin(),untouchedList.end(),p);
					if(itr != untouchedList.end())
					{
						untouchedList.erase(itr);
					}
				}
			}
		}
	}
	/*######################################################################*/
	cout << "Initial count of range queries: "<< rangeQueryPerformed.size() << "\n";
	set<int> donePoints;
	int point;
	map<int,set<int> >::iterator edgeNoItr;
	map<int,set<int> >::iterator edgeYesItr;
	map<int,set<int> >::iterator edgeWeakItr;
	map<int,set<int> >::iterator edgeUnknownItr;
	for(itr_clus1 = clusters.begin(); itr_clus1 != clusters.end(); itr_clus1++)
	{
		point = itr_clus1->first;
		donePoints.insert(point);
		int neiPoint;
		for(itr_clus2 = clusters.begin(); itr_clus2 != clusters.end(); itr_clus2++)
		{
			neiPoint = itr_clus2->first;
			itr_done = find(donePoints.begin(),donePoints.end(),neiPoint);
			if(itr_done != donePoints.end())
			{
				continue;
			}
			int stat = ddcBetPCIR(point,neiPoint,dimension);
			if(stat == 1)
			{
				edgeNoItr = edgeNo.find(point);
				if(edgeNoItr == edgeNo.end())
				{
					set<int> s;
					edgeNo.insert(pair<int,set<int> >(point,s));
				}
				edgeNo[point].insert(neiPoint);
				edgeNoItr = edgeNo.find(neiPoint);
				if(edgeNoItr == edgeNo.end())
				{
					set<int> s;
					edgeNo.insert(pair<int,set<int> >(neiPoint,s));
				}
				edgeNo[neiPoint].insert(point); 
			}
			else if(stat == 0)
			{
				edgeYesItr = edgeYes.find(point);
				if(edgeYesItr == edgeYes.end())
				{
					set<int> s;
					edgeYes.insert(pair<int,set<int> >(point,s));
				}
				edgeYes[point].insert(neiPoint);
				edgeYesItr = edgeYes.find(neiPoint);
				if(edgeYesItr == edgeYes.end())
				{
					set<int> s;
					edgeYes.insert(pair<int,set<int> >(neiPoint,s));
				}
				edgeYes[neiPoint].insert(point); 
			}
			else if(stat == 2)
			{
				edgeWeakItr = edgeWeak.find(point);
				if(edgeWeakItr == edgeWeak.end())
				{
					set<int> s;
					edgeWeak.insert(pair<int,set<int> >(point,s));
				}
				edgeWeak[point].insert(neiPoint);
				edgeWeakItr = edgeWeak.find(neiPoint);
				if(edgeWeakItr == edgeWeak.end())
				{
					set<int> s;
					edgeWeak.insert(pair<int,set<int> >(neiPoint,s));
				}
				edgeWeak[neiPoint].insert(point); 
			}
			else if(stat == 3)
			{
				edgeUnknownItr = edgeUnknown.find(point);
				if(edgeUnknownItr == edgeUnknown.end())
				{
					set<int> s;
					edgeUnknown.insert(pair<int,set<int> >(point,s));
				}
				edgeUnknown[point].insert(neiPoint);
				edgeUnknownItr = edgeUnknown.find(neiPoint);
				if(edgeUnknownItr == edgeUnknown.end())
				{
					set<int> s;
					edgeUnknown.insert(pair<int,set<int> >(neiPoint,s));
				}
				edgeUnknown[neiPoint].insert(point); 				
			}
		}
	}
	/*######################################################################*/
	if(edgeYes.size() > 0)
	{
		connComp();
		//mergeAssignNewEdge(); //implement it
	}
	int iteration = 0;
	while(true)
	{
		cout << "Range queries : " << rangeQueryPerformed.size() << "\n";
		cout << "Iteration : " << iteration << "    No. of clusters : " << clusters.size() << "\n";
       	cout << "CoreList : " << coreList.size() << "\n";
       	cout << "BorderList : " << borderList.size() << "\n";
       	cout << "NoiseList : " << noiseList.size() + neiNoise.size();
		if(stoppingCondition()) //implement it
		{
			/*######################################################################*/
			calculateStatDegree(num_records);
			vector<int> betaPoints = calculateScore(num_records);
			if(betaPoints.size() == 0)
			{
				break;
			}
			/*######################################################################*/
			vector<int>::iterator beta_itr;
			for(beta_itr = betaPoints.begin(); beta_itr != betaPoints.end(); beta_itr++)
			{

			}
		}
		else
		{
			break;
		}
		iteration++;
	}
	processOutliers();

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