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
#include <ctime>
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
map<int, set<int> > edgeWeakBckup;
map<int, set<int> > edgeNoBckup;
map<int, int> visitedNode;


map<int, int> usizeList;
map<int, double> statList;
map<int, double> degList;

clock_t start_time = clock();

vector<int> popFromNoise;
map<std::pair<int, int>, int> noStatusLastIteration;

int dimension = 0;
int num_records = 0;
int alpha = 512;
int beta = 512;
int minDist = 2;
int minPts = 5;

double** readData(char* file_name);
vector<int> pick(int N, int k);
std::unordered_set<int> pickSet(int N, int k, std::mt19937& gen);
void anyDBC();
vector<int> getRandomPoints(vector<int> untouchedList);
double distance(double* a,double* b);
set<int> performRangeQuery(int point);
vector<int> assignStateNei(int point, set<int> listOfNeighbors);
void createPCIR(int point);
int ddcBetPCIR(int core1,int core2);
void connComp();
void DFS(int u, int rep);
void calculateStatDegree();
vector<int> calculateScore();
void dccBetPCLU(int repClust1, int repClust2);
void updateStates();
bool stoppingCondition();
void processNoise(int p);
void processOutliers();
void mergeAssignNewEdge();
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

double** readData(char* file_name)
{
	ifstream fp;
	fp.open(file_name);
	
	fp >> num_records;
	fp >> dimension;
	dataSet = (double**)malloc((num_records)*sizeof(double*));
	int record_num,field_num;
	for(record_num = 0; record_num < num_records; record_num++)
	{
		dataSet[record_num] = (double*)malloc((dimension)*sizeof(double));
	}

	for(record_num = 0; record_num < num_records; record_num++)
	{
		for(field_num = 0; field_num < dimension; field_num++)
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

double distance(double* a,double* b)
{
	double distance = 0;
	int i;
	for(i = 0; i < dimension; i++)
	{
		distance = distance + pow((a[i]-b[i]),2);
	}
	return sqrt(distance);
}

set<int> performRangeQuery(int point)
{
	double* point_a_coordinates = dataSet[point];
	double* point_b_coordinates;
	rangeQueryPerformed.insert(point);
	set<int> neighbourPoints;
	int i;
	for(i = 0; i < num_records; i++)
	{
		point_b_coordinates = dataSet[i];
		if(point_b_coordinates == point_a_coordinates) // point==i
		{
			continue;
		}
		double dist = distance(point_a_coordinates,point_b_coordinates);
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
	neighbourMap[point] = listOfNeighbors;
	int i;
	if(listOfNeighbors.size() >= minPts)
	{
		coreList[point] = "PROCESSED";
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
				coreForPointMap[nei] = s;
			}
			coreForPointMap[nei].insert(point);
			it_rqp = find(rangeQueryPerformed.begin(),rangeQueryPerformed.end(),nei);
			if(it_rqp == rangeQueryPerformed.end())
			{
				it_nm = neighbourMap.find(nei);
				if(it_nm == neighbourMap.end())
				{
					set<int> s;
					neighbourMap[nei] = s;
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
							borderList[nei] = "PROCESSED";
						}
						else
						{
							borderList[nei] = "UNPROCESSED";
						}
					}
				}
				else
				{
					itr_core = coreList.find(nei);
					if(itr_core == coreList.end())
					{
						coreList[nei] = "UNPROCESSED";
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
			neiNoise.erase(it_nn);
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
			if((itr_core == coreList.end()) && (itr == borderList.end()))
			{
				neiNoise.insert(nei);
			}
			it_nm = neighbourMap.find(nei);
			if(it_nm == neighbourMap.end())
			{
				set<int> s;
				neighbourMap[nei] = s;
			}
			neighbourMap[nei].insert(point);
		}
	}
	return touchList;
}


void createPCIR(int point)
{
	set<int> s;
	clusters[point] = s;
	clusters[point].insert(point);
}

int ddcBetPCIR(int core1,int core2)
{
	double dist = distance(dataSet[core1],dataSet[core2]);
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
	visitedNode[u] = rep;
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

void calculateStatDegree()
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
		usizeList[k] = 0;
		statList[k] = 0;
		numBorderPoints[k] = 0;
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
		degList[u] = 0;
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

vector<int> calculateScore()
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
			scoreSet[i] = score;
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

void dccBetPCLU(int repClust1, int repClust2)
{
	int noOfSubClusters = clusters[repClust1].size()*clusters[repClust2].size();
	int noCount = 0;
	bool formedYesEdge = false;
	bool formedWeakEdge = false;
	set<int>::iterator rep1_itr;
	set<int>::iterator rep2_itr;
	map<int,set<int> >::iterator edgeNoItr;
	map<int,set<int> >::iterator edgeYesItr;
	map<int,set<int> >::iterator edgeWeakItr;
	map<int,set<int> >::iterator edgeUnknownItr;
	set<int>::iterator weak_itr;
	for(rep1_itr = clusters[repClust1].begin(); rep1_itr != clusters[repClust1].end(); rep1_itr++)
	{
		int core1 = *rep1_itr;
		for(rep2_itr = clusters[repClust2].begin(); rep2_itr != clusters[repClust2].end(); rep2_itr++)
		{
			int core2 = *rep2_itr;
			if(core1 == core2)
			{
				continue;
			}
			int stat = ddcBetPCIR(core1,core2);
			if(stat == 0)
			{
				formedYesEdge = true;
				edgeYesItr = edgeYes.find(repClust1);
				if(edgeYesItr == edgeYes.end())
				{
					set<int> s;
					edgeYes[repClust1] = s;
				}
				edgeYes[repClust1].insert(repClust2);
				edgeYesItr = edgeYes.find(repClust2);
				if(edgeYesItr == edgeYes.end())
				{
					set<int> s;
					edgeYes[repClust2] = s;
				}
				edgeYes[repClust2].insert(repClust1);
				if(formedWeakEdge)
				{
					weak_itr = find(edgeWeak[repClust1].begin(),edgeWeak[repClust1].end(),repClust2);
					if(weak_itr != edgeWeak[repClust1].end())
					{
						edgeWeak[repClust1].erase(weak_itr);
					}
					if(edgeWeak[repClust1].size() == 0)
					{
						edgeWeakItr = edgeWeak.find(repClust1);
						if(edgeWeakItr != edgeWeak.end())
						{
							edgeWeak.erase(edgeWeakItr);
						}
					}
					weak_itr = find(edgeWeak[repClust2].begin(),edgeWeak[repClust2].end(),repClust1);
					if(weak_itr != edgeWeak[repClust2].end())
					{
						edgeWeak[repClust2].erase(weak_itr);
					}
					if(edgeWeak[repClust2].size() == 0)
					{
						edgeWeakItr = edgeWeak.find(repClust2);
						if(edgeWeakItr != edgeWeak.end())
						{
							edgeWeak.erase(edgeWeakItr);
						}
					}
				}
				return;
			}
			else if(stat == 1)
			{
				noCount++;
			}
			else if(stat == 2 && !(formedYesEdge))
			{
				formedWeakEdge = true;
				edgeWeakItr = edgeWeak.find(repClust1);
				if(edgeWeakItr == edgeWeak.end())
				{
					set<int> s;
					edgeWeak[repClust1] = s;
				}
				edgeWeak[repClust1].insert(repClust2);
				edgeWeakItr = edgeWeak.find(repClust2);
				if(edgeWeakItr == edgeWeak.end())
				{
					set<int> s;
					edgeWeak[repClust2] = s;
				}
				edgeWeak[repClust2].insert(repClust1);
			}
		}
	}
	if(!(formedWeakEdge) && !(formedYesEdge) && (noCount == noOfSubClusters))
	{
		edgeNoItr = edgeNo.find(repClust1);
		if(edgeNoItr == edgeNo.end())
		{
			set<int> s;
			edgeNo[repClust1] = s;
		}
		edgeNo[repClust1].insert(repClust2);
		edgeNoItr = edgeNo.find(repClust2);
		if(edgeNoItr == edgeNo.end())
		{
			set<int> s;
			edgeNo[repClust2] = s;
		}
		edgeNo[repClust2].insert(repClust1);
	}
	else if(!(formedWeakEdge) && !(formedYesEdge) && (noCount != noOfSubClusters))
	{
		edgeUnknownItr = edgeUnknown.find(repClust1);
		if(edgeUnknownItr == edgeUnknown.end())
		{
			set<int> s;
			edgeUnknown[repClust1] = s;
		}
		edgeUnknown[repClust1].insert(repClust2);
		edgeUnknownItr = edgeUnknown.find(repClust2);
		if(edgeUnknownItr == edgeUnknown.end())
		{
			set<int> s;
			edgeUnknown[repClust2] = s;
		}
		edgeUnknown[repClust2].insert(repClust1);	
	}
	return;
}

void updateStates()
{
	vector<int> popWeakUnknown;
	map<int, set<int> >::iterator clus_itr;
	set<int>::iterator v_itr;
	set<int>::iterator nei_itr;
	set<int>::iterator range_itr;
	map<int,int>::iterator usize_itr;
	map<int,set<int> >::iterator edgeNoItr;
	map<int,set<int> >::iterator edgeWeakItr;
	map<int,set<int> >::iterator edgeUnknownItr;
	set<int>::iterator node_itr;
	set<int>::iterator weak_itr;
	set<int>::iterator unknown_itr;
	vector<int>::iterator pwu_itr;
	for(clus_itr = clusters.begin(); clus_itr != clusters.end(); clus_itr++)
	{
		int k = clus_itr->first;
		set<int> v = clus_itr->second;
		usizeList[k] = 0;
		for(v_itr = v.begin(); v_itr != v.end(); v_itr++)
		{
			int p = *v_itr;
			for(nei_itr = neighbourMap[p].begin(); nei_itr != neighbourMap[p].end(); nei_itr++)
			{
				int x = *nei_itr;
				range_itr = find(rangeQueryPerformed.begin(),rangeQueryPerformed.end(),x);
				if(range_itr == rangeQueryPerformed.end())
				{
					usizeList[k]++;
				}
			}
		}
		usize_itr = usizeList.find(k);
		if(usize_itr != usizeList.end())
		{
			if(usizeList[k] == 0)
			{
				for(edgeWeakItr = edgeWeak.begin(); edgeWeakItr != edgeWeak.end(); edgeWeakItr++)
				{
					int node = edgeWeakItr->first;
					set<int> neiNodes = edgeWeakItr->second;
					node_itr = find(neiNodes.begin(),neiNodes.end(),k);
					if(k == node)
					{
						for(node_itr = neiNodes.begin(); node_itr != neiNodes.end(); node_itr++)
						{
							int knode = *node_itr;
							edgeNoItr = edgeNo.find(k);
							if(edgeNoItr == edgeNo.end())
							{
								set<int> s;
								edgeNo[k] = s;
							}
							edgeNo[k].insert(knode);
							edgeNoItr = edgeNo.find(knode);
							if(edgeNoItr == edgeNo.end())
							{
								set<int> s;
								edgeNo[knode] = s;
							}
							edgeNo[knode].insert(k);
						}
					}
					else if(node_itr != neiNodes.end())
					{
						weak_itr = find(edgeWeak[node].begin(),edgeWeak[node].end(),k);
						if(weak_itr != edgeWeak[node].end())
						{
							edgeWeak[node].erase(weak_itr);
						}
						if(edgeWeak[node].size() == 0)
						{
							popWeakUnknown.push_back(node);
						}
					}
				}
				for(edgeUnknownItr = edgeUnknown.begin(); edgeUnknownItr != edgeUnknown.end(); edgeUnknownItr++)
				{
					int node = edgeUnknownItr->first;
					set<int> neiNodes = edgeUnknownItr->second;
					node_itr = find(neiNodes.begin(),neiNodes.end(),k);
					if(k == node)
					{
						for(node_itr = neiNodes.begin(); node_itr != neiNodes.end(); node_itr++)
						{
							int knode = *node_itr;
							edgeNoItr = edgeNo.find(k);
							if(edgeNoItr == edgeNo.end())
							{
								set<int> s;
								edgeNo[k] = s;
							}
							edgeNo[k].insert(knode);
							edgeNoItr = edgeNo.find(knode);
							if(edgeNoItr == edgeNo.end())
							{
								set<int> s;
								edgeNo[knode] = s;
							}
							edgeNo[knode].insert(k);
						}
					}
					else if(node_itr != neiNodes.end())
					{
						unknown_itr = find(edgeUnknown[node].begin(),edgeUnknown[node].end(),k);
						if(unknown_itr != edgeUnknown[node].end())
						{
							edgeUnknown[node].erase(unknown_itr);
						}
						if(edgeUnknown[node].size() == 0)
						{
							popWeakUnknown.push_back(node);
						}
					}
				}
				edgeWeakItr = edgeWeak.find(k);
				if(edgeWeakItr != edgeWeak.end())
				{
					edgeWeak.erase(edgeWeakItr);
				}
				edgeUnknownItr = edgeUnknown.find(k);
				if(edgeUnknownItr != edgeUnknown.end())
				{
					edgeUnknown.erase(edgeUnknownItr);
				}
			}
		}
	}
	for(pwu_itr = popWeakUnknown.begin(); pwu_itr != popWeakUnknown.end(); pwu_itr++)
	{
		int node = *pwu_itr;
		edgeWeakItr = edgeWeak.find(node);
		if(edgeWeakItr != edgeWeak.end())
		{
			edgeWeak.erase(edgeWeakItr);
		}
		edgeUnknownItr = edgeUnknown.find(node);
		if(edgeUnknownItr != edgeUnknown.end())
		{
			edgeUnknown.erase(edgeUnknownItr);
		}
	}
}

void processNoise(int p)
{

	map<int, string>::iterator core_itr;
	map<int, string>::iterator border_itr;
	core_itr = coreList.find(p);
	border_itr = borderList.find(p);
	set<int> coreListKeys;
	map<int, string>::iterator map_itr;
	for(map_itr = coreList.begin(); map_itr != coreList.end(); map_itr++)
	{
		coreListKeys.insert(map_itr->first);
	}
	set<int> listOfNeighbors = neighbourMap[p];
	set<int> intersect;
	set_intersection(listOfNeighbors.begin(),listOfNeighbors.end(),coreListKeys.begin(),coreListKeys.end(),std::inserter(intersect,intersect.begin()));
	
	if((core_itr != coreList.end()) || (border_itr != borderList.end()))
	{
		popFromNoise.push_back(p);
	}
	else if(intersect.size() > 0)
	{
		borderList[p] = "PROCESSED";
		popFromNoise.push_back(p);
	}
	else
	{
		set<int>::iterator range_itr;
		set<int>::iterator nei_itr;
		set<int>::iterator lon_itr;
		map<int,set<int> >::iterator nm_itr;
		for(nei_itr = neighbourMap[p].begin(); nei_itr != neighbourMap[p].end(); nei_itr++)
		{
			int nei = *nei_itr;
			range_itr = find(rangeQueryPerformed.begin(),rangeQueryPerformed.end(),nei);
			if(range_itr != rangeQueryPerformed.end())
			{
				listOfNeighbors = neighbourMap[nei];
			}
			else
			{
				listOfNeighbors = performRangeQuery(nei);
				neighbourMap[nei] = listOfNeighbors;
			}
			for(lon_itr = listOfNeighbors.begin(); lon_itr != listOfNeighbors.end(); lon_itr++)
			{
				int neiN = *lon_itr;
				range_itr = find(rangeQueryPerformed.begin(),rangeQueryPerformed.end(),neiN);
				if(range_itr == rangeQueryPerformed.end())
				{
					nm_itr = neighbourMap.find(neiN);
					if(nm_itr == neighbourMap.end())
					{
						set<int> s;
						neighbourMap[neiN] = s;
					}
					neighbourMap[neiN].insert(nei);
					if(neighbourMap[neiN].size() >= minPts)
					{
						coreList[neiN] = "UNPROCESSED";
						border_itr = borderList.find(neiN);
						if(border_itr != borderList.end())
						{
							borderList.erase(border_itr);
						}
					}
				}
			}
			if(listOfNeighbors.size() >= minPts)
			{
				coreList[nei] = "PROCESSED";
				borderList[p] = "PROCESSED";
				popFromNoise.push_back(nei);
				popFromNoise.push_back(p);
				border_itr = borderList.find(nei);
				if(border_itr != borderList.end())
				{
					borderList.erase(border_itr);
				}
				createPCIR(nei);
				map<int,set<int> >::iterator clus_itr;
				for(clus_itr = clusters.begin(); clus_itr != clusters.end(); clus_itr++)
				{
					int repCore = clus_itr->first;
					if(nei != repCore)
					{
						dccBetPCLU(repCore, nei);
					}
				}
				break;
			}
		}
	}
}

void processOutliers()
{
	set<int>::iterator noise_itr;
	vector<int>::iterator pfn_itr; //popFromNoise iterator
	set<int>::iterator nn_itr; //neiNoise iterator
	map<int, string>::iterator core_itr;
	map<int, string>::iterator border_itr;
	set<int>::iterator range_itr;
	map<int, set<int> >::iterator clus_itr;
	map<int, set<int> >::iterator nm_itr; //neighbourMap iterator
	for(noise_itr = noiseList.begin(); noise_itr != noiseList.end(); noise_itr++)
	{
		int p = *noise_itr;
		processNoise(p);
	}
	for(pfn_itr = popFromNoise.begin(); pfn_itr != popFromNoise.end(); pfn_itr++)
	{
		int p = *pfn_itr;
		noise_itr = find(noiseList.begin(),noiseList.end(),p);
		if(noise_itr != noiseList.end())
		{
			noiseList.erase(noise_itr);
		}
		nn_itr = find(neiNoise.begin(),neiNoise.end(),p);
		if(nn_itr != neiNoise.end())
		{
			neiNoise.erase(nn_itr);
		}
	}
	popFromNoise.clear();
	set<int> listNei;
	set<int>::iterator ln_itr; // listNei iterator
	for(nn_itr = neiNoise.begin(); nn_itr != neiNoise.end(); nn_itr++)
	{
		int p = *nn_itr;
		core_itr = coreList.find(p);
		border_itr = borderList.find(p);
		if((core_itr != coreList.end()) || (border_itr != borderList.end()))
		{
			popFromNoise.push_back(p);
			continue;
		}
		range_itr = find(rangeQueryPerformed.begin(),rangeQueryPerformed.end(),p);
		if(range_itr != rangeQueryPerformed.end())
		{
			listNei = neighbourMap[p];
		}
		else
		{
			listNei = performRangeQuery(p);
			neighbourMap[p] = listNei;
		}
		bool isCore = false;
		if(listNei.size() >= minPts)
		{
			isCore = true;
			coreList[p] = "PROCESSED";
			border_itr = borderList.find(p);
			if(border_itr != borderList.end())
			{
				borderList.erase(border_itr);
			}
			popFromNoise.push_back(p);
			createPCIR(p);
			for(clus_itr = clusters.begin(); clus_itr != clusters.end(); clus_itr++)
			{
				int repCore = clus_itr->first;
				if(p != repCore)
				{
					dccBetPCLU(repCore, p);
				}
			}
		}
		for(ln_itr = listNei.begin(); ln_itr != listNei.end(); ln_itr++)
		{
			int nei = *ln_itr;
			range_itr = find(rangeQueryPerformed.begin(),rangeQueryPerformed.end(),nei);
			if(range_itr == rangeQueryPerformed.end())
			{
				nm_itr = neighbourMap.find(nei);
				if(nm_itr == neighbourMap.end())
				{
					set<int> s;
					neighbourMap[nei] = s;
				}
				neighbourMap[nei].insert(p);
				if(neighbourMap[nei].size() >= minPts)
				{
					coreList[nei] = "UNPROCESSED";
					border_itr = borderList.find(nei);
					if(border_itr != borderList.end())
					{
						borderList.erase(border_itr);
					}
				}
			}
			if(isCore)
			{
				borderList[nei] = "UNPROCESSED";
			}
		}
		processNoise(p);
	}
	for(pfn_itr = popFromNoise.begin(); pfn_itr != popFromNoise.end(); pfn_itr++)
	{
		int p = *pfn_itr;
		noise_itr = find(noiseList.begin(),noiseList.end(),p);
		if(noise_itr != noiseList.end())
		{
			noiseList.erase(noise_itr);
		}
		nn_itr = find(neiNoise.begin(),neiNoise.end(),p);
		if(nn_itr != neiNoise.end())
		{
			neiNoise.erase(nn_itr);
		}
	}
	popFromNoise.clear();
}

void mergeAssignNewEdge()
{
	edgeUnknown.clear();
	set<int> eWN;
	set<int> edgeWeakBckupKeys;
	set<int> edgeNoBckupKeys;
	map<int, set<int> >::iterator ewbk_itr;
	map<int, set<int> >::iterator enbk_itr;
	for(ewbk_itr = edgeWeakBckup.begin(); ewbk_itr != edgeWeakBckup.end(); ewbk_itr++)
	{
		edgeWeakBckupKeys.insert(ewbk_itr->first);
	}
	for(enbk_itr = edgeNoBckup.begin(); enbk_itr != edgeNoBckup.end(); enbk_itr++)
	{
		edgeNoBckupKeys.insert(enbk_itr->first);
	}
	set_union(edgeWeakBckupKeys.begin(),edgeWeakBckupKeys.end(),edgeNoBckupKeys.begin(),edgeNoBckupKeys.end(),std::inserter(eWN,eWN.begin()));
	set<int> doneClusters;
	map<int,set<int> >::iterator itr_clus1;
	map<int,set<int> >::iterator itr_clus2;
	
	for(itr_clus1 = clusters.begin(); itr_clus1 != clusters.end(); itr_clus1++)
	{
        int k = itr_clus1->first;
		set<int> v1_1;
		set<int> v1;
		v1_1 = clusters[k];
        set_difference(v1_1.begin(),v1_1.end(),eWN.begin(),eWN.end(),std::inserter(v1,v1.begin()));
        doneClusters.insert(k);
        for(itr_clus2 = clusters.begin(); itr_clus2 != clusters.end(); itr_clus2++)
		{
			int k2 = itr_clus2->first;
            set<int>::iterator doneClusters_itr;
			doneClusters_itr = find(doneClusters.begin(),doneClusters.end(),k2);
			if(doneClusters_itr != doneClusters.end())
			{
					continue;
			}
            set<int> v2_1;
			set<int> v2;
			v2_1 = clusters[k2];
			set_difference(v2_1.begin(),v2_1.end(),eWN.begin(),eWN.end(),std::inserter(v2,v2.begin()));
			int weakPresent = 0;
			int noCount = 0;
			int noOfSubClusters = (v1_1.size())*(v2.size());
			set<int>::iterator itr_v1_1;
			set<int>::iterator itr_v2;
			for(itr_v1_1 = v1_1.begin(); itr_v1_1 != v1_1.end(); itr_v1_1++)
			{
				int point = *itr_v1_1;
				for(itr_v2 = v2.begin(); itr_v2 != v2.end(); itr_v2++)
				{
					int neipoint = *itr_v2;
					map<int, set<int> >::iterator itr_point,itr_pointe;
					itr_point = edgeWeak.find(point);
					itr_pointe = edgeWeakBckup.find(point);
					if(itr_point != edgeWeak.end())
					{
						set<int>::iterator itr_neipoint;
						itr_neipoint = find(edgeWeak[point].begin(),edgeWeak[point].end(),neipoint);
						if(itr_neipoint != edgeWeak[point].end())
						{
							weakPresent = 1;
							edgeWeak[point].erase(itr_neipoint);
							set<int>::iterator itr_point2;
							itr_point2 = find(edgeWeak[neipoint].begin(),edgeWeak[neipoint].end(),point);
							if(itr_point2 != edgeWeak[neipoint].end())
							{
								edgeWeak[neipoint].erase(itr_point2);
							}
							if(edgeWeak[point].size()==0)
							{
								edgeWeak.erase(itr_point);
							}
							if(edgeWeak[neipoint].size()==0)
							{
								map<int, set<int> >::iterator itr_neipoint2;
								itr_neipoint2 =edgeWeak.find(neipoint);
								if(itr_neipoint2 != edgeWeak.end())
								{
									edgeWeak.erase(itr_neipoint2);
								}
							}
							map<int, set<int> >::iterator itr_point3;
							itr_point3 = edgeWeakBckup.find(point);
							if(itr_point3==edgeWeakBckup.end())
							{
								set<int> s;
								edgeWeakBckup[point] = s;
							}
							edgeWeakBckup[point].insert(neipoint);
							
							map<int, set<int> >::iterator itr_neipoint3;
							itr_neipoint3 = edgeWeakBckup.find(neipoint);
							if(itr_neipoint3==edgeWeakBckup.end())
							{
								set<int> s;
								edgeWeakBckup[neipoint] = s;
							}
							edgeWeakBckup[neipoint].insert(point);	
						}
					}
					else if( itr_pointe!=edgeWeakBckup.end())
					{
						set<int>::iterator itr_neipoint;
						itr_neipoint = find(edgeWeakBckup[point].begin(),edgeWeakBckup[point].end(),neipoint);
						if(itr_neipoint != edgeWeakBckup[point].end())
						{
							weakPresent = 1;
						}	
					}
					itr_point = edgeNo.find(point);
					itr_pointe = edgeNoBckup.find(point);	
					if(itr_point!=edgeNo.end())
					{
						set<int>::iterator itr_neipoint;
						itr_neipoint = find(edgeNo[point].begin(),edgeNo[point].end(),neipoint);
						if(itr_neipoint != edgeNo[point].end())
						{
							noCount+=1;
							edgeNo[point].erase(itr_neipoint);
							set<int>::iterator itr_point2;
							itr_point2 = find(edgeNo[neipoint].begin(),edgeNo[neipoint].end(),point);
							if(itr_point2 != edgeNo[neipoint].end())
							{
								edgeNo[neipoint].erase(itr_point2);
							}
							if(edgeNo[point].size()==0)
							{
								edgeNo.erase(itr_point);
							}
							if(edgeNo[neipoint].size()==0)
							{
								map<int, set<int> >::iterator itr_neipoint2;
								itr_neipoint2 = edgeNo.find(neipoint);
								if(itr_neipoint2 != edgeNo.end())
								{
									edgeNo.erase(itr_neipoint2);
								}
							}
							map<int, set<int> >::iterator itr_point3;
							itr_point3 = edgeNoBckup.find(point);
							if(itr_point3==edgeNoBckup.end())
							{
								set<int> s;
								edgeNoBckup[point] = s;
							}
							edgeNoBckup[point].insert(neipoint);
							
							map<int, set<int> >::iterator itr_neipoint3;
							itr_neipoint3 = edgeNoBckup.find(neipoint);
							if(itr_neipoint3==edgeNoBckup.end())
							{
								set<int> s;
								edgeNoBckup[neipoint] = s;
							}
							edgeNoBckup[neipoint].insert(point);							
						}
					}				
					else if(itr_pointe!=edgeNoBckup.end())
					{
						set<int>::iterator itr_neipoint;
						itr_neipoint = find(edgeNoBckup[point].begin(),edgeNoBckup[point].end(),neipoint);
						if(itr_neipoint != edgeNoBckup[point].end())
						{
							noCount+=1;
						}
					}
				}	
			}
			if(weakPresent==1)
			{
				map<int, set<int> >::iterator edgeWeakItr;
				edgeWeakItr = edgeWeak.find(k);
				if(edgeWeakItr == edgeWeak.end())
				{
					set<int> s;
					edgeWeak[k] = s;
				}
				edgeWeak[k].insert(k2);
				
				map<int, set<int> >::iterator edgeWeakItr2;
				edgeWeakItr2 = edgeWeak.find(k2);
				if(edgeWeakItr2 == edgeWeak.end())
				{
					set<int> s;
					edgeWeak[k2] = s;
				}
				edgeWeak[k2].insert(k);
			}
			else if(noCount==noOfSubClusters)
			{
				if(v2==v2_1 && v1==v1_1)
				{
					map<int, set<int> >::iterator edgeNoItr;
					edgeNoItr = edgeNo.find(k);
					if(edgeNoItr == edgeNo.end())
					{
						set<int> s;
						edgeNo[k] = s;
					}
					edgeNo[k].insert(k2);
				
					map<int, set<int> >::iterator edgeNoItr2;
					edgeNoItr2 = edgeNo.find(k2);
					if(edgeNoItr2 == edgeNo.end())
					{
						set<int> s;
						edgeNo[k2] = s;
					}
					edgeNo[k2].insert(k);
					std::pair<int, int> var;
					var.first = k;
					var.second = k2;
					noStatusLastIteration[var] = noCount;
				}
				else
				{
					int noCountNew = 0;
					set<int> old2;
					set<int> old1;
					set_difference(v2_1.begin(),v2_1.end(),v2.begin(),v2.end(),std::inserter(old2,old2.begin()));
					set_difference(v1_1.begin(),v1_1.end(),v1.begin(),v1.end(),std::inserter(old1,old1.begin()));
					set<int>::iterator itr_old2;
					set<int>::iterator itr_v1;
					for(itr_old2 = old2.begin(); itr_old2 != old2.end(); itr_old2++)
					{
						int point = *itr_old2;
						for(itr_v1 = v1.begin(); itr_v1 != v1.end(); itr_v1++)
						{
							int neipoint = *itr_v1;
							map<int, set<int> >::iterator itr_point;
							map<int, set<int> >::iterator itr_pointe;
							itr_point = edgeNo.find(point);
							itr_pointe = edgeNoBckup.find(point);
							if(itr_point != edgeNo.end())
							{
								set<int>::iterator itr_neipoint;
								itr_neipoint = find(edgeNo[point].begin(),edgeNo[point].end(),neipoint);
								if(itr_neipoint != edgeNo[point].end())
								{
									noCountNew+=1;
									edgeNo[point].erase(itr_neipoint);
									set<int>::iterator itr_point2;
									itr_point2 = find(edgeNo[neipoint].begin(),edgeNo[neipoint].end(),point);
									if(itr_point2 != edgeNo[neipoint].end())
									{
										edgeNo[neipoint].erase(itr_point2);
									}
									if(edgeNo[point].size()==0)
									{
										edgeNo.erase(itr_point);
									}
									if(edgeNo[neipoint].size()==0)
									{
										map<int, set<int> >::iterator itr_neipoint2;
										itr_neipoint2 = edgeNo.find(neipoint);
										if(itr_neipoint2 != edgeNo.end())
										{
											edgeNo.erase(itr_neipoint2);
										}
									}
									map<int, set<int> >::iterator itr_point3;
									itr_point3 = edgeNoBckup.find(point);
									if(itr_point3==edgeNoBckup.end())
									{
										set<int> s;
										edgeNoBckup[point] = s;
									}
									edgeNoBckup[point].insert(neipoint);
									
									map<int, set<int> >::iterator itr_neipoint3;
									itr_neipoint3 = edgeNoBckup.find(neipoint);
									if(itr_neipoint3==edgeNoBckup.end())
									{
										set<int> s;
										edgeNoBckup[neipoint] = s;
									}
									edgeNoBckup[neipoint].insert(point);
								}
							}				
							else if(itr_pointe!=edgeNoBckup.end())
							{
								set<int>::iterator itr_neipoint;
								itr_neipoint = find(edgeNoBckup[point].begin(),edgeNoBckup[point].end(),neipoint);
								if(itr_neipoint != edgeNoBckup[point].end())
								{
									noCountNew+=1;
								}
							}
							set<int>::iterator itr_old1;
							for(itr_old1 = old1.begin(); itr_old1 != old1.end(); itr_old1++)
							{
								int p = *itr_old1;
								pair<int,int> pair1,pair2;
								pair1.first = point;
								pair1.second = p;
								pair2.first = p;
								pair2.second = point;
								auto nslt_itr1 = noStatusLastIteration.find(pair1);
								auto nslt_itr2 = noStatusLastIteration.find(pair2);
								if(nslt_itr1 != noStatusLastIteration.end())
								{
									noCount = noCount + noStatusLastIteration[pair1];
									noStatusLastIteration.erase(nslt_itr1);
								}
								else if(nslt_itr2 != noStatusLastIteration.end())
								{
									noCount = noCount + noStatusLastIteration[pair2];
									noStatusLastIteration.erase(nslt_itr2);
								}
							}
						}
					}
					if(noCountNew==(v1.size()*old2.size()))
					{
						int totalClust = noCount + noCountNew;
						if(totalClust == (v1_1.size()*v2_1.size()))
						{
							map<int, set<int> >::iterator edgeNoItr;
							edgeNoItr = edgeNo.find(k);
							if(edgeNoItr == edgeNo.end())
							{
								set<int> s;
								edgeNo[k] = s;
							}
							edgeNo[k].insert(k2);
						
							map<int, set<int> >::iterator edgeNoItr2;
							edgeNoItr2 = edgeNo.find(k2);
							if(edgeNoItr2 == edgeNo.end())
							{
								set<int> s;
								edgeNo[k2] = s;
							}
							edgeNo[k2].insert(k);
							std::pair<int, int> var;
							var.first = k;
							var.second = k2;
							noStatusLastIteration[var] = totalClust;
						}
					}
				}
			}
			else
			{
				map<int, set<int> >::iterator edgeUnknownItr;
				edgeUnknownItr = edgeUnknown.find(k);
				if(edgeUnknownItr == edgeUnknown.end())
				{
					set<int> s;
					edgeUnknown[k] = s;
				}
				edgeUnknown[k].insert(k2);
			
				map<int, set<int> >::iterator edgeUnknownItr2;
				edgeUnknownItr2 = edgeUnknown.find(k2);
				if(edgeUnknownItr2 == edgeUnknown.end())
				{
					set<int> s;
					edgeUnknown[k2] = s;
				}
				edgeUnknown[k2].insert(k);
			}
		}
	}
	set<int> removeKeys;
	set<int>::iterator remove_itr;
	map<int, set<int> >::iterator edgeWeak_itr;
	set<int> edgeWeakKeys,edgeUnknownKeys,edgeNoKeys,clustersKeys;
	map<int, set<int> >::iterator map_itr;
	for(map_itr = edgeWeak.begin(); map_itr != edgeWeak.end(); map_itr++)
	{
		edgeWeakKeys.insert(map_itr->first);
	}
	for(map_itr = edgeNo.begin(); map_itr != edgeNo.end(); map_itr++)
	{
		edgeNoKeys.insert(map_itr->first);
	}
	for(map_itr = edgeUnknown.begin(); map_itr != edgeUnknown.end(); map_itr++)
	{
		edgeUnknownKeys.insert(map_itr->first);
	}
	for(map_itr = clusters.begin(); map_itr != clusters.end(); map_itr++)
	{
		clustersKeys.insert(map_itr->first);
	}
	set_difference(edgeWeakKeys.begin(),edgeWeakKeys.end(),clustersKeys.begin(),clustersKeys.end(),std::inserter(removeKeys,removeKeys.begin()));
	for(remove_itr = removeKeys.begin(); remove_itr != removeKeys.end(); remove_itr++)
	{
		int u = *remove_itr;
		edgeWeak_itr = edgeWeak.find(u);
		if(edgeWeak_itr != edgeWeak.end())
		{
			edgeWeak.erase(edgeWeak_itr);
		}
	}

	map<int, set<int> >::iterator edgeNo_itr;
	set_difference(edgeNoKeys.begin(),edgeNoKeys.end(),clustersKeys.begin(),clustersKeys.end(),std::inserter(removeKeys,removeKeys.begin()));
	for(remove_itr = removeKeys.begin(); remove_itr != removeKeys.end(); remove_itr++)
	{
		int u = *remove_itr;
		edgeNo_itr = edgeNo.find(u);
		if(edgeNo_itr != edgeNo.end())
		{
			edgeNo.erase(edgeNo_itr);
		}
	}
	
	map<int, set<int> >::iterator edgeUnknown_itr;
	set_difference(edgeUnknownKeys.begin(),edgeUnknownKeys.end(),clustersKeys.begin(),clustersKeys.end(),std::inserter(removeKeys,removeKeys.begin()));
	for(remove_itr = removeKeys.begin(); remove_itr != removeKeys.end(); remove_itr++)
	{
		int u = *remove_itr;
		edgeUnknown_itr = edgeUnknown.find(u);
		if(edgeUnknown_itr != edgeUnknown.end())
		{
			edgeUnknown.erase(edgeUnknown_itr);
		}
	}
}

bool stoppingCondition()
{
	if(edgeWeak.size()>0||edgeUnknown.size()>0)
	{
		return true;
	}
	else
	{
		return false;
	}
}

void anyDBC()
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
				set<int> neighbours_of_point = performRangeQuery(point);
				vector<int> touch = assignStateNei(point,neighbours_of_point);
				itr_core = coreList.find(point);
				if(itr_core != coreList.end())
				{
					createPCIR(point);
				}
				int j;
				for(j = 0; j < touch.size(); j++)
				{
					int p = touch[j];
					itr = find(untouchedList.begin(),untouchedList.end(),p);
					if(itr != untouchedList.end())
					{
						untouchedList.erase(itr);
					}
				}
			}
		}
	}
	cout << "Step1 rangeQuery on intial alpha points---" << clock() - start_time << "seconds ---";
	start_time = clock();
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
			int stat = ddcBetPCIR(point,neiPoint);
			if(stat == 1)
			{
				edgeNoItr = edgeNo.find(point);
				if(edgeNoItr == edgeNo.end())
				{
					set<int> s;
					edgeNo[point] = s;
				}
				edgeNo[point].insert(neiPoint);
				edgeNoItr = edgeNo.find(neiPoint);
				if(edgeNoItr == edgeNo.end())
				{
					set<int> s;
					edgeNo[neiPoint] = s;
				}
				edgeNo[neiPoint].insert(point); 
			}
			else if(stat == 0)
			{
				edgeYesItr = edgeYes.find(point);
				if(edgeYesItr == edgeYes.end())
				{
					set<int> s;
					edgeYes[point] = s;
				}
				edgeYes[point].insert(neiPoint);
				edgeYesItr = edgeYes.find(neiPoint);
				if(edgeYesItr == edgeYes.end())
				{
					set<int> s;
					edgeYes[neiPoint] = s;
				}
				edgeYes[neiPoint].insert(point); 
			}
			else if(stat == 2)
			{
				edgeWeakItr = edgeWeak.find(point);
				if(edgeWeakItr == edgeWeak.end())
				{
					set<int> s;
					edgeWeak[point] = s;
				}
				edgeWeak[point].insert(neiPoint);
				edgeWeakItr = edgeWeak.find(neiPoint);
				if(edgeWeakItr == edgeWeak.end())
				{
					set<int> s;
					edgeWeak[neiPoint] = s;
				}
				edgeWeak[neiPoint].insert(point); 
			}
			else if(stat == 3)
			{
				edgeUnknownItr = edgeUnknown.find(point);
				if(edgeUnknownItr == edgeUnknown.end())
				{
					set<int> s;
					edgeUnknown[point] = s;
				}
				edgeUnknown[point].insert(neiPoint);
				edgeUnknownItr = edgeUnknown.find(neiPoint);
				if(edgeUnknownItr == edgeUnknown.end())
				{
					set<int> s;
					edgeUnknown[neiPoint] = s;
				}
				edgeUnknown[neiPoint].insert(point); 				
			}
		}
	}
	cout << "Step2 Find edges between clusters---" << clock() - start_time << "seconds ---";
	start_time = clock();
	/*######################################################################*/
	if(edgeYes.size() > 0)
	{
		connComp();
		cout << "Step3 Found connected components---" << clock() - start_time << "seconds ---";
		start_time = clock();
		mergeAssignNewEdge(); //implement it
		cout << "Step4 Merge edges and assign new edge---" << clock() - start_time << "seconds ---";
		start_time = clock();
	}
	int iteration = 0;
	while(true)
	{
		cout << "Range queries : " << rangeQueryPerformed.size() << "\n";
		cout << "Iteration : " << iteration << "    No. of clusters : " << clusters.size() << "\n";
       	cout << "CoreList : " << coreList.size() << "\n";
       	cout << "BorderList : " << borderList.size() << "\n";
       	cout << "NoiseList : " << noiseList.size() + neiNoise.size() << "\n";
		if(stoppingCondition()) //implement it
		{
			/*######################################################################*/
			cout << "Step5 Checked stopping condition---" << clock() - start_time << "seconds ---";
			start_time = clock();
			calculateStatDegree();
			vector<int> betaPoints = calculateScore();
			if(betaPoints.size() == 0)
			{
				break;
			}
			cout << "Step6 get beta points for further range queries---" << clock() - start_time << "seconds ---";
			start_time = clock();
			/*######################################################################*/
			vector<int>::iterator beta_itr;
			set<int>::iterator range_itr;
			for(beta_itr = betaPoints.begin(); beta_itr != betaPoints.end(); beta_itr++)
			{
				int point = *beta_itr;
				range_itr = find(rangeQueryPerformed.begin(),rangeQueryPerformed.end(),point);
				set<int> listOfNeighbors;
				set<int>::iterator lon_itr;
				map<int,set<int> >::iterator nm_itr; //iterator for neighbourMap
				map<int, string>::iterator border_itr;
				map<int, string>::iterator core_itr;
				map<int,set<int> >::iterator cfpm_itr;
				if(range_itr != rangeQueryPerformed.end())
				{
					listOfNeighbors = neighbourMap[point];
				}
				else
				{
					listOfNeighbors = performRangeQuery(point);
				}
				if(listOfNeighbors.size() < minPts)
				{
					borderList[point] = "PROCESSED";
					for(lon_itr = listOfNeighbors.begin(); lon_itr != listOfNeighbors.end(); lon_itr++)
					{
						int nei = *lon_itr;
						range_itr = find(rangeQueryPerformed.begin(),rangeQueryPerformed.end(),nei);
						if(range_itr == rangeQueryPerformed.end())
						{
							nm_itr = neighbourMap.find(nei);
							if(nm_itr == neighbourMap.end())
							{
								set<int> s;
								neighbourMap[nei] = s;
							}
							neighbourMap[nei].insert(point);
						}
					}
				}
				else
				{
					coreList[point] = "PROCESSED";
					border_itr = borderList.find(point);
					if(border_itr != borderList.end())
					{
						borderList.erase(border_itr);
					}
					for(lon_itr = listOfNeighbors.begin(); lon_itr != listOfNeighbors.end(); lon_itr++)
					{
						int nei = *lon_itr;
						cfpm_itr = coreForPointMap.find(nei);
						if(cfpm_itr == coreForPointMap.end())
						{
							set<int> s;
							coreForPointMap[nei] = s;
						}
						coreForPointMap[nei].insert(point);
						range_itr = find(rangeQueryPerformed.begin(),rangeQueryPerformed.end(),nei);
						if(range_itr == rangeQueryPerformed.end())
						{
							nm_itr = neighbourMap.find(nei);
							if(nm_itr == neighbourMap.end())
							{
								set<int> s;
								neighbourMap[nei] = s;
							}
							neighbourMap[nei].insert(point);
							if(neighbourMap[nei].size() < minPts)
							{
								border_itr = borderList.find(nei);
								if(border_itr == borderList.end())
								{
									range_itr = find(rangeQueryPerformed.begin(),rangeQueryPerformed.end(),nei);
									if(range_itr != rangeQueryPerformed.end())
									{
										borderList[nei] = "PROCESSED";
									}
									else
									{
										borderList[nei] = "UNPROCESSED";
									}
								}
							}
							else
							{
								core_itr = coreList.find(nei);
								if(core_itr == coreList.end())
								{
									coreList[nei] = "UNPROCESSED";
									border_itr = borderList.find(nei);
									if(border_itr != borderList.end())
									{
										borderList.erase(border_itr);
									}
								}
							}
						}
					}
					createPCIR(point);
					map<int, set<int> >::iterator clus_itr;
					for(clus_itr = clusters.begin(); clus_itr != clusters.end(); clus_itr++)
					{
						int repCore = clus_itr->first;
						if(point != repCore)
						{
							dccBetPCLU(repCore,point);
						}
					}
				}
			}
			cout << "Step7 rangeQuery on beta points, and create clusters and edges---" << clock() - start_time << "seconds ---";
			start_time = clock();
			/*###############################################################*/
			if(edgeYes.size() > 0)
			{
				connComp();
				cout << "Step7 merge clusters : Conn Comp---" << clock() - start_time << "seconds ---";
				start_time = clock();
				mergeAssignNewEdge(); //implement it
				cout << "Step7 merge clusters : Merge Comp---" << clock() - start_time << "seconds ---";
				start_time = clock();
			}
			updateStates();
			cout << "Step8 update states of clusters if possible---" << clock() - start_time << "seconds ---";
			start_time = clock();
			if(edgeYes.size() > 0)
			{
				connComp();
				mergeAssignNewEdge();
				cout << "Step8 Merge any new edges that were formed---" << clock() - start_time << "seconds ---";
				start_time = clock();
			}
		}
		else
		{
			break;
		}
		iteration++;
	}
	processOutliers();
	cout << "Step9 Process outliers---" << clock() - start_time << "seconds ---";
	start_time = clock();
	if(edgeYes.size() > 0)
	{
		connComp();
		mergeAssignNewEdge();
		cout << "Step9 Process outliers, merge any new edges formed---" << clock() - start_time << "seconds ---";
		start_time = clock();
	}
	cout << "Range queries :" << rangeQueryPerformed.size() << "\n";
	cout << "Core list :" << coreList.size() << "\n";
	cout << "Border List :" << borderList.size() << "\n";
	cout << "Noise List :" << noiseList.size()+neiNoise.size() << "\n";
	set<int> range;
	set<int> diff;
	for(i = 0; i < num_records; i++)
	{
		range.insert(i);
	}
	set_difference(range.begin(),range.end(),rangeQueryPerformed.begin(),rangeQueryPerformed.end(),std::inserter(diff,diff.begin()));
	cout << "Unprocessed Points :" << diff.size() << "\n";
	cout << "len(clusters) :" << clusters.size() << "\n";
}

int main(int argc, char* argv[])
{
	clock_t s_time = clock();
	char* file_name = argv[1];
	dataSet = readData(file_name);
	cout << "Reading data---" << clock() - start_time << "seconds ---";
	start_time = clock();
	anyDBC();
	cout << "---" << clock() - s_time << "seconds ---";
	return 0;
}