#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <cfloat>
#include <queue>

using namespace std;

vector<pair<vector<int>, double>> memo;
vector<pair<vector<int>, double>> rest;
int elementsNum = 20;

pair<vector<int>, double> max_two({}, -FLT_MAX);
pair<vector<int>, double> max_three({}, -FLT_MAX);
pair<vector<int>, double> max_four({}, -FLT_MAX);

double Upper = 0;

void insert_combs(vector<int> &arr, vector<int> &data, vector<vector<int>> &combs, int start, int end, int index, int r)  
{  
    if (index == r)  
    {  
    	vector<int> v;
        for (int j = 0; j < r; j++)  
            v.push_back(data.at(j));

        combs.push_back(v);
        return;  
    }  
  
    for (int i = start; i <= end && end - i + 1 >= r - index; i++)  
    {  
        data.at(index) = arr.at(i);  
        insert_combs(arr, data, combs, i+1, end, index+1, r);  
    }  
}

vector<vector<int>> units_pairs(vector<int> &arr, int n, int r)
{
	vector<vector<int>> combs;
	vector<int> data(r);
	insert_combs(arr, data, combs, 0, n-1, 0, r-1);
	insert_combs(arr, data, combs, 0, n-1, 0, r);;
	return combs;
}

vector<int> search_indexes(vector<vector<int>> &v)
{
	vector<int> units;
	vector<vector<int>> pairs;
	vector<int> indexes;


	for (int i = 0; i < v.size(); ++i)
	{
		if (v.at(i).size() == 1)
			units.push_back(v.at(i).at(0));
		else
			pairs.push_back(v.at(i));
	}

	for (int i = 0; i < memo.size(); ++i)
	{
		if (memo.at(i).first.size() == 1)
		{
			if (count(units.begin(), units.end(), memo.at(i).first.at(0)))
				indexes.push_back(i);
		}
		else if ( memo.at(i).first.size() == 2 )
		{
			if ( count(pairs.begin(), pairs.end(), memo.at(i).first))
				indexes.push_back(i);
		}
	}

	return indexes;
}


double calc_value(vector<int> &v)
{
	double x, value = 0;;

	if (v.size() == 1)
		value = rand() % 6 + 10; 
	else if (v.size() == 2)
	{
		int digit1 = v.at(0);
		int digit2 = v.at(1);

		vector<vector<int>> digits{{digit1}, {digit2}};
		vector<int> indexes{search_indexes(digits)};

		int index1 = indexes.at(0);
		int index2 = indexes.at(1);

		value = rand() % 5;

		if(value > max_two.second)
		{
			max_two.first.clear();
			for (int i = 0; i < v.size(); ++i)
				max_two.first.push_back(v.at(i));
			max_two.second = value;
			
		}
	}
	else if (v.size() == 3)
	{
		vector<vector<int>> combs = units_pairs(v, v.size(), 2);
		vector<int> indexes{search_indexes(combs)};	// size of indexes is 6. First 3 are indexes of units. Last 3 are indexes of pairs
		value = memo.at(indexes.at(3)).second + memo.at(indexes.at(4)).second + memo.at(indexes.at(5)).second
				+ memo.at(indexes.at(0)).second + memo.at(indexes.at(1)).second + memo.at(indexes.at(2)).second ; 

		if(value > max_three.second)
		{
			max_three.first.clear();
			for (int i = 0; i <v.size(); ++i)
				max_three.first.push_back(v.at(i));
			max_three.second = value;
			
		}
	}
	else
	{
		vector<vector<int>> combs = units_pairs(v, v.size(), 2);
		vector<int> indexes{search_indexes(combs)};	// size of indexes is 10. First 4 are indexes of units. Last 6 are indexes of pairs
		value = memo.at(indexes.at(4)).second + memo.at(indexes.at(5)).second + memo.at(indexes.at(6)).second
				+ memo.at(indexes.at(7)).second + memo.at(indexes.at(8)).second + memo.at(indexes.at(9)).second
				+ memo.at(indexes.at(0)).second + memo.at(indexes.at(1)).second + memo.at(indexes.at(2)).second + memo.at(indexes.at(3)).second;

		if(value > max_four.second)
		{
			max_four.first.clear();
			for (int i = 0; i <v.size(); ++i)
				max_four.first.push_back(v.at(i));
			max_four.second = value;
			
		}
	}

	return value;
}

/********** Functions to find Combinations of possible Coalitions **********/
void formCombination(vector<int> &arr, vector<int> &data, int start, int end, int index, int r)  
{  
    if (index == r)  
    {  
    	vector<int> v;
        for (int j = 0; j < r; j++)  
            v.push_back(data.at(j));

       if (r < 3)
        	memo.push_back(make_pair(v, calc_value(v)));
        else
        	rest.push_back(make_pair(v, calc_value(v)));
        return;  
    }  
  
    for (int i = start; i <= end && end - i + 1 >= r - index; i++)  
    {  
        data.at(index) = arr.at(i);  
        formCombination(arr, data, i+1, end, index+1, r);  
    }  
}
    
void coalCombinations(vector<int> &arr, int n, int r)  
{  
    vector<int> data(r);
    formCombination(arr, data, 0, n-1, 0, r);  
}  
/********** Functions to find Combinations of possible Coalitions **********/

void exaustive_search(vector<int> &elements)
{
	coalCombinations(elements, elements.size(), 1);
	coalCombinations(elements, elements.size(), 2);
	coalCombinations(elements, elements.size(), 3);
	coalCombinations(elements, elements.size(), 4);

	// printf("========================\n");
	// for (int i = 0; i < memo.size(); ++i)
	// {
	// 	for (int j = 0; j < memo.at(i).first.size(); ++j)
	// 		printf("%d ", memo.at(i).first.at(j));
	// 	printf("\t | value: %.2lf \n", memo.at(i).second);
	// }
	// printf("========================\n");

	// printf("\nrest:\n");
	// for (int i = 0; i < rest.size(); ++i)
	// {
	// 	for (int j = 0; j < rest.at(i).first.size(); ++j)
	// 		printf("%d ", rest.at(i).first.at(j));
	// 	printf("   | value: %.2lf \n", rest.at(i).second);
	// }

	printf("\nMax utility pair: \n");
	for (int i = 0; i < max_two.first.size(); ++i)
		printf("%d ", max_two.first.at(i));
	printf(" , value = %.2lf", max_two.second);

	// printf("\nMax utility threesome: \n");
	// for (int i = 0; i < max_three.first.size(); ++i)
	// 	printf("%d ", max_three.first.at(i));
	// printf(" , value = %.2lf", max_three.second);

	printf("\nMax utility Quarteto: \n");
	for (int i = 0; i < max_four.first.size(); ++i)
		printf("%d ", max_four.first.at(i));
	printf(" , value = %.2lf", max_four.second);

	printf("\n\n\n\n");
}

/************ For Exaustive Search *************/

vector<int> find_indexes(vector<vector<int>> &v)
{
	vector<int> units;
	vector<vector<int>> pairs;
	vector<int> indexes;


	for (int i = 0; i < v.size(); ++i)
	{
		if (v.at(i).size() == 1)
			units.push_back(v.at(i).at(0));
		else
			pairs.push_back(v.at(i));
	}

	for (int i = 0; i < memo.size(); ++i)
	{
		if (memo.at(i).first.size() == 1)
		{
			if (count(units.begin(), units.end(), memo.at(i).first.at(0)))
				indexes.push_back(i);
		}
		else
		{
			if (count(pairs.begin(), pairs.end(), memo.at(i).first))
				indexes.push_back(i);
			else
			{
				vector<int> p{memo.at(i).first.at(1)};
				p.push_back(memo.at(i).first.at(0));
				if (count(pairs.begin(), pairs.end(), p))
					indexes.push_back(i);
			}
		}
	}

	return indexes;
}

struct treeNode {
    treeNode* parent;
    int ID;
    int element;
    double value;
    pair<vector<int>, int> U;	
    vector<int> excluders;
    vector<int> members;
    bool blocked;
};


// This function takes a set of agents as an input. And returns the value of their coalition
double coal_value(vector<int> v)
{
	double value = 0;

	if (v.empty())
		return value;

	if ( v.size() == 1 )
	{
		for (unsigned i = 0; i < elementsNum; ++i)
		{
			if ( memo.at(i).first == v )
				return memo.at(i).second;
		}
	}
	else if ( v.size() == 2 )
	{
		vector<int> p{v.at(1)};
		p.push_back(v.at(0));

		for (unsigned i = elementsNum; i < memo.size(); ++i)
		{
			if ( memo.at(i).first == v || memo.at(i).first == p )
				return memo.at(i).second;
		}
	}
	else if ( v.size() == 3 )
	{
		vector<vector<int>> combs = units_pairs(v, v.size(), 2);
		vector<int> indexes{find_indexes(combs)};	// size of indexes is 6. First 3 are indexes of units. Last 3 are indexes of pairs
		value = memo.at(indexes.at(3)).second + memo.at(indexes.at(4)).second + memo.at(indexes.at(5)).second
				- ( memo.at(indexes.at(0)).second + memo.at(indexes.at(1)).second + memo.at(indexes.at(2)).second );
	}
	else if (v.size() == 4 )
	{
		vector<vector<int>> combs = units_pairs(v, v.size(), 2);
		vector<int> indexes{find_indexes(combs)};	// size of indexes is 10. First 4 are indexes of units. Last 6 are indexes of pairs
		value = memo.at(indexes.at(4)).second + memo.at(indexes.at(5)).second + memo.at(indexes.at(6)).second
				+ memo.at(indexes.at(7)).second + memo.at(indexes.at(8)).second + memo.at(indexes.at(9)).second
				- 2*( memo.at(indexes.at(0)).second + memo.at(indexes.at(1)).second + memo.at(indexes.at(2)).second + memo.at(indexes.at(3)).second );
	}

	return value;
}

/*
This function takes two sets of agents as an input. The first set includes agents who have to be members of a coalition.
The second set includes agents who have to be left outside Coalition. It return an estimation for the maximum value that
can be reached
*/
pair<vector<int>, int> upperEstimation(vector<int> in, vector<int> out, vector<int> &s)
{
	pair<vector<int>, int> estimation({}, 0);

	vector<int> p;
	int max;

	if (in.empty())
	{
		for (unsigned i = 0; i < s.size(); ++i)
		{
			if ( !count(out.begin(), out.end(), s.at(i)) )
			{
				in.push_back(s.at(i));
				break;
			}
		}
	}

	if ( in.size() == 1 )
	{
		int el = in.at(0);
		vector<int> v;
		int sec;

		max = 0;

		for (unsigned i = elementsNum; i < memo.size(); ++i)	// Check only Pairs - Skip Units
		{
			int elem1 = memo.at(i).first.at(0);
			int elem2 = memo.at(i).first.at(1);

			if ( el != elem1 && el != elem2 )
				continue;

			if ( count(out.begin(), out.end(), elem1) || count(out.begin(), out.end(), elem2) )
				continue;

			if ( memo.at(i).second > max )
			{
				max = memo.at(i).second;
				v = memo.at(i).first;
			}
		}

		max = 0;

		for (unsigned i = elementsNum; i < memo.size(); ++i)	// Check only Pairs - Skip Units
		{
			int elem1 = memo.at(i).first.at(0);
			int elem2 = memo.at(i).first.at(1);

			if ( count(v.begin(), v.end(), elem1) || count(v.begin(), v.end(), elem2) )
				continue;

			if ( count(out.begin(), out.end(), elem1) || count(out.begin(), out.end(), elem2) )
				continue;

			p = memo.at(i).first;
			p.insert(p.end(), v.begin(), v.end());
			int val = coal_value(p);

			if ( val > max )
			{
				max = val;
				estimation.first = p;
			}
		}

		estimation.second = max;
	}
	else if ( in.size() == 2 )
	{
		max = 0;

		for (unsigned i = elementsNum; i < memo.size(); ++i)	// Check only Pairs - Skip Units
		{
			int elem1 = memo.at(i).first.at(0);
			int elem2 = memo.at(i).first.at(1);

			if ( count(in.begin(), in.end(), elem1) || count(in.begin(), in.end(), elem2) )
				continue;

			if ( count(out.begin(), out.end(), elem1) || count(out.begin(), out.end(), elem2) )
				continue;

			p = memo.at(i).first;
			p.insert(p.end(), in.begin(), in.end());
			int val = coal_value(p);

			if ( val > max )
			{
				max = val;
				estimation.first = p;
			}
		}

		estimation.second = max;
	}
	else if ( in.size() == 3 ) 
	{
		max = 0;

		for (unsigned i = 0; i < elementsNum; ++i)
		{
			int elem = memo.at(i).first.at(0);

			if ( count(in.begin(), in.end(), elem) || count(out.begin(), out.end(), elem) )
				continue;

			p = in;
			p.push_back(elem);
			int val = coal_value(p);

			if ( val > max )
			{
				max = val;
				estimation.first = p;
			}
		}

		estimation.second = max;
	}

	return estimation;
}

/*
This function creates a new node. We determine it's parent node, ID, current list of nodes included and nodes excluded. 
If parameter included = true, we include current node in solution. Otherwise, this current node will be excluded
*/
struct treeNode* newNode(int id, int memb, treeNode* p, bool included)
{
    treeNode* n = new treeNode;
    (*n).ID = id;
    (*n).blocked = false;
    (*n).parent = p;
    (*n).element = memb;

    if (p == 0)	// For root only
    {
    	(*n).members;
    	(*n).excluders;
    }			
    else
    {
    	(*n).members = (*p).members;
    	(*n).excluders = (*p).excluders;
    }

    if( included && p != 0 )
    	(*n).members.push_back(memb);
    else if ( !included && p != 0 )
    	(*n).excluders.push_back(memb);

    return n;
}

void nodes_parameters(treeNode* n, vector<int> &s)
{
	if ( (*n).members.size() > 4 )
	{
		(*n).blocked = true;
		return;
	}

	(*n).value = coal_value((*n).members);
	
	treeNode* parent = (*n).parent;

	if (parent == 0)
	{
		(*n).U = upperEstimation((*n).members, (*n).excluders, s);
		return;
	}

	if( (*n).members.size() == 4 )
	{
		(*n).U.first = (*n).members;
		(*n).U.second = (*n).value;
		return;
	}

	bool skip = true;

	for (unsigned i = 0; i < (*n).members.size(); ++i)
	{
		if ( !count( (*parent).U.first.begin(), (*parent).U.first.end(), (*n).members.at(i) ) )
			skip = false;
	}

	for (unsigned i = 0; i < (*n).excluders.size(); ++i)
	{
		if ( count( (*parent).U.first.begin(), (*parent).U.first.end(), (*n).excluders.at(i) ) )
			skip = false;
	}

	if (skip)
		(*n).U = (*parent).U;
	else
		(*n).U = upperEstimation((*n).members, (*n).excluders, s);
}

void print_node(treeNode* n)
{
	printf("\nNode: \t %d \n", (*n).ID);
	if ( (*n).parent != 0)
		printf("Parent:  %d \n", (*n).parent->ID);
	printf("Element: %d \n", (*n).element);
	printf("Value: \t %f \n", (*n).value);
	printf("Upper: \t %d { ", (*n).U.second);
	for (int i = 0; i < (*n).U.first.size(); ++i)
		printf("%d ", (*n).U.first.at(i));
	printf("} \n");
	printf("Including: [ ");
	for (unsigned j = 0; j < (*n).members.size(); ++j)
		printf("%d ", (*n).members.at(j));
	printf("]\n");
	printf("Excluding: [ ");
	for (unsigned j = 0; j < (*n).excluders.size(); ++j)
		printf("%d ", (*n).excluders.at(j));
	printf("]\n");
	string result = ((*n).blocked) ? "Blocked" : "Free";
	cout << result << endl;
}


vector<int> order_by_shapley(vector<int> &elements)
{
	double shapley[elements.size()];
	double shapley_copy[elements.size()];

	for (unsigned i = 0; i < elements.size(); ++i)
	{
		int el = elements.at(i);
		int edges = 0;

		for (unsigned i = 0; i < memo.size(); ++i)
		{
			if ( memo.at(i).first.size() == 1 )
			{
				if ( memo.at(i).first.at(0) == el ){
					shapley[i] = memo.at(i).second;
					shapley_copy[i] = memo.at(i).second;
				}
			}
			else
			{
				if ( memo.at(i).first.at(0) == el || memo.at(i).first.at(1) == el )
					edges += memo.at(i).second;
			}
		}

		shapley[i] += 0.5*edges;
		shapley_copy[i] += 0.5*edges;
	}

	// for (int i = 0; i < elements.size(); ++i)
	// 	printf("\nElement %d has shapley: %f", elements.at(i), shapley[i]);

	vector<int> sorted;

	for (unsigned i = 0; i < elements.size(); ++i)
	{
		double max = 0;
		int idx = 0;

		for (int j = 0; j < elements.size(); ++j)
		{
			if (shapley_copy[j] > max)
			{
				max = shapley_copy[j];
				idx = j;
			}
		}

		// sorted.push_back(elements.at(idx));
		sorted.push_back(idx);
		shapley_copy[idx] = -1;
	}

	for (unsigned i = 0; i < sorted.size()-1; i++)    
    {
    	for (unsigned j = 0; j < sorted.size()-i-1; j++) 
    	{
    		if (shapley[sorted.at(j)] == shapley[sorted.at(j+1)])
			{
				if ( memo.at(sorted.at(j)).second < memo.at(sorted.at(j+1)).second )
					iter_swap(sorted.begin() + j, sorted.begin() + (j+1));
			}
    	} 
    }  

    int idx;

    for (int i = 0; i < sorted.size(); ++i)
    {
    	idx = sorted.at(i);
    	sorted.at(i) = elements.at(idx);
    }

	return sorted;
}

void recalculatePairs()
{
	for (unsigned i = 0; i < memo.size(); ++i)
	{
		if (memo.at(i).first.size() == 1)
			continue;

		int digit1 = memo.at(i).first.at(0);
		int digit2 = memo.at(i).first.at(1);

		vector<vector<int>> digits{{digit1}, {digit2}};
		vector<int> indexes{find_indexes(digits)};

		int index1 = indexes.at(0);
		int index2 = indexes.at(1);

		memo.at(i).second += memo.at(index1).second + memo.at(index2).second; 
	}
}



void branchNbound(vector<int> &elements)
{
	vector<int> sorted = order_by_shapley(elements);	// Elements sorted by ther shapley value
	recalculatePairs();

	printf("\nSorted:\n");
	for (int i = 0; i < sorted.size(); ++i)
		printf("%d ", sorted.at(i));

	printf("\n\n========================\n");
	for (int i = 0; i < memo.size(); ++i)
	{
		for (int j = 0; j < memo.at(i).first.size(); ++j)
			printf("%d ", memo.at(i).first.at(j));
		printf("\t | value: %.2lf \n", memo.at(i).second);
	}
	printf("========================\n\n");

	vector<pair<int, int>> line{make_pair(-1, sorted.at(0))};
	for (int i = 0; i < sorted.size(); ++i)
    {
    	if ( i == sorted.size() - 1)
    	{
    		line.push_back(make_pair(sorted.at(i), -1));
    		break;
    	}

    	line.push_back(make_pair(sorted.at(i), sorted.at(i+1)));
    }

	int idx = 0;

	treeNode* root = newNode(idx, -1, NULL, false);
	queue<treeNode*> extend;
    extend.push(root);

    vector<int> bestCoal;
    int maxVal = -1;

    while( !extend.empty() )
    {
    	treeNode* N = extend.front();
    	extend.pop();

    	printf("\nextend size: %ld\n", extend.size());

    	nodes_parameters(N, sorted);
    	print_node(N);

    	if ( !(*N).blocked && (*N).U.second >= maxVal )
    	{
    		printf("\nIn if.. (*N).U.second: %d , maxVal: %d\n", (*N).U.second, maxVal);
    		if ( (*N).value > maxVal )
    		{
    			bestCoal = (*N).members;
    			maxVal = (*N).value;
    		}

    		int next;
    		for (unsigned i = 0; i < line.size(); ++i)
    		{
    			if (line.at(i).first == (*N).element )
    				next = line.at(i).second;
    		}

			if (next == -1)
    			break;

		    idx++;
		    treeNode* n1 = newNode(idx, next, N, true);
		    idx++;
		    treeNode* n2 = newNode(idx, next, N, false);

		    extend.push(n1);
		    extend.push(n2);
    	}
    }

    printf("\nBest coalition: \n{ ");
    for (int i = 0; i < bestCoal.size(); ++i)
    	printf("%d ", bestCoal.at(i));
    printf("}  , value: %d\n", maxVal);
}

int main() 
{
	srand ( unsigned ( time(0) ) );

	vector<int> elements;

	for (int i = 0; i < elementsNum; ++i)
		elements.push_back(i);

	random_shuffle ( elements.begin(), elements.end() );

	exaustive_search(elements);
	branchNbound(elements);
}