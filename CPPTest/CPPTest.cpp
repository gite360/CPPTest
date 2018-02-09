// CPPTest.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <IOSTREAM>
#include<vector>  
#include<queue>  
#include<functional> 
using namespace std;


struct Node{

	int key;
	int value;

};

struct priorityDecreasing {

	bool operator ()(const Node & a, const Node & b){
		return a.key < b.key;
	}
};

struct priorityIncrement{

	bool operator ()(const Node & a, const Node & b){
		return a.key > b.key;
	}
};


int main()
{
	const int len = 5;
	int i;


	priority_queue<Node, vector<Node>, priorityIncrement > pq2;
	
	vector<Node> b(len);

	b[0].key = 6; b[0].value = 1;
	b[1].key = 9; b[1].value = 5;
	b[2].key = 2; b[2].value = 3;
	b[3].key = 8; b[3].value = 2;
	b[4].key = 1; b[4].value = 4;

	for (i = 0; i < len; i++) {
		pq2.push(b[i]);
	}
		

	
	for (i = 0; i < len; i++)
	{
		cout << pq2.top().key << '\t' << pq2.top().value << endl;
		pq2.pop();
		
	}

	system("pause");
	return 0;
}

