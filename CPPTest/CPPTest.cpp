// CPPTest.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <algorithm>
#include<vector>  
#include<queue>  
#include<functional>
#include<time.h>
#include<stdlib.h>
#include<stdio.h>

#include <fstream>
#include <sstream>
#include <string>


using namespace std;


//struct Node{
//
//	int key;
//	int value;
//
//};
//
//struct priorityDecreasing {
//
//	bool operator ()(const Node & a, const Node & b){
//		return a.key < b.key;
//	}
//};
//
//struct priorityIncrement{
//
//	bool operator ()(const Node & a, const Node & b){
//		return a.key > b.key;
//	}
//};

int N = 2;

int main()
{
	/*const int len = 5;
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
	}*/

	//static int n = 6;
	//int *test_array=new int;
	//test_array[4] = 5;
	//cout<< test_array[4]<<endl;

	//srand((unsigned)time(NULL));
	//double a[1000];
	////vector<int> a(100);
	//printf("(The range is :0~100).\n");
	//for (int i = 0; i < 1000; i++) {
	//	a[i] = double(rand()%101);    //用rand函数生成0-100的随机数，并赋值给数组a[i]
	//	printf("%10f", a[i]);
	//	if (i % 10 == 0 && i != 0)
	//		printf("\n");
	//}
	//cout << endl;
	/*int num = 2;

	if (num & 1 == true) {
		cout <<num<< " is odd" << endl;
	}
	else {
		cout <<num<< " is even" << endl;
	}

	char* arry = new char[100];
	delete[]arry;
	arry = nullptr;
	delete[]arry;
	cout << std::min(10, 10) << endl;*/

	double a[100][97];
	double* pointer = &a[0][0];
	string str;
	double temp;
	string d = ",";
	string ptr;

	for (int i = 0; i < 100; i++) {
		for (int j = 0; j < 97; j++) {
			a[i][j] = NULL;
		}
	}

	ifstream fs = ifstream("TRAIN");

	while (!fs.eof()) {
		//getline(fs, str, ',');
		fs >> str;
		stringstream sstr(str);
		while (getline(sstr, ptr, ',')) {
			*pointer = stod(ptr);
			 pointer++;
		}
		
	}

	fs.close();

	for (int i = 0; i < 100; i++) {
		for (int j = 0; j < 97; j++) {
			cout << a[i][j] << " ";
		}
		cout << endl << endl;
	}

	/*while (!fs.eof()) {
		fs >> temp >> str;
		*pointer = temp;
		
		pointer++;
	}*/

	/*
	while (!fs.eof())
	{

		fs >> *pointer;
		pointer++;
	}*/

	/*while (fs >> *pointer) {
		pointer++;
	}*/

	/*double value=0;
	while (fs >> value)
	{
		*pointer = value;
		pointer++;
	}*/

	/*for (int i = 0; i < 100;i++) {
		for (int j = 0; j < 96;j++) {
			fs >> a[i][j];
		}
	}*/
	
	//getchar();
	system("pause");
	return 0;
}

