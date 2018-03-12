// CPPTest.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <algorithm>
#include<vector>  
#include<list> 
#include<queue>  
#include<functional>
#include<time.h>
#include<stdlib.h>
#include<stdio.h>

#include <fstream>
#include <sstream>
#include <string>


using namespace std;

typedef double DataType;

struct APCA {
	double* v = nullptr;
	double* r = nullptr;
	int segmentNum = NULL; //the dimension of the index(eg. MBR, imput parameter)

	/*~APCA() {
		if (v!=nullptr) {
			delete[] v;
			v = nullptr;
		}
		if (r != nullptr) {
			delete[] r;
			r = nullptr;
		}
	}*/
};


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



//DataType getNextPowerOf2(unsigned int array_length) {
//
//	if (array_length <= 0) return 0;
//
//	array_length--;
//	array_length |= array_length >> 1;
//	array_length |= array_length >> 2;
//	array_length |= array_length >> 4;
//	array_length |= array_length >> 8;
//	array_length |= array_length >> 16;
//
//	return DataType(array_length+1);
//}

unsigned int getNextPowerOf2(unsigned int array_length) {

	if (array_length <= 0) return 0;

	array_length--;
	array_length |= array_length >> 1;
	array_length |= array_length >> 2;
	array_length |= array_length >> 4;
	array_length |= array_length >> 8;
	array_length |= array_length >> 16;

	return array_length + 1;
}

void getHDWT( DataType* fa_temp_original_array, unsigned int power_of_2, DataType* wavelet_transform_time_series) {

	//cout << power_of_2 << endl;
	int f_interation_times = int(log2(power_of_2));//3
	//cout << f_interation_times << endl;
	int temp_power_of_2 = power_of_2>>1;

	DataType *temp = new DataType[power_of_2];


	for (int i = 0; i < power_of_2;i++) {
		temp[i] = fa_temp_original_array[i];
	}

	for (int i = f_interation_times; i >= 0; i--) {//3

		for (int j = 0; j < temp_power_of_2; j++) {
			wavelet_transform_time_series[j] = (temp[j << 1] + temp[1 + (j << 1)]) / 2;
			temp[j] = (temp[j << 1] + temp[1 + (j << 1)]) / 2;
			wavelet_transform_time_series[j + temp_power_of_2] = wavelet_transform_time_series[j] - temp[1 + (j << 1)];
		}

		//cout << temp_power_of_2 << endl;
		temp_power_of_2>>=1 ;
	}
	
	delete[] temp;
}


void getHDWT1(DataType* fa_temp_original_array,const unsigned int &power_of_2, list<DataType>& wavelet_transform_time_series) {

	int f_interation_times = int(log2(power_of_2));//3
	int f_loop_times= power_of_2 >> 1;
	DataType *f_average_arrayp = new DataType[f_loop_times];

	for (int i = 0; i < f_loop_times; i++) {
		cout << (fa_temp_original_array[i << 1] + fa_temp_original_array[(i << 1) + 1]) / 2 - fa_temp_original_array[(i << 1) + 1] << endl;
		wavelet_transform_time_series.push_back((fa_temp_original_array[i<<1]+ fa_temp_original_array[(i << 1)+1])/2- fa_temp_original_array[(i << 1) + 1]);
	}

}



void getAPCAPoint(const DataType* orginal_time_series, APCA &italicC, DataType &n, const int &N) {

	if (n < N) {                                          //n must longer than N
		cout << "ERROR: N has to smaller than n!!!!!!!";
	}
	int i = 0, j = 0;

	unsigned int power_of_2 = getNextPowerOf2(n);

	int zero_number = power_of_2 - n;

	DataType* fa_temp_original_array = new DataType[power_of_2];


	while (i<power_of_2) {
		if (i<n) {
			fa_temp_original_array[j] = orginal_time_series[j];
		}
		else {
			fa_temp_original_array[j] = 0.0;
		}
		i++;
	}

	DataType* wavelet_transform_time_series = new DataType[power_of_2];

	delete[] fa_temp_original_array;
	delete[] wavelet_transform_time_series;


}



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

	//double a[100][97];
	//double* pointer = &a[0][0];
	//string str;
	//string d = ",";
	//string ptr;

	//for (int i = 0; i < 100; i++) {
	//	for (int j = 0; j < 97; j++) {
	//		a[i][j] = NULL;
	//	}
	//}

	//ifstream fs = ifstream("TRAIN");

	//while (!fs.eof()) {
	//	//getline(fs, str, ',');
	//	fs >> str;
	//	stringstream sstr(str);
	//	while (getline(sstr, ptr, ',')) {
	//		*pointer = stod(ptr);
	//		 pointer++;
	//	}
	//	
	//}

	//fs.close();

	//for (int i = 0; i < 100; i++) {
	//	for (int j = 0; j < 97; j++) {
	//		cout << a[i][j] << " ";
	//	}
	//	cout << endl << endl;
	//}

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


	DataType test_array[8] = {7,5,5,3,3,3,4,6};
	DataType* mp_HDWT_time_series = new DataType[8];
	
	getHDWT(test_array,getNextPowerOf2(8),mp_HDWT_time_series);

	for (int i = 0; i < 8;i++) {
		cout << mp_HDWT_time_series[i] << " ";
	}

	delete[] mp_HDWT_time_series;
	system("pause");
	return 0;
}

