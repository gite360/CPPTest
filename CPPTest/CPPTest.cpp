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
#include <math.h>
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

struct CoeffcientPair {
	int HDWT_id = NULL;
	DataType HDWT_original_value = NULL;
	DataType HDWT_normalized_value = NULL;
	DataType HDWT_fabs_value = NULL;

};

struct cmpLess {
	bool operator () (const CoeffcientPair &a, const CoeffcientPair &b) {
		return a.HDWT_fabs_value < b.HDWT_fabs_value;             // 从大到小  
	}
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

inline bool isPowerOf2(const DataType& old_length) {
	return getNextPowerOf2(unsigned int(old_length)) == old_length;
}

bool padZero(DataType* temp_original_array, const DataType& old_length, unsigned int& new_length, DataType* after_padding_array) {

	int i = 0;
	while (i < new_length) {
		if (i < old_length) {
			after_padding_array[i] = temp_original_array[i];
		}
		else {
			after_padding_array[i] = 0.0;
		}
		i++;
	}

	return true;
}

bool mergeSegments(DataType* wavelet_transform_time_series, unsigned int power_of_2, DataType M, int &retained_coeffs_length) {

	if (power_of_2 < M) {
		cout << "Something wrong! the number of segments is greater than M!!!" << endl;
		return false;
	}
	else if (power_of_2 == M) {
		return true;
	}
	else {
		return true;
	}
}


void getHDWT(unsigned int& power_of_2, DataType M, DataType* wavelet_transform_time_series, priority_queue<DataType>& fq_truncate_index) {

	//cout << power_of_2 << endl;
	int f_interation_times = int(log2(power_of_2)) - 1;//3
	//cout << f_interation_times << endl;
	int temp_power_of_2 = power_of_2 >> 1;

	DataType *fd_temp = new DataType[power_of_2];
	priority_queue<CoeffcientPair, vector<CoeffcientPair>, cmpLess> fq_HDWT_coefficients;

	CoeffcientPair temp_Coefficient;

	for (int i = 0; i < power_of_2; i++) {
		fd_temp[i] = wavelet_transform_time_series[i];
	}

	for (int i = f_interation_times; i >= 0; i--) {//3

		for (int j = 0; j < temp_power_of_2; j++) {
			wavelet_transform_time_series[j] = (fd_temp[j << 1] + fd_temp[1 + (j << 1)]) / 2;
			fd_temp[j] = (fd_temp[j << 1] + fd_temp[1 + (j << 1)]) / 2;
			wavelet_transform_time_series[j + temp_power_of_2] = wavelet_transform_time_series[j] - fd_temp[1 + (j << 1)];

			temp_Coefficient.HDWT_original_value = wavelet_transform_time_series[j + temp_power_of_2];
			temp_Coefficient.HDWT_normalized_value = temp_Coefficient.HDWT_original_value / pow(2, double(i) / 2);
			temp_Coefficient.HDWT_fabs_value = fabs(temp_Coefficient.HDWT_normalized_value);
			temp_Coefficient.HDWT_id = j + temp_power_of_2;
			fq_HDWT_coefficients.push(temp_Coefficient);
			//cout << double(i) / 2 <<" "<< pow(2, double(i) / 2) << endl;
		}
		temp_power_of_2 >>= 1;
	}

	temp_Coefficient.HDWT_original_value = *wavelet_transform_time_series;
	temp_Coefficient.HDWT_normalized_value = *wavelet_transform_time_series;
	temp_Coefficient.HDWT_fabs_value = abs(temp_Coefficient.HDWT_normalized_value);
	temp_Coefficient.HDWT_id = 0;
	fq_HDWT_coefficients.push(temp_Coefficient);

	for (int i = 0; i < M; i++) {
		fq_truncate_index.push(fq_HDWT_coefficients.top().HDWT_id);
		fq_HDWT_coefficients.pop();
	}

	/*while (!fq_HDWT_coefficients.empty()) {
		cout << fq_HDWT_coefficients.top().HDWT_id << " " << fq_HDWT_coefficients.top().HDWT_original_value << "    " << fq_HDWT_coefficients.top().HDWT_normalized_value << "    " << fq_HDWT_coefficients.top().HDWT_fabs_value << endl;
		fq_HDWT_coefficients.pop();
	}*/

	delete[] fd_temp;
}

void reconstructAPCA(DataType* wavelet_transform_time_series, const int& retained_coeffs_length, DataType* apca_presentation) {

	DataType temp_apca_value = NULL;
	int loop_times = log2(retained_coeffs_length);
	int inner_loop_times = 1;
	queue<DataType> temp_original_queue;

	for (int i = 0; i < retained_coeffs_length; i++) {
		apca_presentation[i] = wavelet_transform_time_series[i];
		//temp_original_queue.push(wavelet_transform_time_series[i]);
	}

	temp_original_queue.push(*wavelet_transform_time_series);

	for (int i = 0; i < loop_times; i++) {
		for (int j = 0; j < inner_loop_times; j++) {

			apca_presentation[(j << 1) + 1] = temp_original_queue.front() - apca_presentation[j + inner_loop_times];
			apca_presentation[j << 1] = 2 * temp_original_queue.front() - apca_presentation[(j << 1) + 1];

			temp_original_queue.pop();
			temp_original_queue.push(apca_presentation[j << 1]);
			temp_original_queue.push(apca_presentation[(j << 1) + 1]);

		}

		inner_loop_times <<= 1;
	}

	/*for (int i = 0; i < retained_coeffs_length; i++) {
		cout << apca_presentation[i] << endl;
	}*/
}


void reconstructApproximateAPCA(DataType* wavelet_transform_time_series, unsigned int& power_of_2, const int& retained_coeffs_length, priority_queue<DataType> fq_truncate_index, DataType* apca_presentation) {

	DataType temp_apca_value = NULL;
	unsigned int* fpa_Cr_value = new unsigned int[retained_coeffs_length];
	int loop_times = log2(retained_coeffs_length);
	int inner_loop_times = 1;
	queue<DataType> temp_original_queue;

	for (int i = 0; i < retained_coeffs_length; i++) {
		apca_presentation[i] = 0.0;
		fpa_Cr_value[i] = power_of_2;
		//temp_original_queue.push(wavelet_transform_time_series[i]);
	}

	for (int i = 0; i < retained_coeffs_length - 1; i++) {
		apca_presentation[int(fq_truncate_index.top())] = wavelet_transform_time_series[int(fq_truncate_index.top())];
		fq_truncate_index.pop();
		//temp_original_queue.push(wavelet_transform_time_series[i]);
	}

	//for (int i = 0; i < retained_coeffs_length; i++) {
	//	apca_presentation[i] = wavelet_transform_time_series[i];
	//	//temp_original_queue.push(wavelet_transform_time_series[i]);
	//}

	temp_original_queue.push(*wavelet_transform_time_series);

	for (int i = 0; i < loop_times; i++) {
		for (int j = 0; j < inner_loop_times; j++) {

			apca_presentation[(j << 1) + 1] = temp_original_queue.front() - apca_presentation[j + inner_loop_times];
			apca_presentation[j << 1] = 2 * temp_original_queue.front() - apca_presentation[(j << 1) + 1];
			
			fpa_Cr_value[j << 1] = ((j << 1) + 1)*power_of_2 / (inner_loop_times<<1)-1;
			fpa_Cr_value[(j << 1) + 1] = fpa_Cr_value[j << 1]+ power_of_2 / (inner_loop_times << 1);
			

			temp_original_queue.pop();
			temp_original_queue.push(apca_presentation[j << 1]);
			temp_original_queue.push(apca_presentation[(j << 1) + 1]);

		}

		inner_loop_times <<= 1;
	}

	for (int i = 0; i < retained_coeffs_length; i++) {
		cout << fpa_Cr_value[i] <<" "<< apca_presentation[i] << endl;
	}


	delete[] fpa_Cr_value;
}





bool getExactMeanValue(DataType* wavelet_transform_time_series, DataType &n) {


	return true;

}

void getAPCAPoint(DataType* orginal_time_series, APCA &italicC, DataType &n, const int &N) {

	if (n < N) {                                          //n must longer than N
		cout << "ERROR: N has to smaller than n!!!!!!!";
	}

	int i = 0, j = 0, retained_coeffs_length = NULL;
	unsigned int power_of_2 = NULL;
	DataType* wavelet_transform_time_series;
	DataType* apca_presentation;
	priority_queue<DataType> fq_truncate_index;

	power_of_2 = getNextPowerOf2(n);
	wavelet_transform_time_series = new DataType[power_of_2];
	padZero(orginal_time_series, n, power_of_2, wavelet_transform_time_series);
	getHDWT(power_of_2, N, wavelet_transform_time_series, fq_truncate_index);
	retained_coeffs_length = fq_truncate_index.top() + 1;
	apca_presentation = new DataType[retained_coeffs_length];
	reconstructApproximateAPCA(wavelet_transform_time_series, power_of_2, retained_coeffs_length, fq_truncate_index, apca_presentation);

	delete[] apca_presentation;
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


	DataType test_array[8] = { 7,5,5,3,3,3,4,6 };
	//DataType* mp_HDWT_time_series = new DataType[8];
	DataType n = 8, M = 3;
	APCA test_apca;
	unsigned int after_padding_length = NULL;
	int retained_coeffs_length = NULL;


	//getHDWT(test_array, getNextPowerOf2(8), M, mp_HDWT_time_series, retained_coeffs_length);
	getAPCAPoint(test_array, test_apca, n, M);

	for (int i = 0; i < 8; i++) {
		//cout << mp_HDWT_time_series[i] << " ";
	}
	cout << endl;

	//delete[] mp_HDWT_time_series;
	system("pause");
	return 0;
}

