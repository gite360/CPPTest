﻿// CPPTest.cpp : Defines the entry point for the console application.
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

#include <numeric>      // std::accumulate


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
		return a.HDWT_fabs_value < b.HDWT_fabs_value;             // from big to small  
	}
};

struct MERGE_INDEX {
	DataType distance_AE = NULL;
	int* index_array = nullptr;
	int index_array_length = NULL;
};

struct cmpMore {
	bool operator () (const MERGE_INDEX &a, const MERGE_INDEX &b) {
		return a.distance_AE > b.distance_AE;             // decrease
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


template<typename T>
void initialArray(T* array, int array_length) {

	for (int i = 0; i < array_length; i++) {
		array[i] = NULL;
	}
}

inline unsigned int getNextPowerOf2(unsigned int array_length) {

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

bool truncateZero(const DataType& old_length, unsigned int& new_length, APCA& apca_presentation) {

	if (new_length == old_length) {
		return true;
	}
	else if (new_length > old_length) {
		for (int i = apca_presentation.segmentNum - 1; i >= 0; i--) {
			if (old_length< apca_presentation.r[i] && new_length > apca_presentation.r[i - 1]) {
				apca_presentation.r[i] = old_length;
				return true;
			}
		}
		cout << "Something Wrong!!!" << endl;
		return false;
	}
	else if (new_length < old_length) {
		cout << "Something wrong! new_length is greater than old_length!!!" << endl;
		return false;
	}
	else if (getNextPowerOf2(unsigned int(old_length)) != new_length) {
		cout << "Something wrong! new_length is not Power of 2!!!" << endl;
		return false;
	}
	else {
		cout << "Something Wrong!!!" << endl;
		return false;
	}

}



void getHDWT(unsigned int& power_of_2, const int& M, DataType* wavelet_transform_time_series, priority_queue<DataType>& fq_truncate_index) {

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
		cout << " Big magnitude: " << fq_HDWT_coefficients.top().HDWT_id << endl;
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


void reconstructApproximateAPCA(DataType* wavelet_transform_time_series, unsigned int& power_of_2, const int& retained_coeffs_length, priority_queue<DataType> fq_truncate_index, APCA& apca_presentation) {

	DataType temp_apca_value = NULL;
	int loop_times = log2(retained_coeffs_length);
	int inner_loop_times = 1;
	queue<DataType> temp_original_queue;

	for (int i = 0; i < retained_coeffs_length; i++) {
		apca_presentation.v[i] = 0.0;
		apca_presentation.r[i] = power_of_2;
		//temp_original_queue.push(wavelet_transform_time_series[i]);
	}

	for (int i = 0; i < fq_truncate_index.size(); i++) {
		apca_presentation.v[int(fq_truncate_index.top())] = wavelet_transform_time_series[int(fq_truncate_index.top())];
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

			apca_presentation.v[(j << 1) + 1] = temp_original_queue.front() - apca_presentation.v[j + inner_loop_times];
			apca_presentation.v[j << 1] = 2 * temp_original_queue.front() - apca_presentation.v[(j << 1) + 1];

			apca_presentation.r[j << 1] = ((j << 1) + 1)*power_of_2 / (inner_loop_times << 1) - 1;
			apca_presentation.r[(j << 1) + 1] = apca_presentation.r[j << 1] + power_of_2 / (inner_loop_times << 1);

			temp_original_queue.pop();
			temp_original_queue.push(apca_presentation.v[j << 1]);
			temp_original_queue.push(apca_presentation.v[(j << 1) + 1]);

		}
		inner_loop_times <<= 1;
	}

	for (int i = 0; i < apca_presentation.segmentNum; i++) {
		cout << apca_presentation.r[i] << " " << apca_presentation.v[i] << endl;
	}

}


template <typename T>
T inline getAve(T* originalArray, const T &arrayLength) {
	return accumulate(originalArray, originalArray + int(arrayLength), 0) / arrayLength;
}

void getExactMeanValue(DataType* orginal_time_series, APCA& apca_presentation) {

	double fd_segment_length = NULL;
	int i = 0;


	apca_presentation.v[i] = getAve(orginal_time_series, apca_presentation.r[0] + 1);
	cout << apca_presentation.r[i] << " " << apca_presentation.v[i] << endl;

	for (i = 1; i < apca_presentation.segmentNum; i++) {
		apca_presentation.v[i] = getAve(orginal_time_series + int(apca_presentation.r[i - 1] + 1), double(apca_presentation.r[i] - apca_presentation.r[i - 1]));
		cout << apca_presentation.r[i] << " " << apca_presentation.v[i] << endl;;
	}

}

double distanceAE(const DataType* orginal_time_series, const APCA& merge_apca_presentation) {
	double reconstruction_error = 0;

	int i = 0;
	for (int j = 0; j <= merge_apca_presentation.r[0]; j++) {
		reconstruction_error += pow(merge_apca_presentation.v[i] - orginal_time_series[j], 2);
	}

	for (i = 1; i < merge_apca_presentation.segmentNum; i++) {
		for (int j = merge_apca_presentation.r[i - 1] + 1; j <= merge_apca_presentation.r[i]; j++) {
			reconstruction_error += pow(merge_apca_presentation.v[i] - orginal_time_series[j], 2);
		}
	}
	//cout << "distance AE: " << sqrt(reconstruction_error) << endl;
	return sqrt(reconstruction_error);
}


bool combinateSegments(const DataType* orginal_time_series, const APCA& merge_apca_presentation, int segment_number, const int& n, const int& M, int merge_number, int merge_index, const int& last_merge_point, vector<double>& result_vector) {

	if (merge_number == 0) {

		for (int i = 0; i < M; i++) {
			result_vector.push_back(merge_apca_presentation.r[i]);
		}
		result_vector.push_back(distanceAE(orginal_time_series, merge_apca_presentation));
		return true;

	}
	else {

		APCA temp_apca;
		temp_apca.r = new DataType[segment_number];
		temp_apca.v = new DataType[segment_number];
		temp_apca.segmentNum = merge_apca_presentation.segmentNum - 1;


		while (merge_index < merge_apca_presentation.segmentNum - 1) {

			if (merge_apca_presentation.r[merge_index] > last_merge_point) {
				for (int i = 0; i < merge_index; i++) {
					temp_apca.r[i] = merge_apca_presentation.r[i];
					temp_apca.v[i] = merge_apca_presentation.v[i];
				}

				temp_apca.v[merge_index] = (merge_apca_presentation.v[merge_index] + merge_apca_presentation.v[merge_index + 1]) / 2;
				temp_apca.r[merge_index] = merge_apca_presentation.r[merge_index + 1];

				int i = merge_index + 1;

				while (i < temp_apca.segmentNum) {
					temp_apca.r[i] = merge_apca_presentation.r[i + 1];
					temp_apca.v[i] = merge_apca_presentation.v[i + 1];
					i++;
				}


				int last_merge_point = merge_apca_presentation.r[merge_index];
				combinateSegments(orginal_time_series, temp_apca, segment_number - 1, n, M, merge_number - 1, 0, last_merge_point, result_vector);

			}
			merge_index++;
		}
		delete[] temp_apca.r;
		delete[] temp_apca.v;
	}

}

int getCombinationSum(int segment_number, const int& merge_frequency) {
	segment_number--;
	int combination_sum = 1;
	for (int i = 0; i < merge_frequency; i++) {
		combination_sum *= segment_number;
		segment_number--;
	}
	return combination_sum;
}


int getCombinationNumber(int a, int b) {
	int record[100][2];
	if (a == 0 && b == 0) return false;
	int V1, V2;
	V1 = 0, V2 = 1;

	record[0][V1] = 1; //C(1,0)=1
	record[1][V1] = 1; //C(1,1)=1

	for (int i = 2; i <= a; i++){
		record[0][V2] = record[i][V2] = 1;//C(i,0)=1,C(i,i)=1
		for (int j = 1; j < i; j++){
			record[j][V2] = record[j - 1][V1] + record[j][V1];
		}

		swap(V1, V2);
	}

	cout << record[b][V1] << endl;
	return record[b][V1];
}

bool mergeSegments(DataType* orginal_time_series,const APCA& apca_presentation,const DataType& n, const int& M) {

	if (apca_presentation.segmentNum > M) {
		int merge_frequency = apca_presentation.segmentNum - M;

		int segment_number = apca_presentation.segmentNum;
		int combination_sum = getCombinationNumber(apca_presentation.segmentNum - 1, merge_frequency);

		priority_queue<MERGE_INDEX, vector<MERGE_INDEX>, cmpMore> queue_reconstrution_error;
		vector<double> result_vector;

		combinateSegments(orginal_time_series, apca_presentation, segment_number, n, M, merge_frequency, 0, 0, result_vector);

		MERGE_INDEX* best_merge_index;
		MERGE_INDEX* p_array_begin_index;
		best_merge_index = new MERGE_INDEX[combination_sum];
		p_array_begin_index = best_merge_index;

		for (int i = 0; i < combination_sum; i++) {
			best_merge_index[i].index_array_length = M;
			best_merge_index[i].index_array = new int[M];
			initialArray(best_merge_index[i].index_array, M);
		}

		int  count = 0;

		for (vector<double>::iterator it = result_vector.begin(); it != result_vector.end(); ++it) {
			//std::cout << ' ' << *it;
			while (count < M) {
				best_merge_index->index_array[count] = *(it + count);
				count++;
			}
			best_merge_index->distance_AE = *(it + count);
			queue_reconstrution_error.push(*best_merge_index);
			best_merge_index++;
			count = 0;
			it += M;
		}

		while (!queue_reconstrution_error.empty()) {
			for (int i = 0; i < M; i++) {
				cout << queue_reconstrution_error.top().index_array[i] << " ";
			}
			cout << queue_reconstrution_error.top().distance_AE << endl;
			queue_reconstrution_error.pop();
		}

		best_merge_index = p_array_begin_index;

		for (int i = 0; i < combination_sum; i++) {
			best_merge_index[i].index_array = nullptr;
			delete[] best_merge_index[i].index_array;
		}
		best_merge_index = nullptr;
		delete[] best_merge_index;

	}
	else if (apca_presentation.segmentNum == M) {
		return true;
	}
	else if (apca_presentation.segmentNum < M) {
		cout << "Something wrong! the number of segments is greater than M!!!" << endl;
		return false;
	}
	else {
		cout << "Something Wrong!!!" << endl;
		return false;
	}

}


void getAPCAPoint(DataType* orginal_time_series, APCA &italicC, DataType &n, const int &N) {

	if (n < N) {                                          //n must longer than N
		cout << "ERROR: N has to smaller than n!!!!!!!";
	}

	int i = 0, j = 0, retained_coeffs_length = NULL;
	unsigned int power_of_2 = NULL;
	DataType* wavelet_transform_time_series;
	APCA general_apca_presentation;
	APCA best_apca_presentation;


	priority_queue<DataType> fq_truncate_index;


	power_of_2 = getNextPowerOf2(n);
	wavelet_transform_time_series = new DataType[power_of_2];

	padZero(orginal_time_series, n, power_of_2, wavelet_transform_time_series);
	getHDWT(power_of_2, N, wavelet_transform_time_series, fq_truncate_index);

	retained_coeffs_length = fq_truncate_index.top() + 1;
	general_apca_presentation.segmentNum = retained_coeffs_length;
	general_apca_presentation.v = new DataType[retained_coeffs_length];
	general_apca_presentation.r = new DataType[retained_coeffs_length];

	reconstructApproximateAPCA(wavelet_transform_time_series, power_of_2, retained_coeffs_length, fq_truncate_index, general_apca_presentation);
	truncateZero(n, power_of_2, general_apca_presentation);
	getExactMeanValue(orginal_time_series, general_apca_presentation);
	mergeSegments(orginal_time_series, general_apca_presentation, n, N);

	delete[] general_apca_presentation.r;
	delete[] general_apca_presentation.v;
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
	//	a[i] = double(rand()%101);   
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


	getCombinationNumber(7, 2);

	//delete[] mp_HDWT_time_series;
	system("pause");
	return 0;
}

