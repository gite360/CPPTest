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
//#include <cstring>

#include <numeric>      // std::accumulate


using namespace std;

typedef double DataType;

#define MAXMIUM 999999

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

struct COEFFICIENT_PAIR {
	int HDWT_id = NULL;
	DataType HDWT_original_value = NULL;
	DataType HDWT_normalized_value = NULL;
	DataType HDWT_magnitude = NULL;

};

struct cmpLess {
	bool operator () (const COEFFICIENT_PAIR &a, const COEFFICIENT_PAIR &b) {
		return a.HDWT_magnitude < b.HDWT_magnitude;             // from big to small  
	}
};

struct APCA_MERGE {
	DataType distance_AE = NULL;
	int* r = nullptr;
	DataType* v = nullptr;
};

struct cmpMore {
	bool operator () (const APCA_MERGE &a, const APCA_MERGE &b) {
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



void initialArray(APCA_MERGE& array, const int& array_length) {

	for (int i = 0; i < array_length; i++) {
		array.r[i] = NULL;
		array.v[i] = NULL;
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

bool truncateZero(const DataType& old_length, const unsigned int& new_length, APCA& apca_presentation) {
	cout << "segmentNum: " << apca_presentation.segmentNum;
	cout << ", new length: " << new_length << ", old length: " << old_length << endl;
	if (new_length == old_length) {
		return true;
	}
	else if (new_length > old_length) {
		int f_segments_last_point = apca_presentation.segmentNum - 1;

		int f_old_last_point = old_length - 1;

		for (int i = f_segments_last_point; i >= 0; i--) {
			if (f_old_last_point <= apca_presentation.r[i] && f_old_last_point > apca_presentation.r[i - 1]) {
				//cout <<"***"<<i<<" "<< apca_presentation.r[i] << endl;
				apca_presentation.r[i] = f_old_last_point;
				apca_presentation.segmentNum = i + 1;
				//cout << apca_presentation.r[i] << endl;
				return true;
			}
		}
		cout << "truncateZero(): Something Wrong!!!!!!!!!!!!!!!!!!!!!!" << endl;
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



void getHDWT(unsigned int& length_power_of_2, const int& M, DataType* wavelet_transform_time_series, priority_queue<DataType>& fq_truncate_index) {

	//cout << power_of_2 << endl;
	int f_interation_times = int(log2(length_power_of_2)) - 1;//3
	cout << "f_interation_times: " << f_interation_times << endl;
	int temp_power_of_2 = length_power_of_2 >> 1;

	DataType* fd_temp = new DataType[length_power_of_2];
	priority_queue<COEFFICIENT_PAIR, vector<COEFFICIENT_PAIR>, cmpLess> fq_HDWT_coefficients;

	COEFFICIENT_PAIR temp_Coefficient;

	for (int i = 0; i < length_power_of_2; i++) {
		fd_temp[i] = wavelet_transform_time_series[i];
	}

	for (int i = f_interation_times; i >= 0; i--) {//3

		for (int j = 0; j < temp_power_of_2; j++) {
			wavelet_transform_time_series[j] = (fd_temp[j << 1] + fd_temp[1 + (j << 1)]) / 2;
			fd_temp[j] = (fd_temp[j << 1] + fd_temp[1 + (j << 1)]) / 2;
			wavelet_transform_time_series[j + temp_power_of_2] = wavelet_transform_time_series[j] - fd_temp[1 + (j << 1)];

			temp_Coefficient.HDWT_original_value = wavelet_transform_time_series[j + temp_power_of_2];
			temp_Coefficient.HDWT_normalized_value = temp_Coefficient.HDWT_original_value / pow(2, double(i) / 2);
			temp_Coefficient.HDWT_magnitude = fabs(temp_Coefficient.HDWT_normalized_value);
			temp_Coefficient.HDWT_id = j + temp_power_of_2;

			cout << "(" << temp_Coefficient.HDWT_original_value << "," << temp_Coefficient.HDWT_magnitude << ") ";
			/*for (int i = 0; i < length_power_of_2; i++) {
				cout<<wavelet_transform_time_series[i] << " ";
			}
			cout << endl;*/
			/*for (int i = 0; i < length_power_of_2; i++) {
				cout << fd_temp[i] << " ";
			}
			cout << endl;*/

			fq_HDWT_coefficients.push(temp_Coefficient);
		}
		temp_power_of_2 >>= 1;
		cout << endl;
	}

	cout << "wavelet_transform_time_series: ";
	for (int i = 0; i < length_power_of_2; i++) {
		cout << wavelet_transform_time_series[i] << " ";
	}
	cout << endl;

	temp_Coefficient.HDWT_original_value = *wavelet_transform_time_series;
	temp_Coefficient.HDWT_normalized_value = *wavelet_transform_time_series;
	temp_Coefficient.HDWT_magnitude = abs(temp_Coefficient.HDWT_normalized_value);
	temp_Coefficient.HDWT_id = 0;
	fq_HDWT_coefficients.push(temp_Coefficient);

	for (int i = 0; i < M; i++) {
		cout << " Big magnitude(ID): " << fq_HDWT_coefficients.top().HDWT_id << endl;
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

	cout << " approximate length : " << endl;
	for (int i = 0; i < apca_presentation.segmentNum; i++) {
		cout << apca_presentation.r[i] << " " << apca_presentation.v[i] << endl;
	}

}


template <typename T>
T inline getAve(const T* originalArray, const T &arrayLength) {
	return accumulate(originalArray, originalArray + int(arrayLength), 0.0) / T(arrayLength);
}

void getExactMeanValue(const DataType* orginal_time_series, APCA& apca_presentation) {

	double fd_segment_length = NULL;
	int i = 0;

	cout << "Exact mean value: " << endl;
	apca_presentation.v[i] = getAve(orginal_time_series, apca_presentation.r[0] + 1);
	cout << apca_presentation.r[i] << " " << apca_presentation.v[i] << endl;

	for (i = 1; i < apca_presentation.segmentNum; i++) {
		apca_presentation.v[i] = getAve(orginal_time_series + int(apca_presentation.r[i - 1]+1), double(apca_presentation.r[i] - apca_presentation.r[i - 1]));
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
	cout << "distance AE: " << sqrt(reconstruction_error) << endl;
	return sqrt(reconstruction_error);
}


bool combinateSegments(const DataType* orginal_time_series, const APCA& merge_apca_presentation, int segment_number, const int& n, const int& M, int merge_number, int merge_index, const int& last_merge_point, vector<DataType>& result_vector) {

	if (merge_number == 0) {

		for (int i = 0; i < M; i++) {
			cout << merge_apca_presentation.r[i] << " ";
			result_vector.push_back(merge_apca_presentation.r[i]);
			result_vector.push_back(merge_apca_presentation.v[i]);
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
			//cout << merge_apca_presentation.r[merge_index] << " " << last_merge_point << endl;
			if (merge_apca_presentation.r[merge_index] >= last_merge_point) {
				//cout << ">>>>>>>" << endl;
				for (int i = 0; i < merge_index; i++) {
					temp_apca.r[i] = merge_apca_presentation.r[i];
					temp_apca.v[i] = merge_apca_presentation.v[i];
				}

				//temp_apca.v[merge_index] = (merge_apca_presentation.v[merge_index] + merge_apca_presentation.v[merge_index + 1]) / 2;
				temp_apca.v[merge_index] = getAve(orginal_time_series + int(temp_apca.r[merge_index - 1] + 1), DataType(merge_apca_presentation.r[merge_index + 1] - merge_apca_presentation.r[merge_index - 1]));
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

	for (int i = 2; i <= a; i++) {
		record[0][V2] = record[i][V2] = 1;//C(i,0)=1,C(i,i)=1
		for (int j = 1; j < i; j++) {
			record[j][V2] = record[j - 1][V1] + record[j][V1];
		}

		swap(V1, V2);
	}

	cout << "Combination Numbers: " << record[b][V1] << endl;
	return record[b][V1];
}


bool mergeSegmentsRecursivly(const DataType* orginal_time_series, const APCA& apca_presentation, const DataType& n, const int& M, APCA& italicC) {

	if (apca_presentation.segmentNum > M) {
		int merge_frequency = apca_presentation.segmentNum - M;
		cout << "merge_frequency: " << merge_frequency << endl;
		int combination_sum = getCombinationNumber(apca_presentation.segmentNum - 1, merge_frequency);
		priority_queue<APCA_MERGE, vector<APCA_MERGE>, cmpMore> queue_reconstrution_error;
		vector<DataType> result_vector;
		APCA_MERGE* best_merge_index;
		APCA_MERGE* p_array_begin_index;
		best_merge_index = new APCA_MERGE[combination_sum];

		combinateSegments(orginal_time_series, apca_presentation, apca_presentation.segmentNum, n, M, merge_frequency, 0, 0, result_vector);

		p_array_begin_index = best_merge_index;

		for (int i = 0; i < combination_sum; i++) {
			best_merge_index[i].r = new int[M];
			best_merge_index[i].v = new DataType[M];
			initialArray(best_merge_index[i], M);
		}

		/*for (vector<DataType>::iterator it = result_vector.begin(); it != result_vector.end(); ++it) {
			cout << *it << " ";
		}
		cout << endl;*/

		int i = 0;
		for (vector<DataType>::iterator it = result_vector.begin(); it != result_vector.end(); ++it) {

			while (i < M) {
				best_merge_index->r[i] = *(it + (i << 1));
				best_merge_index->v[i] = *(it + (i << 1) + 1);
				i++;
			}

			best_merge_index->distance_AE = *(it + (i << 1));
			i = 0;
			queue_reconstrution_error.push(*best_merge_index);
			best_merge_index++;
			it += (M << 1);
		}

		/*while (!queue_reconstrution_error.empty()) {
			for (int i = 0; i < M;i++) {
				cout << queue_reconstrution_error.top().r[i] << " " << queue_reconstrution_error.top().v[i] << endl;
			}
			cout << queue_reconstrution_error.top().distance_AE << endl;
			queue_reconstrution_error.pop();
			cout << endl;
		}*/

		for (int i = 0; i < M; i++) {
			italicC.r[i] = queue_reconstrution_error.top().r[i];
			italicC.v[i] = queue_reconstrution_error.top().v[i];
		}

		getExactMeanValue(orginal_time_series, italicC);

		best_merge_index = p_array_begin_index;
		for (int i = 0; i < combination_sum; i++) {
			//best_merge_index[i].v = nullptr;
			delete[] best_merge_index[i].v;
			//best_merge_index[i].r = nullptr;
			delete[] best_merge_index[i].r;

		}
		//p_array_begin_index = nullptr;
		//best_merge_index = nullptr;
		delete[] best_merge_index;

	}
	else if (apca_presentation.segmentNum == M) {

		for (int i = 0; i < M; i++) {
			italicC.r[i] = apca_presentation.r[i];
			italicC.v[i] = apca_presentation.v[i];
		}

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


void getMergeIndex(const DataType* orginal_time_series, const int& merge_frequency, APCA& apca_presentation) {
	DataType sefments_difference = NULL;
	int merge_start_index = NULL;
	double initial_segment_length = apca_presentation.r[0] + 1.0;
	//cout <<"initial_segment_length:  " <<initial_segment_length << endl;
	for (int i = 0; i < merge_frequency; i++) { //merge times
		sefments_difference = MAXMIUM;

		apca_presentation.segmentNum--;
		cout << "segment Number: " << apca_presentation.segmentNum << endl;
		for (int j = 0; j < apca_presentation.segmentNum; j++) {
			cout << "sefments_difference: " << sefments_difference << " sgetments diff: " << fabs(apca_presentation.v[j + 1] - apca_presentation.v[j]) << endl;
			if (sefments_difference > fabs(apca_presentation.v[j + 1] - apca_presentation.v[j])) {

				sefments_difference = fabs(apca_presentation.v[j + 1] - apca_presentation.v[j]);
				merge_start_index = j;
			}
		}
		cout << "merge_start_index: " << merge_start_index << endl;
		//apca_presentation.v[merge_start_index] = (apca_presentation.v[merge_start_index] + apca_presentation.v[merge_start_index + 1]) / 2;
		//apca_presentation.v[merge_start_index] = getAve(orginal_time_series + int(apca_presentation.r[merge_start_index - 1] + 1), DataType(apca_presentation.r[merge_start_index + 1] - apca_presentation.r[merge_start_index - 1]));
		
		if (merge_start_index != 0) {
			//double temp1 = getAve(orginal_time_series + int(apca_presentation.r[merge_start_index - 1] + 1), DataType(apca_presentation.r[merge_start_index + 1] - apca_presentation.r[merge_start_index - 1]));
			//apca_presentation.v[merge_start_index] = getAve(orginal_time_series + int(apca_presentation.r[merge_start_index - 1] + 1), DataType(apca_presentation.r[merge_start_index + 1] - apca_presentation.r[merge_start_index - 1]));
			cout << apca_presentation.r[merge_start_index] - apca_presentation.r[merge_start_index - 1] << endl;
			cout<<apca_presentation.v[merge_start_index] << endl;

			apca_presentation.v[merge_start_index] = (apca_presentation.v[merge_start_index] * ((apca_presentation.r[merge_start_index] - apca_presentation.r[merge_start_index - 1]) / initial_segment_length) + apca_presentation.v[merge_start_index + 1] * ((apca_presentation.r[merge_start_index + 1] - apca_presentation.r[merge_start_index]) / initial_segment_length)) / ((apca_presentation.r[merge_start_index + 1] - apca_presentation.r[merge_start_index - 1])/ initial_segment_length);
			//cout << temp1 << " *************** " << apca_presentation.v[merge_start_index] << endl;
			//if (temp != apca_presentation.v[merge_start_index]) cout << "ERROROROROROOROROROROROO! " << endl;
		}
		else {
			apca_presentation.v[merge_start_index] = getAve(orginal_time_series, DataType(apca_presentation.r[merge_start_index + 1]));
			double temp = apca_presentation.v[merge_start_index];
			apca_presentation.v[merge_start_index] = (apca_presentation.v[merge_start_index] * ((apca_presentation.r[merge_start_index] + 1) / initial_segment_length) + apca_presentation.v[merge_start_index + 1] * ((apca_presentation.r[merge_start_index + 1] - apca_presentation.r[merge_start_index]) / initial_segment_length)) / (apca_presentation.r[merge_start_index + 1]);
			if (temp != apca_presentation.v[merge_start_index]) cout << "ERROROROROROOROROROROROO!!!!! " << endl;
		}
		
		
		
		
		
		apca_presentation.r[merge_start_index] = apca_presentation.r[merge_start_index + 1];

		merge_start_index++;
		cout << "segment Number: " << apca_presentation.segmentNum << endl;
		cout << "merge_start_index: " << merge_start_index << endl;
		while (merge_start_index < apca_presentation.segmentNum) {
			apca_presentation.r[merge_start_index] = apca_presentation.r[merge_start_index + 1];
			apca_presentation.v[merge_start_index] = apca_presentation.v[merge_start_index + 1];
			merge_start_index++;
		}
		for (int i = 0; i < apca_presentation.segmentNum; i++) {
			cout << apca_presentation.r[i] << " ";
		}
	}
}


bool mergeSegmentsIterator(const DataType* orginal_time_series, APCA& apca_presentation, const DataType& n, const int& M, APCA& italicC) {

	if (apca_presentation.segmentNum > M) {
		int merge_frequency = apca_presentation.segmentNum - M;

		getMergeIndex(orginal_time_series, merge_frequency, apca_presentation);

		cout << "Segment Number: " << apca_presentation.segmentNum << endl;
		for (int i = 0; i < apca_presentation.segmentNum; i++) {
			cout << apca_presentation.r[i] << " " << apca_presentation.v[i] << endl;
		}
		cout << endl;

		memcpy(italicC.r, apca_presentation.r, sizeof(double)*italicC.segmentNum);
		memcpy(italicC.v, apca_presentation.v, sizeof(double)*italicC.segmentNum);

		//getExactMeanValue(orginal_time_series, apca_presentation);
	}
	else if (apca_presentation.segmentNum == M) {

		for (int i = 0; i < M; i++) {
			memcpy(italicC.r, apca_presentation.r, sizeof(double)*italicC.segmentNum);
			memcpy(italicC.v, apca_presentation.v, sizeof(double)*italicC.segmentNum);
		}
		return true;
	}
	else if (apca_presentation.segmentNum < M) {
		cout << "mergeSegmentsIterator(): Something wrong! the number of segments is greater than M!!!" << endl;
		return false;
	}
	else {
		cout << "Something Wrong!!!" << endl;
		return false;
	}
}


void getAPCAPoint(DataType* orginal_time_series, const double &n, const int& N, APCA& italicC) {

	if (n < N) {                                          //n must longer than N
		cout << "ERROR: N has to smaller than n!!!!!!!";
	}

	int i = 0, j = 0, retained_coeffs_length = NULL;
	unsigned int length_power_of_2 = NULL;
	DataType* wavelet_transform_time_series;
	APCA general_apca_presentation;

	priority_queue<DataType> fq_truncate_index;

	length_power_of_2 = getNextPowerOf2(n);
	wavelet_transform_time_series = new DataType[length_power_of_2];

	padZero(orginal_time_series, n, length_power_of_2, wavelet_transform_time_series);

	cout << "Original time series: ";
	for (int i = 0; i < length_power_of_2; i++) {
		cout << wavelet_transform_time_series[i] << " ";
	}
	cout << endl;

	getHDWT(length_power_of_2, N, wavelet_transform_time_series, fq_truncate_index);
	//retained_coeffs_length = fq_truncate_index.top() + 1;
	//cout<< fq_truncate_index.top()<<endl;
	retained_coeffs_length = getNextPowerOf2(fq_truncate_index.top() + 1);
	cout << "retained_coeffs_length: " << retained_coeffs_length << endl;
	general_apca_presentation.segmentNum = retained_coeffs_length;
	general_apca_presentation.v = new DataType[retained_coeffs_length];
	general_apca_presentation.r = new DataType[retained_coeffs_length];

	reconstructApproximateAPCA(wavelet_transform_time_series, length_power_of_2, retained_coeffs_length, fq_truncate_index, general_apca_presentation);
	truncateZero(n, length_power_of_2, general_apca_presentation);
	getExactMeanValue(orginal_time_series, general_apca_presentation);
	mergeSegmentsIterator(orginal_time_series, general_apca_presentation, n, N, italicC);
	//mergeSegmentsRecursivly(orginal_time_series, general_apca_presentation, n, N, italicC);

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
#define ARRAY_LENGTH 8
	DataType test_array[ARRAY_LENGTH] = { 7,5,5,3,3,3,4,6 };
	//DataType test_array[ARRAY_LENGTH] = { 0.80467, 0.36777, 0.24389, 0.026614, -0.2744, 0.096731, -0.74773, -1.6098 };
	//DataType* mp_HDWT_time_series = new DataType[8];
	DataType n = ARRAY_LENGTH, M = 2;
	APCA test_apca;
	test_apca.segmentNum = M;
	test_apca.r = new DataType[M];
	test_apca.v = new DataType[M];
	unsigned int after_padding_length = NULL;
	int retained_coeffs_length = NULL;

	//getHDWT(test_array, getNextPowerOf2(8), M, mp_HDWT_time_series, retained_coeffs_length);
	getAPCAPoint(test_array, n, M, test_apca);
	cout << endl;
	for (int i = 0; i < M; i++) {
		cout << test_apca.r[i] << " " << test_apca.v[i] << endl;
	}
	cout << endl;

	distanceAE(test_array, test_apca);

	//getCombinationNumber(7, 2);
	//test_apca.r = nullptr;
	delete[] test_apca.r;
	//test_apca.v = nullptr;
	delete[] test_apca.v;
	system("pause");
	return 0;
}

