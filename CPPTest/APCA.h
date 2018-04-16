#pragma once
#ifndef APCA_H
#define APCA_H

#include <algorithm>
#include <queue>
#include <string>
#include <cassert>
#include <iostream>
#include <functional> 
#include <numeric>

#include <time.h>
#include <Windows.h> 
#include <stdlib.h>
#include <stdio.h>
#include <list>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstdint>
#include "RTree.h"
#include "saxquantizer.hpp"

using namespace std;
#define APCA_TEMPLATE template<typename DataType>
#define APCA_QUAL  CAPCA<DataType> 
#define APCA_KNN_QUAL CAPCA_KNN<DataType>

typedef double DataType;
typedef double ElementType;
typedef RTree<DataType, ElementType> ApcaRTree;


//DataType static g_index_point_number;
//double* g_query_time_series;


double g_d_time_whole_first_run = NULL;
double g_n_time_getMBR_count = NULL;
double g_n_time_apca_point_count = NULL;
double g_n_time_count1 = NULL;
double g_n_time_count2 = NULL;
double g_n_time_count_while_first_part = NULL;
double g_n_time_count4 = NULL;
double g_n_time_loop_result_count = NULL;
double g_n_time_queue_pop_count = NULL;
double g_n_time_leaf_node_push_count = NULL;
double g_n_time_leaf_node_distanceLB_count = NULL;
double g_n_time_leaf_node_assignment_count = NULL;
double g_n_time_child_node_MINDIST_count = NULL;
double g_n_time_child_node_push_count = NULL;
double g_n_time_base_KNN_push_count = NULL;

int g_n_account_get_result = NULL;
double g_n_account_apca_point = NULL;
int g_n_account_leaf_node = NULL;
int g_n_account_child_node = NULL;
int memory_account = NULL;

template<typename DataType>
class CAPCA
{
public:

	struct APCA;
	struct COEFFICIENT_PAIR;
	struct cmpLess;
	struct APCA_MERGE;
	struct cmpMoreDistance;
	struct cmpLessDistance;
	struct cmpMoreIndex;
	struct APCA_ARRAY;
	

	inline unsigned int getNextPowerOf2(unsigned int array_length);
	bool padZero(DataType* temp_original_array, const DataType& old_length, unsigned int& new_length, DataType* after_padding_array);
	bool truncateZero(const DataType& old_length, const unsigned int& new_length, APCA& apca_presentation);
	void getHDWT(unsigned int& length_power_of_2, const int& M, DataType* wavelet_transform_time_series, priority_queue<DataType>& fq_truncate_index);
	void reconstructApproximateAPCA(DataType* wavelet_transform_time_series, unsigned int& power_of_2, const int& retained_coeffs_length, priority_queue<DataType> fq_truncate_index, APCA& apca_presentation);
	DataType inline getAve(const DataType* originalArray, const DataType &arrayLength);
	bool getExactMeanValue(const DataType* orginal_time_series, APCA& apca_presentation);
	bool mergeSegmentsIterator(APCA& apca_presentation, const DataType& n, const int& M, APCA& italicC);
	void getMergeIndex(const int& merge_frequency, APCA& apca_presentation);
	void getAPCAPoint(DataType* orginal_time_series, const double &n, const int& N, APCA& italicC);

	APCA& getCmax(double* C, const APCA &italicC, APCA &Cmax);
	APCA& getCmin(double* C, const APCA &italicC, APCA &Cmin);
	APCA& getMBR(const APCA &Cmin, const APCA &Cmax, APCA &MBR);

	RTree<DataType, ElementType>& buidRTreeIndex(RTree<DataType, ElementType> &APCARTree, double* g_query_time_series, const double& g_time_series_length, const double& g_index_point_number, APCA_ARRAY *APCALinkOriginal, const string& file_name);

	//friend class CPACA_KNN;
private:

};

template<typename DataType>
class CAPCA_KNN : public CAPCA<DataType>, public RTree<DataType, DataType> {
public:

	struct regionG;
	struct APCA_NODE_PAIR;
	struct ORIGINAL_TIME_SERIES_PAIR;
	struct priorityDecreasing;
	struct priorityIncrement;
	struct priorityDistanceEUC;
	struct incrementCMP;

	double distanceEUC(const double* Q, const DataType& n, const double* C, const DataType& m);
	double distanceLB(const APCA_QUAL::APCA &QProjection, const APCA_QUAL::APCA &italicC);
	void getFileStreamByID(const string& file_name, const double& g_time_series_length, const int& original_time_series_id, DataType* original_time_series);
	regionG& getRegionG(const ApcaRTree::Rect &MBR, regionG &G);
	double minDistQRt(double *Q, const regionG &G, const int &i);
	double MINDISTQR(double *Q, const DataType &n, const regionG &G);
	APCA_QUAL::APCA& QAPCAProjection(double *Q, const double &n, const APCA_QUAL::APCA &italicC, APCA_QUAL::APCA &QProjection);
	bool APCAKNNSearch2(DataType* g_query_time_series, const DataType& g_index_point_number, const DataType& g_time_series_length, const ApcaRTree &apcaRTree, APCA_QUAL::APCA_ARRAY* APCALinkOriginal, const int &K, const string& file_name);

private:


};


APCA_TEMPLATE
struct APCA_QUAL::APCA {
	DataType* v = nullptr;
	DataType* r = nullptr;
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

APCA_TEMPLATE
struct APCA_QUAL::COEFFICIENT_PAIR {
	int HDWT_id = NULL;
	DataType HDWT_original_value = NULL;
	DataType HDWT_normalized_value = NULL;
	DataType HDWT_magnitude = NULL;

};

APCA_TEMPLATE
struct APCA_QUAL::cmpLess {
	bool operator () (const COEFFICIENT_PAIR &a, const COEFFICIENT_PAIR &b) {
		return a.HDWT_magnitude < b.HDWT_magnitude;             // from big to small  
	}
};

APCA_TEMPLATE
struct APCA_QUAL::APCA_MERGE {
	DataType segment_differences = NULL;
	DataType* r = nullptr;
	DataType* v = nullptr;
};

APCA_TEMPLATE
struct APCA_QUAL::cmpMoreDistance {
	bool operator () (const APCA_MERGE &a, const APCA_MERGE &b) {
		return a.segment_differences > b.segment_differences;             // small to big
	}
};

APCA_TEMPLATE
struct APCA_QUAL::cmpLessDistance {
	bool operator () (const APCA_MERGE &a, const APCA_MERGE &b) {
		return a.segment_differences < b.segment_differences;             // big to small
	}
};

APCA_TEMPLATE
struct APCA_QUAL::cmpMoreIndex {
	bool operator () (const APCA_MERGE &a, const APCA_MERGE &b) {
		return *a.r > *b.r;             // small to big
	}
};

APCA_TEMPLATE
struct APCA_QUAL::APCA_ARRAY {
	//DataType *originalLink = nullptr;
	//int original_time_series_id = NULL;
	APCA APCALink;
};


APCA_TEMPLATE
inline unsigned int APCA_QUAL::getNextPowerOf2(unsigned int array_length) {

	if (array_length <= 0) return 0;

	array_length--;
	array_length |= array_length >> 1;
	array_length |= array_length >> 2;
	array_length |= array_length >> 4;
	array_length |= array_length >> 8;
	array_length |= array_length >> 16;

	return array_length + 1;
}

APCA_TEMPLATE
bool APCA_QUAL::padZero(DataType* temp_original_array, const DataType& old_length, unsigned int& new_length, DataType* after_padding_array) {

	int i = 0;
	while (i < int(new_length)) {
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

APCA_TEMPLATE
bool APCA_QUAL::truncateZero(const DataType& old_length, const unsigned int& new_length, APCA& apca_presentation) {
	//cout << "segmentNum: " << apca_presentation.segmentNum;
	//cout << ", new length: " << new_length << ", old length: " << old_length << endl;
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

APCA_TEMPLATE
void APCA_QUAL::getHDWT(unsigned int& length_power_of_2, const int& M, DataType* wavelet_transform_time_series, priority_queue<DataType>& fq_truncate_index) {

	//cout << power_of_2 << endl;
	int f_interation_times = int(log2(length_power_of_2)) - 1;//3
															  //cout << "f_interation_times: " << f_interation_times << endl;
	int temp_power_of_2 = length_power_of_2 >> 1;

	DataType* fd_temp = new DataType[length_power_of_2];
	memory_account = sizeof(fd_temp);
	priority_queue<COEFFICIENT_PAIR, vector<COEFFICIENT_PAIR>, cmpLess> fq_HDWT_coefficients;

	COEFFICIENT_PAIR temp_Coefficient;

	for (int i = 0; i < int(length_power_of_2); i++) {
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

			//cout << "(" << temp_Coefficient.HDWT_original_value << "," << temp_Coefficient.HDWT_magnitude << ") ";
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
		//cout << endl;
	}

	/*cout << "wavelet_transform_time_series: ";
	for (int i = 0; i < length_power_of_2; i++) {
	cout << wavelet_transform_time_series[i] << " ";
	}
	cout << endl;*/

	temp_Coefficient.HDWT_original_value = *wavelet_transform_time_series;
	temp_Coefficient.HDWT_normalized_value = *wavelet_transform_time_series;
	temp_Coefficient.HDWT_magnitude = abs(temp_Coefficient.HDWT_normalized_value);
	temp_Coefficient.HDWT_id = 0;
	fq_HDWT_coefficients.push(temp_Coefficient);

	for (int i = 0; i < M; i++) {
		//cout << " Big magnitude(ID): " << fq_HDWT_coefficients.top().HDWT_id << endl;
		fq_truncate_index.push(fq_HDWT_coefficients.top().HDWT_id);
		fq_HDWT_coefficients.pop();
	}

	/*while (!fq_HDWT_coefficients.empty()) {
	cout << fq_HDWT_coefficients.top().HDWT_id << " " << fq_HDWT_coefficients.top().HDWT_original_value << "    " << fq_HDWT_coefficients.top().HDWT_normalized_value << "    " << fq_HDWT_coefficients.top().HDWT_fabs_value << endl;
	fq_HDWT_coefficients.pop();
	}*/

	delete[] fd_temp;
}

APCA_TEMPLATE
void APCA_QUAL::reconstructApproximateAPCA(DataType* wavelet_transform_time_series, unsigned int& power_of_2, const int& retained_coeffs_length, priority_queue<DataType> fq_truncate_index, APCA& apca_presentation) {

	DataType temp_apca_value = NULL;
	int loop_times = log2(retained_coeffs_length);
	int inner_loop_times = 1;
	queue<DataType> temp_original_queue;

	for (int i = 0; i < retained_coeffs_length; i++) {
		apca_presentation.v[i] = 0.0;
		apca_presentation.r[i] = power_of_2;
		//temp_original_queue.push(wavelet_transform_time_series[i]);
	}

	for (int i = 0; i < int(fq_truncate_index.size()); i++) {
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

	/*cout << " approximate length : " << endl;
	for (int i = 0; i < apca_presentation.segmentNum; i++) {
	cout << apca_presentation.r[i] << " " << apca_presentation.v[i] << endl;
	}*/

}

APCA_TEMPLATE
DataType inline APCA_QUAL::getAve(const DataType* originalArray, const DataType &arrayLength) {
	//cout << "accumulate: " << accumulate(originalArray, originalArray + int(arrayLength), 0) << ", " << arrayLength << endl;
	return accumulate(originalArray, originalArray + int(arrayLength), 0.0) / DataType(arrayLength);
}


APCA_TEMPLATE
bool APCA_QUAL::getExactMeanValue(const DataType* orginal_time_series, APCA& apca_presentation) {

	double fd_segment_length = NULL;
	int i = 0;
	if (apca_presentation.segmentNum == 1) {
		return true;
	}
	//cout << "Exact mean value: " << endl;
	apca_presentation.v[i] = getAve(orginal_time_series, apca_presentation.r[0] + 1);
	//cout << apca_presentation.r[i] << " " << apca_presentation.v[i] << endl;

	for (i = 1; i < apca_presentation.segmentNum; i++) {
		apca_presentation.v[i] = getAve(orginal_time_series + int(apca_presentation.r[i - 1] + 1), double(apca_presentation.r[i] - apca_presentation.r[i - 1]));
		//cout << apca_presentation.r[i] << " " << apca_presentation.v[i] << endl;;
	}
	return true;
}


APCA_TEMPLATE
void APCA_QUAL::getMergeIndex(const int& merge_frequency, APCA& apca_presentation) {
	//void getMergeIndex(const DataType* orginal_time_series, const int& merge_frequency, APCA& apca_presentation) {
	DataType sefments_difference = NULL;
	int merge_start_index = NULL;
	double initial_segment_length = apca_presentation.r[0] + 1.0;

	for (int i = 0; i < merge_frequency; i++) { //merge times
		sefments_difference = DBL_MAX;

		apca_presentation.segmentNum--;
		//cout << "segment Number: " << apca_presentation.segmentNum << endl;
		for (int j = 0; j < apca_presentation.segmentNum; j++) {
			//cout << "sefments_difference: " << sefments_difference << " sgetments diff: " << fabs(apca_presentation.v[j + 1] - apca_presentation.v[j]) << endl;
			if (sefments_difference > fabs(apca_presentation.v[j + 1] - apca_presentation.v[j])) {
				sefments_difference = fabs(apca_presentation.v[j + 1] - apca_presentation.v[j]);
				merge_start_index = j;
			}
		}
		//cout << "merge_start_index: " << merge_start_index << endl;
		//apca_presentation.v[merge_start_index] = (apca_presentation.v[merge_start_index] + apca_presentation.v[merge_start_index + 1]) / 2;
		if (merge_start_index != 0) {

			//double temp1 = getAve(orginal_time_series + int(apca_presentation.r[merge_start_index - 1] + 1), DataType(apca_presentation.r[merge_start_index + 1] - apca_presentation.r[merge_start_index - 1]));
			//apca_presentation.v[merge_start_index] = getAve(orginal_time_series + int(apca_presentation.r[merge_start_index - 1] + 1), DataType(apca_presentation.r[merge_start_index + 1] - apca_presentation.r[merge_start_index - 1]));

			apca_presentation.v[merge_start_index] = (apca_presentation.v[merge_start_index] * (apca_presentation.r[merge_start_index] - apca_presentation.r[merge_start_index - 1]) + apca_presentation.v[merge_start_index + 1] * (apca_presentation.r[merge_start_index + 1] - apca_presentation.r[merge_start_index])) / (apca_presentation.r[merge_start_index + 1] - apca_presentation.r[merge_start_index - 1]);

			/*if (temp1 != apca_presentation.v[merge_start_index]) {
			cout << "ERROROROROROOROROROROROO! " << endl;
			cout << temp1 <<"  "<<apca_presentation.v[merge_start_index] << endl;
			}*/
		}
		else {
			//apca_presentation.v[merge_start_index] = getAve(orginal_time_series, DataType(apca_presentation.r[merge_start_index + 1]));
			/*double temp2 = apca_presentation.v[merge_start_index];
			double temp4 = apca_presentation.r[merge_start_index];
			double temp3 = apca_presentation.v[merge_start_index+1];
			double temp5 = apca_presentation.r[merge_start_index + 1];*/
			//double temp1 = getAve(orginal_time_series, DataType(apca_presentation.r[merge_start_index + 1]+1));
			apca_presentation.v[merge_start_index] = (apca_presentation.v[merge_start_index] * double(apca_presentation.r[merge_start_index] + 1.0) + apca_presentation.v[merge_start_index + 1] * (apca_presentation.r[merge_start_index + 1] - apca_presentation.r[merge_start_index])) / double(apca_presentation.r[merge_start_index + 1] + 1);
			/*if (temp1 != apca_presentation.v[merge_start_index]) {
			cout << "ERROROROROROOROROROROROO! " << endl;
			cout << "initial_segment_length: " << initial_segment_length << endl;
			cout<<temp4 <<" "<<temp2 << endl;

			cout <<temp5 <<" "<<temp3 << endl;
			cout << merge_start_index << endl;
			cout << temp1 <<"    "<< apca_presentation.v[merge_start_index] << endl;
			}*/
		}

		apca_presentation.r[merge_start_index] = apca_presentation.r[merge_start_index + 1];

		merge_start_index++;
		//cout << "segment Number: " << apca_presentation.segmentNum << endl;
		//cout << "merge_start_index: " << merge_start_index << endl;
		while (merge_start_index < apca_presentation.segmentNum) {
			apca_presentation.r[merge_start_index] = apca_presentation.r[merge_start_index + 1];
			apca_presentation.v[merge_start_index] = apca_presentation.v[merge_start_index + 1];
			merge_start_index++;
		}
		/*for (int i = 0; i < apca_presentation.segmentNum; i++) {
		cout << apca_presentation.r[i] << " ";
		}*/
	}
}

APCA_TEMPLATE
bool APCA_QUAL::mergeSegmentsIterator(APCA& apca_presentation, const DataType& n, const int& M, APCA& italicC) {
	//bool mergeSegmentsIterator(const DataType* orginal_time_series, APCA& apca_presentation, const DataType& n, const int& M, APCA& italicC) {

	if (apca_presentation.segmentNum > M) {
		int merge_frequency = apca_presentation.segmentNum - M;

		//getMergeIndexQueue(merge_frequency, apca_presentation); //2 for 1
		getMergeIndex(merge_frequency, apca_presentation); //2 for 1

														   //cout << "Segment Number: " << apca_presentation.segmentNum << endl;
														   /*for (int i = 0; i < apca_presentation.segmentNum; i++) {
														   cout << apca_presentation.r[i] << " " << apca_presentation.v[i] << endl;
														   }
														   cout << endl;*/

														   //memcpy(italicC.r, apca_presentation.r, sizeof(double)*italicC.segmentNum);
														   //memcpy(italicC.v, apca_presentation.v, sizeof(double)*italicC.segmentNum);
		assert(apca_presentation.segmentNum == italicC.segmentNum);
		copy(apca_presentation.r, apca_presentation.r + apca_presentation.segmentNum, italicC.r);
		copy(apca_presentation.v, apca_presentation.v + apca_presentation.segmentNum, italicC.v);

		return true;
		/*for (int i = 0; i < italicC.segmentNum; i++) {
		cout << italicC.r[i] << " " << italicC.v[i] << endl;
		}
		cout << endl;*/

		//getExactMeanValue(orginal_time_series, apca_presentation);
	}
	else if (apca_presentation.segmentNum == M) {
		//memcpy(italicC.r, apca_presentation.r, sizeof(double)*italicC.segmentNum);
		//memcpy(italicC.v, apca_presentation.v, sizeof(double)*italicC.segmentNum);

		copy(apca_presentation.r, apca_presentation.r + apca_presentation.segmentNum, italicC.r);
		copy(apca_presentation.v, apca_presentation.v + apca_presentation.segmentNum, italicC.v);

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

APCA_TEMPLATE
void APCA_QUAL::getAPCAPoint(DataType* orginal_time_series, const double &n, const int& N, APCA& italicC) {

	assert(n > N);

	int i = 0, j = 0, retained_coeffs_length = NULL;
	unsigned int length_power_of_2 = NULL;
	DataType* wavelet_transform_time_series;
	APCA general_apca_presentation;

	priority_queue<DataType> fq_truncate_index;

	length_power_of_2 = getNextPowerOf2(n);
	wavelet_transform_time_series = new DataType[length_power_of_2];

	padZero(orginal_time_series, n, length_power_of_2, wavelet_transform_time_series);

	/*cout << "Add zero Original time series: ";
	for (int i = 0; i < length_power_of_2; i++) {
	cout << wavelet_transform_time_series[i] << " ";
	}
	cout << endl;*/

	getHDWT(length_power_of_2, N, wavelet_transform_time_series, fq_truncate_index);
	//retained_coeffs_length = fq_truncate_index.top() + 1;
	//cout<< fq_truncate_index.top()<<endl;
	retained_coeffs_length = getNextPowerOf2(fq_truncate_index.top() + 1);
	//cout << "retained_coeffs_length: " << retained_coeffs_length << endl;
	general_apca_presentation.segmentNum = retained_coeffs_length;
	general_apca_presentation.v = new DataType[retained_coeffs_length];
	general_apca_presentation.r = new DataType[retained_coeffs_length];

	reconstructApproximateAPCA(wavelet_transform_time_series, length_power_of_2, retained_coeffs_length, fq_truncate_index, general_apca_presentation);
	//cout << "Total Resolution: " << log2(getNextPowerOf2(n)) << ",  Resolution: " << log2(retained_coeffs_length) << endl;

	truncateZero(n, length_power_of_2, general_apca_presentation);
	getExactMeanValue(orginal_time_series, general_apca_presentation);
	mergeSegmentsIterator(general_apca_presentation, n, N, italicC);
	//mergeSegmentsRecursivly0(orginal_time_series, general_apca_presentation, n, N, italicC);

	delete[] general_apca_presentation.r;
	delete[] general_apca_presentation.v;
	delete[] wavelet_transform_time_series;
}

APCA_TEMPLATE
typename APCA_QUAL::APCA &APCA_QUAL::getCmin(double* C, const APCA &italicC, APCA &Cmin) {
	//cout << "getCmin() \n";

	Cmin.segmentNum = italicC.segmentNum;
	//cout << "segmentNum = " << Cmin.segmentNum << endl;

	int i = 0;
	Cmin.r[i] = italicC.r[i];
	Cmin.v[i] = *min_element(C, C + int(italicC.r[i]) + 1);

	//cout << "r[" << i << "] = " << Cmin->r[i] << ", Cmin0[" << i << "] = " << Cmin->v[i] << endl << endl;

	for (i = 1; i < italicC.segmentNum; i++) {
		Cmin.r[i] = italicC.r[i];
		Cmin.v[i] = *min_element(C + int(italicC.r[i - 1]) + 1, C + int(italicC.r[i]) + 1);

		//cout << "r[" << i << "] = " << Cmin.r[i] << ", Cmin0[" << i << "] = " << Cmin.v[i] << endl << endl;
	}
	return Cmin;
}

APCA_TEMPLATE
typename APCA_QUAL::APCA& APCA_QUAL::getCmax(double* C, const APCA &italicC, APCA &Cmax) {
	//cout << "getCmax() \n";
	Cmax.segmentNum = italicC.segmentNum;
	//cout << "segmentNum = " << Cmax.segmentNum << endl;
	int i = 0;
	Cmax.r[i] = italicC.r[i];
	Cmax.v[i] = *max_element(C, C + int(italicC.r[i]) + 1);

	//cout << "r[" << i << "] = " << Cmax.r[i] << ", Cmax0[" << i << "] = " << Cmax.v[i] << endl << endl;

	for (i = 1; i < italicC.segmentNum; i++) {
		Cmax.r[i] = italicC.r[i];
		Cmax.v[i] = *max_element(C + int(italicC.r[i - 1]) + 1, C + int(italicC.r[i]) + 1);

		//cout << "r[" << i << "] = " << Cmax.r[i] << ", Cmax0[" << i << "] = " << Cmax.v[i] << endl << endl;
	}
	return Cmax;
}

APCA_TEMPLATE
typename APCA_QUAL::APCA& APCA_QUAL::getMBR(const APCA &Cmin, const APCA &Cmax, APCA &MBR) {//N=2M
															 //cout << "getMBR()\n";
	MBR.segmentNum = Cmin.segmentNum * 2;

	if (Cmin.segmentNum != Cmax.segmentNum) {
		cout << "Wrong!!!!!!!!!!!" << endl;
	}

	int i = NULL;
	for (i = 0; i < MBR.segmentNum; i++) {
		if (i & 1) {							//odd odd is Cv
			MBR.r[i] = Cmin.v[(i - 1) / 2]; //r is min value
			MBR.v[i] = Cmax.v[(i - 1) / 2];  //v is max value
		}
		else {									//even is Cr
			MBR.r[i] = Cmin.r[i / 2];
			MBR.v[i] = Cmax.r[i / 2];
		}
		//cout << "L[" << i << "] = " << MBR.r[i] << ", H[" << i << "] = " << MBR.v[i] << endl;
	}
	return MBR;
}



APCA_TEMPLATE
RTree<DataType, ElementType>& APCA_QUAL::buidRTreeIndex(RTree<DataType, ElementType> &APCARTree, double* g_query_time_series, const double& g_time_series_length, const double& g_index_point_number, APCA_ARRAY *APCALinkOriginal, const string& file_name) {
	//RTree<DataType, ElementType>& buidRTreeIndex(RTree<DataType, ElementType> &APCARTree, double* g_query_time_series, const double& g_time_series_length, const double& g_index_point_number, double(&test_d_original_time_series)[ROW][COLUMN], Link *APCALinkOriginal) {
	printf(">>>***Build RTree Index***<<<\n");

	int i = NULL, j = NULL, f_insert_count = NULL;
	DataType* original_time_series = new DataType[g_time_series_length + 1];
	APCA CminParameter;//temp  Cmin
	CminParameter.r = new double[APCARTree.NUMDIMS / 2];
	CminParameter.v = new double[APCARTree.NUMDIMS / 2];
	APCA CmaxParameter;//temp Cmax
	CmaxParameter.r = new double[APCARTree.NUMDIMS / 2];
	CmaxParameter.v = new double[APCARTree.NUMDIMS / 2];
	APCA MBRParameter; //Temp MBR
	MBRParameter.r = new double[APCARTree.NUMDIMS];
	MBRParameter.v = new double[APCARTree.NUMDIMS];

	memory_account = sizeof(*CminParameter.r);
	cout<<"Memory account: "<< memory_account<<endl;
	/*cout << "Query Point : ";
	for (i = 0; i < g_time_series_length; i++) cout << g_query_time_series[i] << ", ";
	cout << endl;*/

	string fs_row_string;
	string fs_row_number;
	ifstream file_stream = ifstream(file_name);
	assert(file_stream);

	f_insert_count = 0;
	for (i = 0; i < g_index_point_number; i++) {

		//getRandomAPCAPoint(APCALinkOriginal[i].originalLink, g_time_series_length);
		//APCALinkOriginal[i].originalLink = test_d_original_time_series[i];
		//getNormalArray(APCALinkOriginal[i].originalLink, g_time_series_length);

		/*cout << "Original Point : ";
		for (j = 0; j < g_time_series_length; j++) {
		cout << APCALinkOriginal[i].originalLink[j] << ", ";
		}
		cout << endl;*/

		file_stream >> fs_row_string;
		stringstream sstr(fs_row_string);
		int string_id = 0;
		while (getline(sstr, fs_row_number, ',') && string_id < g_time_series_length + 1) {
			original_time_series[string_id] = stod(fs_row_number);
			//cout << original_time_series[string_id] << ", ";
			string_id++;
		}
		//cout << endl;

		for (int j = 0; j < g_time_series_length; j++) {
			original_time_series[j] = original_time_series[j + 1];
		}

		//divideRemainderPAA(original_time_series, APCALinkOriginal[i].APCALink, g_time_series_length, APCARTree.NUMDIMS / 2);
		getAPCAPoint(original_time_series, g_time_series_length, APCARTree.NUMDIMS / 2, APCALinkOriginal[i].APCALink);
		//divideRemainderPAA(APCALinkOriginal[i].originalLink, APCALinkOriginal[i].APCALink, g_time_series_length, APCARTree.NUMDIMS / 2);
		//getAPCAPoint(APCALinkOriginal[i].originalLink, g_time_series_length, APCARTree.NUMDIMS / 2, APCALinkOriginal[i].APCALink);

		/*cout << "APCA Point : ";
		for (j = 0; j < APCARTree.NUMDIMS / 2; j++) {
		cout << "(" << APCALinkOriginal[i].APCALink.r[j] << ", " << APCALinkOriginal[i].APCALink.v[j] << ")" << " ";
		}
		cout << endl;*/

		getMBR(getCmin(original_time_series, APCALinkOriginal[i].APCALink, CminParameter), getCmax(original_time_series, APCALinkOriginal[i].APCALink, CmaxParameter), MBRParameter);
		//getMBR(getCmin(APCALinkOriginal[i].originalLink, APCALinkOriginal[i].APCALink, CminParameter), getCmax(APCALinkOriginal[i].originalLink, APCALinkOriginal[i].APCALink, CmaxParameter), MBRParameter);

		/*cout << "MBR: ";
		for (j = 0; j < MBRParameter.segmentNum; j++) {
		cout << " (" << MBRParameter.r[j] << ", " << MBRParameter.v[j] << ") ";
		}
		cout << "\n\nBegin Insert " << f_insert_count << ": \n";*/

		APCARTree.Insert(MBRParameter.r, MBRParameter.v, i);

		f_insert_count++;
	}
	file_stream.close();
	cout << "Root Node : sub node number = " << APCARTree.m_root->m_count << " Root level = : " << APCARTree.m_root->m_level << "\n\nBegin to build a RTree:\n";
	cout << "\n RTree conclusion\n The number of RTree Data Point = : " << APCARTree.Count() << endl;

	delete[] original_time_series;
	delete[] CminParameter.r;
	delete[] CminParameter.v;
	delete[] CmaxParameter.r;
	delete[] CmaxParameter.v;
	delete[] MBRParameter.r;
	delete[] MBRParameter.v;

	return APCARTree;
}

APCA_TEMPLATE
struct APCA_KNN_QUAL::regionG {
	int regionNum = NULL;
	double *G1 = nullptr;
	double *G2 = nullptr;
	double *G3 = nullptr;
	double *G4 = nullptr;
};

APCA_TEMPLATE
struct APCA_KNN_QUAL::APCA_NODE_PAIR {								//queue data structure.
	double d_dist = NULL;										//dist.
	int original_time_series_id = NULL;                  // APCA point ID.
	APCA_QUAL::APCA *p_APCA_point = nullptr;						//data point. if NULL, this is internal node.
	RTree<DataType, DataType>::Node *p_rtree_node = nullptr; //subNode, if NULL, this is apca point.
};

APCA_TEMPLATE
struct APCA_KNN_QUAL::ORIGINAL_TIME_SERIES_PAIR {  //for temp queue
	double d_dist = NULL;
	double *p_original_time_series = nullptr;
	int original_time_series_id = NULL;
	/*~originalTimeSeriesPair(){
	if (value!=nullptr)
	{
	delete[] value;
	value = NULL;
	}
	}*/
};

APCA_TEMPLATE
struct APCA_KNN_QUAL::priorityDecreasing {

	bool operator ()(const APCA_NODE_PAIR & a, const APCA_NODE_PAIR & b) {
		return a.d_dist < b.d_dist;
	}
};

APCA_TEMPLATE
struct APCA_KNN_QUAL::priorityIncrement {

	bool operator ()(const APCA_NODE_PAIR & a, const APCA_NODE_PAIR & b) {
		return a.d_dist > b.d_dist;
	}
};



APCA_TEMPLATE
struct APCA_KNN_QUAL::priorityDistanceEUC {

	bool operator ()(const ORIGINAL_TIME_SERIES_PAIR & a, const ORIGINAL_TIME_SERIES_PAIR & b) {
		return a.d_dist > b.d_dist;
	}
};


APCA_TEMPLATE
struct APCA_KNN_QUAL::incrementCMP {

	bool operator ()(const double & a, const double & b) {
		return a > b;
	}
};

APCA_TEMPLATE
double APCA_KNN_QUAL::distanceEUC(const double* Q, const DataType& n, const double* C, const DataType& m) {
	//cout << "\n\n distanceEUC()" << endl;
	int i = 0, j = 0;
	double distance = 0, sum = 0;

	assert(n == m);

	//sum += italicC.r[0] * pow((QProjection.v[0] - italicC.v[0]), 2);

	for (i = 0; i < n; i++) {
		//cout << Q[i]<<" " << C[i]<<" " << Q[i] - C[i] << endl;
		sum += pow((Q[i] - C[i]), 2);
	}
	distance = sqrt(sum);
	//cout << "distanceEUC  : " << distance << endl;
	return distance;
}




APCA_TEMPLATE
double APCA_KNN_QUAL::distanceLB(const APCA_QUAL::APCA &QProjection, const APCA_QUAL::APCA &italicC) {

	//cout << QProjection.v[0]<<"   #####" << italicC.v[0] << endl;
	//cout << "\n\n distanceLB() \n";
	int i = 0, j = 0;
	double distance = 0;

	assert(QProjection.segmentNum == italicC.segmentNum);

	double sum = (italicC.r[0] + 1) * (QProjection.v[0] - italicC.v[0])*(QProjection.v[0] - italicC.v[0]);
	//cout << QProjection.v[i] << " " << italicC.v[i] << " " << QProjection.v[i]+1 << endl;
	//cout << sum << endl;

	for (i = 1; i < italicC.segmentNum; i++) {
		//cout<<i<<" "<< italicC.r[i] - italicC.r[i - 1]<<" " << QProjection.v[i]<<" " << italicC.v[i] <<" "<< QProjection.v[i] - italicC.v[i] << endl;
		sum += (italicC.r[i] - italicC.r[i - 1]) * (QProjection.v[i] - italicC.v[i])*(QProjection.v[i] - italicC.v[i]);
		//cout << sum << endl;
	}

	distance = sqrt(sum);
	//cout << "distanceLB  : " << distance << endl;
	return distance;
}

APCA_TEMPLATE
void APCA_KNN_QUAL::getFileStreamByID(const string& file_name, const double& g_time_series_length, const int& original_time_series_id, DataType* original_time_series) {
	//bool getFileStream(const string& file_name, T(&mg_input_data)[ROW][COLUMN]) {

	/*T* mp_pointer = &mg_input_data[0][0];*/
	string fs_row_string;
	string fs_row_number;

	/*for (int i = 0; i < ROW; i++) {
	for (int j = 0; j <= COLUMN; j++) {
	mg_input_data[i][j] = NULL;
	}
	}*/

	ifstream file_stream = ifstream(file_name);
	assert(file_stream);
	/*if (!file_stream) {
	cout << "Can't open file " << file_name << " !!!!";
	return false;
	}*/
	int i = 0;
	while (!file_stream.eof())
	{
		file_stream >> fs_row_string;

		if (i == original_time_series_id) {
			stringstream sstr(fs_row_string);
			int string_id = 0;
			while (getline(sstr, fs_row_number, ',') && string_id < g_time_series_length) {
				original_time_series[string_id] = stod(fs_row_number);
				string_id++;
			}
			break;
		}
		i++;
	}

	file_stream.close();

	/*mp_pointer = &mg_input_data[0][0];
	for (int i = 0; i < ROW; i++) {
	for (int j = 0; j < COLUMN; j++) {
	cout << mg_input_data[i][j] << " ";
	}
	cout << endl << endl;
	}*/

}

APCA_TEMPLATE
typename APCA_KNN_QUAL::regionG& APCA_KNN_QUAL::getRegionG(const ApcaRTree::Rect &MBR, regionG &G) {

	double run_time0;
	_LARGE_INTEGER time_start0;  //start time 
	_LARGE_INTEGER time_over0;   //finish time  
	double dqFreq0;      //timer frequency  
	LARGE_INTEGER f0;    //timer frequency  
	QueryPerformanceFrequency(&f0);
	dqFreq0 = (double)f0.QuadPart;
	QueryPerformanceCounter(&time_start0);

	int i = 0;
	//cout << "getRegionG()\n";

	/*if (!NUMDIMS % 2) {
	cout << "NUMDIMS is not even    Wrong!!!!!";
	}*/
	G.G1[i] = MBR.m_min[1];
	G.G2[i] = 0;
	G.G3[i] = MBR.m_max[1];
	G.G4[i] = MBR.m_max[0];
	//cout << "G1[" << i << "] = " << G.G1[i] << ", G2[" << i << "] = " << G.G2[i] << ", G3[" << i << "] = " << G.G3[i] << ", G4[" << i << "] = " << G.G4[i] << endl;

	for (i = 1; i < G.regionNum; i++) {
		G.G1[i] = MBR.m_min[i * 2 + 1];
		G.G2[i] = MBR.m_min[i * 2 - 2] + 1;
		G.G3[i] = MBR.m_max[i * 2 + 1];
		G.G4[i] = MBR.m_max[i * 2];
		//cout << "G1[" << i << "] = " << G.G1[i] << ", G2[" << i << "] = " << G.G2[i] << ", G3[" << i << "] = " << G.G3[i] << ", G4[" << i << "] = " << G.G4[i] << endl;
	}

	QueryPerformanceCounter(&time_over0);    //Finish recording time  
	run_time0 = 1000000 * (time_over0.QuadPart - time_start0.QuadPart) / dqFreq0;
	g_n_time_getMBR_count += run_time0;

	//printf("\nrun_time = %f us\n", run_time0);
	return G;
}


APCA_TEMPLATE
double APCA_KNN_QUAL::minDistQRt(double *Q, const regionG &G, const int &i) {

	int j = 0, count = 0;
	double minDistanceQG = 0;
	//cout << "\nG.regionNum = " << G.regionNum << ", q[" << i << "] =" << Q[i] << endl;

	priority_queue<double, vector<double>, incrementCMP > tempMinDist;
	//double* f_d_temp_result_array = new double[G.regionNum];

	for (j = 0; j < G.regionNum; j++) {
		//cout << "G2[" << j << "] = " << G.G2[j] << ", G4[" << j << "] = " << G.G4[j] << ", G1[" << j << "] = " << G.G1[j] << ", G3[" << j << "] = " << G.G3[j] << endl;
		if (i >= G.G2[j] && i <= G.G4[j]) {
			//cout << "Active!" << endl;
			if (Q[i] < G.G1[j]) {
				minDistanceQG = pow((G.G1[j] - Q[i]), 2);

				//minDistanceArray[count] = minDistanceQG;
				//count++;
				tempMinDist.push(minDistanceQG);
				//cout << "q[" << i << "] < G1[" << j << "], MINDIST = " << minDistanceQG << "\n\n";
			}
			else if (Q[i] > G.G3[j]) {
				minDistanceQG = pow((Q[i] - G.G3[j]), 2);
				//minDistanceArray[count] = minDistanceQG;
				//count++;
				tempMinDist.push(minDistanceQG);
				//cout << "q[" << i << "] > G3[" << j << "], MINDIST = " << minDistanceQG << "\n\n";
			}
			else {
				minDistanceQG = 0;
				//cout << "G1 >= qt <= G3, MINDIST = " << minDistanceQG << endl;
				//delete[] f_d_temp_result_array;
				return minDistanceQG;
			}
		}
		else {
			//cout << "G is not active!!!!! \n\n";
			continue;
		}
	}
	//minDistanceQG = *min_element(minDistanceArray, minDistanceArray + count);
	minDistanceQG = tempMinDist.top();

	//while (!tempMinDist.empty()) {
	//tempMinDist.pop();
	//}
	//cout << endl;
	//cout << "The MINDIST = " << minDistanceQG << endl;
	//queue <regionG> activeG;
	//delete[] f_d_temp_result_array;
	return minDistanceQG;
}




APCA_TEMPLATE
double APCA_KNN_QUAL::MINDISTQR(double *Q, const DataType &n, const regionG &G) {
	//cout << "MINDISTQR()\n";
	int i = NULL;
	double minDistanceQR = 0;
	for (i = 0; i < n; i++) {
		minDistanceQR += minDistQRt(Q, G, i);
	}
	//queue <regionG> activeG;
	//cout << "SUM MINDIST(Q,R) = " << sqrt(minDistanceQR) << endl;
	return sqrt(minDistanceQR);
}

APCA_TEMPLATE
typename APCA_QUAL::APCA& APCA_KNN_QUAL::QAPCAProjection(double *Q, const double &n, const APCA_QUAL::APCA &italicC, APCA_QUAL::APCA &QProjection) {
	//cout << "QAPCAProjection()" << endl;
	int i = 0, j = 0;
	double sum = 0;
	QProjection.segmentNum = italicC.segmentNum;
	//cout << "\n" << "Q.segmentNum = " << QProjection.segmentNum << endl;

	if (n != italicC.r[italicC.segmentNum - 1] + 1) {
		cout << "QAPCAProjection()  :  !!!!ERROR:: Q[].segmentNum == C.segmentNum!!!!!!!" << endl;
	}

	QProjection.r[0] = italicC.r[0];
	while (j <= italicC.r[0]) {
		sum += Q[j];
		//cout << "Q[" << j << "] = " << Q[j] << " ";
		j++;
	}
	QProjection.v[i] = sum / (italicC.r[i] + 1);
	//cout << "\nQv[" << i << "] = " << QProjection.v[i] << ", Qr[" << i << "] = " << QProjection.r[i] << "\n\n";

	for (i = 1; i < italicC.segmentNum; i++) {
		QProjection.r[i] = italicC.r[i];
		sum = 0;
		for (j = 1; j <= (italicC.r[i] - italicC.r[i - 1]); j++) {
			sum += Q[j + int(italicC.r[i - 1])];
			//cout << "Q[" << j + italicC.r[i - 1] << "] = " << Q[j + int(italicC.r[i - 1])] << " ";
		}
		QProjection.v[i] = sum / (italicC.r[i] - italicC.r[i - 1]);
		//cout << "\nQv[" << i << "] = " << QProjection.v[i] << ", Qr[" << i << "] = " << QProjection.r[i] << "\n\n";
	}
	return QProjection;
}
APCA_TEMPLATE
bool APCA_KNN_QUAL::APCAKNNSearch2(DataType* g_query_time_series, const DataType& g_index_point_number, const DataType& g_time_series_length, const ApcaRTree &apcaRTree, APCA_QUAL::APCA_ARRAY* APCALinkOriginal, const int &K, const string& file_name) {

	double whole_first_run_time;
	_LARGE_INTEGER whole_first_time_start;  //start time 
	_LARGE_INTEGER whole_first_time_over;   //finish time  
	double whole_first_dqFreq;      //timer frequency  
	LARGE_INTEGER whole_first_f;    //timer frequency  
	QueryPerformanceFrequency(&whole_first_f);
	whole_first_dqFreq = (double)whole_first_f.QuadPart;
	QueryPerformanceCounter(&whole_first_time_start);

	//printf("<///////**  APCA KNN Begin  **////////>\n");
	int i = NULL, j = NULL;

	assert(K < g_index_point_number);

	priority_queue<APCA_NODE_PAIR, vector<APCA_NODE_PAIR>, priorityIncrement > queue;
	APCA_NODE_PAIR f_APCA_Root, f_temp_APCA_Pair;

	list <ORIGINAL_TIME_SERIES_PAIR> temp, result;
	ORIGINAL_TIME_SERIES_PAIR tempOriginalTimeSeriesPair;

	f_APCA_Root.p_rtree_node = apcaRTree.m_root;
	f_APCA_Root.d_dist = 0;
	queue.push(f_APCA_Root);
	//cout << "Queue.top = " << queue.top().key << " " << queue.size() << " " << queue.top().APCAValue << " " << queue.top().id_originalTimeSeries << " " << queue.top().value << endl;
	//printf("<///////**    KNN Begin   **////////>\n");

	int n_data_point_count = 0;
	int n_distanceLB_count = 0;
	int n_push_leaf_node_count = 0;

	g_n_time_getMBR_count = 0;
	g_n_time_apca_point_count = 0;
	g_n_time_count1 = 0;
	g_n_time_count2 = 0;
	g_n_time_count_while_first_part = 0;
	g_n_time_count4 = 0;
	g_n_time_loop_result_count = 0;
	g_n_time_leaf_node_push_count = 0;
	g_n_time_leaf_node_distanceLB_count = 0;
	g_n_time_leaf_node_assignment_count = 0;
	g_n_time_child_node_MINDIST_count = 0;
	g_n_time_child_node_push_count = 0;


	QueryPerformanceCounter(&whole_first_time_over);    //Finish recording time  
	g_d_time_whole_first_run = whole_first_run_time = 1000000 * (whole_first_time_over.QuadPart - whole_first_time_start.QuadPart) / whole_first_dqFreq;

	int n_result_count = 0;
	while (!queue.empty()) {
		/*while (!queue.empty()|| result.size() != K) {*/
		//while (result.size() != K) {
		//cout << "KNN: the " << count << " turn : \n";

		double run_time3;
		_LARGE_INTEGER time_start3;  //start time 
		_LARGE_INTEGER time_over3;   //finish time  
		double dqFreq3;      //timer frequency  
		LARGE_INTEGER f3;    //timer frequency  
		QueryPerformanceFrequency(&f3);
		dqFreq3 = (double)f3.QuadPart;
		QueryPerformanceCounter(&time_start3);

		double run_time4;
		_LARGE_INTEGER time_start4;  //start time 
		_LARGE_INTEGER time_over4;   //finish time  
		double dqFreq4;      //timer frequency  
		LARGE_INTEGER f4;    //timer frequency  
		QueryPerformanceFrequency(&f4);
		dqFreq4 = (double)f4.QuadPart;
		QueryPerformanceCounter(&time_start4);

		//while (queue.empty() && result.size() != K) {
		//	//printf("Queue is Empty! But K is not enough!!!\n");
		//	result.push(temp.top());
		//	//cout << "Result = " << result.size() << endl;
		//	temp.pop();
		//}

		QueryPerformanceCounter(&time_over4);    //Finish recording time  
		run_time4 = 1000000 * (time_over4.QuadPart - time_start4.QuadPart) / dqFreq4;
		g_n_time_count4 += run_time4;
		//printf("\nrun_time = %f us\n", run_time0);

		/*if (count == K) {
		cout << "////*********EMPTY  Found " << K << " result!!!!!!********///" << endl;break;}*/

		double run_time5;
		_LARGE_INTEGER time_start5;  //start time 
		_LARGE_INTEGER time_over5;   //finish time  
		double dqFreq5;      //timer frequency  
		LARGE_INTEGER f5;    //timer frequency  
		QueryPerformanceFrequency(&f5);
		dqFreq5 = (double)f5.QuadPart;
		QueryPerformanceCounter(&time_start5);

		APCA_NODE_PAIR m_temp_queue_top = queue.top();

		cout << "    Begin Loop:     top.dist: " << m_temp_queue_top.d_dist << "    temp.size() = " << temp.size() << ", temp iterator: ";
		for (list<ORIGINAL_TIME_SERIES_PAIR>::iterator it = temp.begin(); it != temp.end(); ++it) cout << it->d_dist << ", ";
		cout << endl;

		for (list<ORIGINAL_TIME_SERIES_PAIR>::iterator plist = temp.begin(); plist != temp.end(); ) {
			//cout << "        Loop: " << plist->d_dist << " vs " << m_temp_queue_top.d_dist << endl;
			if (plist->d_dist <= m_temp_queue_top.d_dist) {
				cout << "           <= " << endl;
				result.push_back(*plist);
				plist = temp.erase(plist);
			}
			else plist++;
			if (K == result.size()) {
				cout << "*****Find result!!!!!!!!!!!!! " << result.size() << endl;
				cout << "Result list: ";
				for (list<ORIGINAL_TIME_SERIES_PAIR>::iterator it = result.begin(); it != result.end(); ++it)
					cout << it->d_dist << ", " << it->original_time_series_id << ", ";

				return true;
			}
		}
		//cout << "**count = " << count << endl;

		QueryPerformanceCounter(&time_over5);    //Finish recording time  
		run_time5 = 1000000 * (time_over5.QuadPart - time_start5.QuadPart) / dqFreq5;
		g_n_time_loop_result_count += run_time5;

		//printf("\nrun_time = %f us\n", run_time0);

		cout << "Before Pop: Queue.size(): " << queue.size() << " Pop: top.dist:" << queue.top().d_dist << " " << queue.top().original_time_series_id << endl;
		queue.pop();
		cout << "After Pop: Queue.size(): " << queue.size();
		if (!queue.empty()) cout << " Pop: top.dist:" << queue.top().d_dist << " " << queue.top().original_time_series_id << endl;
		else cout << endl;

		QueryPerformanceCounter(&time_over3);    //Finish recording time  
		run_time3 = 1000000 * (time_over3.QuadPart - time_start3.QuadPart) / dqFreq3;
		g_n_time_count_while_first_part += run_time3;

		//printf("\nrun_time = %f us\n", run_time0);

		if (m_temp_queue_top.p_rtree_node == nullptr) { //is APCA data point
			cout << "    queue.top is data point\n";

			DataType* original_time_series = new DataType[g_time_series_length + 1];

			g_n_account_apca_point++;
			double run_time0;
			_LARGE_INTEGER time_start0;  //start time 
			_LARGE_INTEGER time_over0;   //finish time  
			double dqFreq0;      //timer frequency  
			LARGE_INTEGER f0;    //timer frequency  
			QueryPerformanceFrequency(&f0);
			dqFreq0 = (double)f0.QuadPart;
			QueryPerformanceCounter(&time_start0);

			//tempOriginalTimeSeriesPair.original_time_series_id = APCALinkOriginal[m_temp_queue_top.original_time_series_id].original_time_series_id;
			//tempOriginalTimeSeriesPair.p_original_time_series = APCALinkOriginal[m_temp_queue_top.original_time_series_id].originalLink;


			tempOriginalTimeSeriesPair.original_time_series_id = m_temp_queue_top.original_time_series_id;
			getFileStreamByID(file_name, g_time_series_length + 1, m_temp_queue_top.original_time_series_id, original_time_series);
			for (int j = 0; j < g_time_series_length; j++) {
				original_time_series[j] = original_time_series[j + 1];
			}
			tempOriginalTimeSeriesPair.d_dist = distanceEUC(g_query_time_series, g_time_series_length, original_time_series, g_time_series_length);
			//tempOriginalTimeSeriesPair.d_dist = distanceEUC(g_query_time_series, g_time_series_length, tempOriginalTimeSeriesPair.p_original_time_series, g_time_series_length);
			delete[] original_time_series;
			cout << "        temp.insert(" << tempOriginalTimeSeriesPair.d_dist << "), DLB: " << m_temp_queue_top.d_dist << endl;
			temp.push_back(tempOriginalTimeSeriesPair);

			QueryPerformanceCounter(&time_over0);    //Finish recording time  
			run_time0 = 1000000 * (time_over0.QuadPart - time_start0.QuadPart) / dqFreq0;
			g_n_time_apca_point_count += run_time0;

			//printf("\nrun_time = %f us\n", run_time0);
		}
		else if (m_temp_queue_top.p_rtree_node->IsLeaf()) {//is Leaf Node  tempQueueTop.value->IsLeaf() || tempQueueTop.swith==1
			cout << "    queue.top is Leaf Node, Leaf Node MINDIST: " << m_temp_queue_top.d_dist << ", Leaf Node level: " << m_temp_queue_top.p_rtree_node->m_level << ", Leaf Node m_account: " << m_temp_queue_top.p_rtree_node->m_count << endl;

			APCA_QUAL::APCA QProjection;
			QProjection.segmentNum = apcaRTree.NUMDIMS / 2;
			//cout << "\n" << "Q.segmentNum = " << QProjection->segmentNum << endl;
			QProjection.r = new DataType[QProjection.segmentNum];
			QProjection.v = new DataType[QProjection.segmentNum];

			double run_time1;
			_LARGE_INTEGER time_start1;  //start time 
			_LARGE_INTEGER time_over1;   //finish time  
			double dqFreq1;      //timer frequency  
			LARGE_INTEGER f1;    //timer frequency  
			QueryPerformanceFrequency(&f1);
			dqFreq1 = (double)f1.QuadPart;
			QueryPerformanceCounter(&time_start1);

			for (int i = 0; i < m_temp_queue_top.p_rtree_node->m_count; i++) {

				g_n_account_leaf_node++;

				double assignment_run_time;
				_LARGE_INTEGER assignment_time_start;  //start time 
				_LARGE_INTEGER assignment_time_over;   //finish time  
				double assignment_dqFreq;      //timer frequency  
				LARGE_INTEGER assignment_f;    //timer frequency  
				QueryPerformanceFrequency(&assignment_f);
				assignment_dqFreq = (double)assignment_f.QuadPart;
				QueryPerformanceCounter(&assignment_time_start);

				f_temp_APCA_Pair.original_time_series_id = int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data);
				//cout << tempAPCAPair.id_originalTimeSeries << endl;
				f_temp_APCA_Pair.p_APCA_point = &APCALinkOriginal[int(m_temp_queue_top.p_rtree_node->m_branch[i].m_data)].APCALink;
				//cout << "KNN : data point ID = " << tempAPCAPair.id_originalTimeSeries << endl;
				//cout << tempAPCAPair.APCAValue << endl;

				n_distanceLB_count++;

				QueryPerformanceCounter(&assignment_time_over);    //Finish recording time  
				assignment_run_time = 1000000 * (assignment_time_over.QuadPart - assignment_time_start.QuadPart) / assignment_dqFreq;
				g_n_time_leaf_node_assignment_count += assignment_run_time;

				//printf("\nrun_time = %f us\n", run_time0);

				double distanceLB_run_time;
				_LARGE_INTEGER distanceLB_time_start;  //start time 
				_LARGE_INTEGER distanceLB_time_over;   //finish time  
				double distanceLB_dqFreq;      //timer frequency  
				LARGE_INTEGER distanceLB_f;    //timer frequency  
				QueryPerformanceFrequency(&distanceLB_f);
				distanceLB_dqFreq = (double)distanceLB_f.QuadPart;
				QueryPerformanceCounter(&distanceLB_time_start);

				//f_temp_APCA_Pair.d_dist = distanceLB(QAPCAProjection(g_query_time_series, g_time_series_length, *f_temp_APCA_Pair.p_APCA_point), *f_temp_APCA_Pair.p_APCA_point);
				f_temp_APCA_Pair.d_dist = distanceLB(QAPCAProjection(g_query_time_series, g_time_series_length, *f_temp_APCA_Pair.p_APCA_point, QProjection), *f_temp_APCA_Pair.p_APCA_point);
				QueryPerformanceCounter(&distanceLB_time_over);    //Finish recording time  
				distanceLB_run_time = 1000000 * (distanceLB_time_over.QuadPart - distanceLB_time_start.QuadPart) / distanceLB_dqFreq;
				g_n_time_leaf_node_distanceLB_count += distanceLB_run_time;

				//printf("\nrun_time = %f us\n", run_time0);

				//cout << tempAPCAPair.key;
				f_temp_APCA_Pair.p_rtree_node = nullptr;
				//	cout << tempAPCAPair.key;

				double queue_run_time;
				_LARGE_INTEGER queue_time_start;  //start time 
				_LARGE_INTEGER queue_time_over;   //finish time  
				double queue_dqFreq;      //timer frequency  
				LARGE_INTEGER queue_f;    //timer frequency  
				QueryPerformanceFrequency(&queue_f);
				queue_dqFreq = (double)queue_f.QuadPart;
				QueryPerformanceCounter(&queue_time_start);

				cout << "            Push apca data. Branch id: " << i << " DLB: " << f_temp_APCA_Pair.d_dist << ", time series ID: " << f_temp_APCA_Pair.original_time_series_id << endl;
				queue.push(f_temp_APCA_Pair);
				n_push_leaf_node_count++;

				QueryPerformanceCounter(&queue_time_over);    //Finish recording time  
				queue_run_time = 1000000 * (queue_time_over.QuadPart - queue_time_start.QuadPart) / queue_dqFreq;
				g_n_time_leaf_node_push_count += queue_run_time;

				//printf("\nrun_time = %f us\n", run_time0);

				//cout << "KNN : queue.size() = " << queue.size() << endl;
				//		//cout << tempAPCAPair.key;
			}
			//	//cout << tempAPCAQueue.top().id_originalTimeSeries << ", " << tempAPCAQueue.top().key << endl;
			//	//cout << tempOriginalTimeSeriesPair.key;

			QueryPerformanceCounter(&time_over1);    //Finish recording time  
			run_time1 = 1000000 * (time_over1.QuadPart - time_start1.QuadPart) / dqFreq1;
			g_n_time_count1 += run_time1;

			//printf("\nrun_time = %f us\n", run_time0);
			delete[] QProjection.r;
			QProjection.r = nullptr;
			delete[] QProjection.v;
			QProjection.v = nullptr;
		}
		else {															//is internal Node
			cout << "    queue.top is Internal Node, MINDIST: " << m_temp_queue_top.d_dist << ", Internal Node level: " << m_temp_queue_top.p_rtree_node->m_level << ", Internal Node m_account: " << m_temp_queue_top.p_rtree_node->m_count << endl;
			double run_time2;
			_LARGE_INTEGER time_start2;  //start time 
			_LARGE_INTEGER time_over2;   //finish time  
			double dqFreq2;      //timer frequency  
			LARGE_INTEGER f2;    //timer frequency  
			QueryPerformanceFrequency(&f2);
			dqFreq2 = (double)f2.QuadPart;
			QueryPerformanceCounter(&time_start2);

			APCA_NODE_PAIR tempApcaPair;

			for (int branch_index = 0; branch_index < m_temp_queue_top.p_rtree_node->m_count; branch_index++) {

				g_n_account_child_node++;

				regionG fs_region_G;
				fs_region_G.regionNum = apcaRTree.NUMDIMS / 2;
				//cout << "regionNum = " << G.regionNum << endl;
				fs_region_G.G1 = new double[apcaRTree.NUMDIMS / 2];
				fs_region_G.G2 = new double[apcaRTree.NUMDIMS / 2];
				fs_region_G.G3 = new double[apcaRTree.NUMDIMS / 2];
				fs_region_G.G4 = new double[apcaRTree.NUMDIMS / 2];

				double mindistQR_run_time;
				_LARGE_INTEGER mindistQR_time_start;  //start time 
				_LARGE_INTEGER mindistQR_time_over;   //finish time  
				double mindistQR_dqFreq;      //timer frequency  
				LARGE_INTEGER mindistQR_f;    //timer frequency  
				QueryPerformanceFrequency(&mindistQR_f);
				mindistQR_dqFreq = (double)mindistQR_f.QuadPart;
				QueryPerformanceCounter(&mindistQR_time_start);


				tempApcaPair.d_dist = MINDISTQR(g_query_time_series, g_time_series_length, getRegionG(m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_rect, fs_region_G));

				/*for (int i = 0; i < fs_region_G.regionNum; i++) {
				cout << "            G1[" << i << "] = " << fs_region_G.G1[i] << ", G2[" << i << "] = " << fs_region_G.G2[i] << ", G3[" << i << "] = " << fs_region_G.G3[i] << ", G4[" << i << "] = " << fs_region_G.G4[i] << endl;
				}*/

				QueryPerformanceCounter(&mindistQR_time_over);    //Finish recording time  
				mindistQR_run_time = 1000000 * (mindistQR_time_over.QuadPart - mindistQR_time_start.QuadPart) / mindistQR_dqFreq;
				g_n_time_child_node_MINDIST_count += mindistQR_run_time;

				//printf("\nrun_time = %f us\n", run_time0);

				tempApcaPair.p_rtree_node = m_temp_queue_top.p_rtree_node->m_branch[branch_index].m_child;

				double push_run_time;
				_LARGE_INTEGER push_time_start;  //start time 
				_LARGE_INTEGER push_time_over;   //finish time  
				double push_dqFreq;      //timer frequency  
				LARGE_INTEGER push_f;    //timer frequency  
				QueryPerformanceFrequency(&push_f);
				push_dqFreq = (double)push_f.QuadPart;
				QueryPerformanceCounter(&push_time_start);

				cout << "            push internal node, Branch id: " << branch_index << " internal node MINDIST: " << tempApcaPair.d_dist << ", internal node level: " << tempApcaPair.p_rtree_node->m_level << ", internal node m_count: " << tempApcaPair.p_rtree_node->m_count << endl;
				queue.push(tempApcaPair);
				//cout << "KNN : queue.size() = " << queue.size() << endl;

				QueryPerformanceCounter(&push_time_over);    //Finish recording time  
				push_run_time = 1000000 * (push_time_over.QuadPart - push_time_start.QuadPart) / push_dqFreq;
				g_n_time_child_node_push_count += push_run_time;

				delete[] fs_region_G.G1;
				fs_region_G.G1 = nullptr;
				delete[] fs_region_G.G2;
				fs_region_G.G2 = nullptr;
				delete[] fs_region_G.G3;
				fs_region_G.G3 = nullptr;
				delete[] fs_region_G.G4;
				fs_region_G.G4 = nullptr;

				//printf("\nrun_time = %f us\n", run_time0);
			}

			QueryPerformanceCounter(&time_over2);    //Finish recording time  
			run_time2 = 1000000 * (time_over2.QuadPart - time_start2.QuadPart) / dqFreq2;
			g_n_time_count2 += run_time2;

			//printf("\nrun_time = %f us\n", run_time0);
		}
	}

	return false;
}




#endif
