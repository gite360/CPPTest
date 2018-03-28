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
#include <cstdint>
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
	DataType* r = nullptr;
	DataType* v = nullptr;
};

struct cmpMoreDistance {
	bool operator () (const APCA_MERGE &a, const APCA_MERGE &b) {
		return a.distance_AE > b.distance_AE;             // decrease
	}
};

struct cmpLessDistance {
	bool operator () (const APCA_MERGE &a, const APCA_MERGE &b) {
		return a.distance_AE < b.distance_AE;             // Increase
	}
};

struct cmpMoreIndex {
	bool operator () (const APCA_MERGE &a, const APCA_MERGE &b) {
		return *a.r > *b.r;             // decrease
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


double getPAA(const DataType* x, const DataType& n, const DataType& N, DataType *y) {

	if (n < N) {                                          //n must longer than N
		cout << "N has to smaller than n!!!!!!!";
		return false;
	}

	double Remainder = int(n) % int(N);
	//cout << "Remainder = " << Remainder << ", n = " << n << endl;
	double integerDividend = n - Remainder;
	double integerQuotient = integerDividend / N;
	double integerSegmentLength = integerQuotient + 1;

	int j = 0;
	int endOfLongSegment = 0, indexOfLongSegment = 0;
	int OnOff = 0;
	for (int i = 1; i <= N; i++) {
		//cout << endl;
		double sum = 0;
		if (i <= Remainder) {
			for (j = (i - 1) * integerSegmentLength; j < i * integerSegmentLength; j++) {
				sum += x[j];
				//cout << "x[" << j << "] = " << x[j] << " ";
			}
			y[i - 1] = sum / integerSegmentLength;
			//cout << endl << "y[" << i - 1 << "] = " << y[i - 1] << endl;
			indexOfLongSegment = i;
			endOfLongSegment = i * integerSegmentLength;
		}
		else {
			for (j = endOfLongSegment + (i - indexOfLongSegment - 1)*integerQuotient; j < endOfLongSegment + (i - indexOfLongSegment)*integerQuotient; j++) {
				sum += x[j];
				cout << "x[" << j << "] = " << x[j] << " ";
			}
			y[i - 1] = sum / integerQuotient;
			//cout << endl << "y[" << i - 1 << "] = " << y[i - 1] << endl;
		}
	}

}


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
	cout << "original: " << *originalArray << ", " << arrayLength << endl;
	return accumulate(originalArray, originalArray + int(arrayLength), 0.0) / T(arrayLength);
}

double getSegmentMeanValue(const int& merge_start_index, const APCA& apca_presentation) {
	if (merge_start_index != 0)
		return  (apca_presentation.v[merge_start_index] * (apca_presentation.r[merge_start_index] - apca_presentation.r[merge_start_index - 1]) + apca_presentation.v[merge_start_index + 1] * (apca_presentation.r[merge_start_index + 1] - apca_presentation.r[merge_start_index])) / (apca_presentation.r[merge_start_index + 1] - apca_presentation.r[merge_start_index - 1]);
	else
		return  (apca_presentation.v[merge_start_index] * (apca_presentation.r[merge_start_index]+1.0) + apca_presentation.v[merge_start_index + 1] * (apca_presentation.r[merge_start_index + 1] - apca_presentation.r[merge_start_index])) / (apca_presentation.r[merge_start_index + 1] +1.0);
}

bool getExactMeanValue(const DataType* orginal_time_series, APCA& apca_presentation) {

	double fd_segment_length = NULL;
	int i = 0;

	cout << "Exact mean value: " << endl;

	if (apca_presentation.segmentNum == 1) {
		return true;
	}

	apca_presentation.v[i] = getAve(orginal_time_series, apca_presentation.r[0] + 1);
	cout << apca_presentation.r[i] << " " << apca_presentation.v[i] << endl;

	for (i = 1; i < apca_presentation.segmentNum; i++) {
		apca_presentation.v[i] = getAve(orginal_time_series + int(apca_presentation.r[i - 1] + 1), double(apca_presentation.r[i] - apca_presentation.r[i - 1]));
		cout << apca_presentation.r[i] << " " << apca_presentation.v[i] << endl;;
	}
	return true;
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

bool combinateSegments0(const DataType* orginal_time_series, const APCA& merge_apca_presentation, int segment_number, const int& n, const int& M, int merge_number, int merge_index, const int& last_merge_point, APCA_MERGE& minimum_distanceAE_merge) {

	if (merge_number == 0) {
		double temp_distance_AE = distanceAE(orginal_time_series, merge_apca_presentation);

		if (minimum_distanceAE_merge.distance_AE > temp_distance_AE) {
			minimum_distanceAE_merge.distance_AE = temp_distance_AE;
			//memcpy(minimum_distanceAE_merge.r, merge_apca_presentation.r, sizeof(double)*merge_apca_presentation.segmentNum);
			//memcpy(minimum_distanceAE_merge.v, merge_apca_presentation.v, sizeof(double)*merge_apca_presentation.segmentNum);

			copy(merge_apca_presentation.r , merge_apca_presentation.r + merge_apca_presentation.segmentNum, minimum_distanceAE_merge.r);
			copy(merge_apca_presentation.v , merge_apca_presentation.v + merge_apca_presentation.segmentNum, minimum_distanceAE_merge.v);
		}

		/*for (int i = 0; i < M; i++) {
			cout << merge_apca_presentation.r[i] << " ";
			result_vector.push_back(merge_apca_presentation.r[i]);
			result_vector.push_back(merge_apca_presentation.v[i]);
		}
		result_vector.push_back(distanceAE(orginal_time_series, merge_apca_presentation));*/

		return true;
	}
	else {

		APCA temp_apca;
		temp_apca.r = new DataType[segment_number];
		temp_apca.v = new DataType[segment_number];
		temp_apca.segmentNum = merge_apca_presentation.segmentNum - 1;

		while (merge_index < merge_apca_presentation.segmentNum - 1) {
			//cout << merge_apca_presentation.r[merge_index] << " " << last_merge_point << endl;
			cout << "merge_index: " << merge_index << endl;
			if (merge_apca_presentation.r[merge_index] >= last_merge_point) {
				//cout << ">>>>>>>" << endl;
				/*for (int i = 0; i < merge_index; i++) {
					temp_apca.r[i] = merge_apca_presentation.r[i];
					temp_apca.v[i] = merge_apca_presentation.v[i];
				}*/

				copy(merge_apca_presentation.r, merge_apca_presentation.r + merge_index, temp_apca.r);
				copy(merge_apca_presentation.v, merge_apca_presentation.v + merge_index, temp_apca.v);

				//temp_apca.v[merge_index] = getAve(orginal_time_series + int(temp_apca.r[merge_index - 1] + 1), DataType(merge_apca_presentation.r[merge_index + 1] - merge_apca_presentation.r[merge_index - 1]));
				temp_apca.v[merge_index] = getSegmentMeanValue(merge_index, merge_apca_presentation);
				temp_apca.r[merge_index] = merge_apca_presentation.r[merge_index + 1];

				int i = merge_index + 1;

				/*while (i < temp_apca.segmentNum) {
					temp_apca.r[i] = merge_apca_presentation.r[i + 1];
					temp_apca.v[i] = merge_apca_presentation.v[i + 1];
					i++;
				}*/
				
				copy(merge_apca_presentation.r + merge_index + 2, merge_apca_presentation.r + merge_apca_presentation.segmentNum, temp_apca.r + merge_index + 1);
				copy(merge_apca_presentation.v + merge_index + 2, merge_apca_presentation.v + merge_apca_presentation.segmentNum, temp_apca.v + merge_index + 1);
				
				int last_merge_point = merge_apca_presentation.r[merge_index];
				combinateSegments0(orginal_time_series, temp_apca, segment_number - 1, n, M, merge_number - 1, 0, last_merge_point, minimum_distanceAE_merge);

			}
			merge_index++;
		}
		delete[] temp_apca.r;
		delete[] temp_apca.v;
	}
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


uint64_t Combinations(unsigned int n, unsigned int r)
{
	if (r > n)
		return 0;

	/** We can use Pascal's triange to determine the amount
	* of combinations. To calculate a single line:
	*
	* v(r) = (n - r) / r
	*
	* Since the triangle is symmetrical, we only need to calculate
	* until r -column.
	*/

	uint64_t v = n--;

	for (unsigned int i = 2; i < r + 1; ++i, --n)
		v = v * n / i;

	return v;
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
	cout << "mergeSegmentsRecursivly()" << endl;
	if (apca_presentation.segmentNum > M) {
		int merge_frequency = apca_presentation.segmentNum - M;
		cout << "merge_frequency: " << merge_frequency << endl;

		APCA_MERGE minimum_distanceAE_merge;
		minimum_distanceAE_merge.r = new DataType[M];
		minimum_distanceAE_merge.v = new DataType[M];
		minimum_distanceAE_merge.distance_AE = DBL_MAX;

		combinateSegments0(orginal_time_series, apca_presentation, apca_presentation.segmentNum, n, M, merge_frequency, 0, 0, minimum_distanceAE_merge);
		
		for (int i = 0; i < M; i++) {
			cout << minimum_distanceAE_merge.r[i] << ", " << minimum_distanceAE_merge.v[i] << endl;
		}
		cout << endl;

		delete[] minimum_distanceAE_merge.r;
		delete[] minimum_distanceAE_merge.v;

	}
	else if (apca_presentation.segmentNum == M) {
		cout << "M = apca_presentation.segmentNum" << endl;
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

void getMergeIndexQueue(const int& merge_frequency, APCA& apca_presentation) {
	cout << "getMergeIndexQueue()" << endl;
	priority_queue<APCA_MERGE, vector<APCA_MERGE>, cmpLessDistance> queue_difference;
	priority_queue<APCA_MERGE, vector<APCA_MERGE>, cmpMoreIndex> queue_index;
	APCA_MERGE temp_merge_index;
	//temp_merge_index.r = new DataType[merge_frequency];
	//temp_merge_index.r = new DataType;
	temp_merge_index.r = new DataType[apca_presentation.segmentNum + 1];
	DataType *pointer_beginer_r = temp_merge_index.r;
	//DataType *pointer_beginer_v = temp_merge_index.v;
	int merge_start_index = NULL;

	cout << "segment Number: " << apca_presentation.segmentNum << endl;
	merge_start_index = 0;
	for (int j = 0; j < apca_presentation.segmentNum - 1; j++) {
		cout << "j: " << j << endl;
		//cout << "sefments_difference: " << sefments_difference << " sgetments diff: " << fabs(apca_presentation.v[j + 1] - apca_presentation.v[j]) << endl;
		temp_merge_index.distance_AE = fabs(apca_presentation.v[j + 1] - apca_presentation.v[j]);
		cout << "distance_AE: " << temp_merge_index.distance_AE << endl;
		if (queue_difference.size() < merge_frequency) {
			//temp_merge_index.r[merge_start_index] = apca_presentation.r[j];
			*temp_merge_index.r = apca_presentation.r[j];
			cout << "r: " << *temp_merge_index.r << endl;
			queue_difference.push(temp_merge_index);
			//merge_start_index++;
			temp_merge_index.r++;
			//queue_difference.emplace(temp_merge_index);
		}
		else if (queue_difference.top().distance_AE > temp_merge_index.distance_AE) {
			queue_difference.pop();
			//temp_merge_index.r[merge_start_index] = apca_presentation.r[j];
			*temp_merge_index.r = apca_presentation.r[j];;
			cout << "********r: " << *temp_merge_index.r << endl;
			//queue_difference.emplace(temp_merge_index);
			queue_difference.push(temp_merge_index);
			//merge_start_index++;
			temp_merge_index.r++;
		}
	}

	if (queue_difference.size() != merge_frequency) {
		cout << "getMergeIndexQueue()： Wrong！！！！！！！！！！！！！！！" << endl;
	}

	while (!queue_difference.empty()) {
		cout << "queue.r: " << *queue_difference.top().r << endl;
		queue_index.emplace(queue_difference.top());
		queue_difference.pop();
	}

	for (int i = 0; i < apca_presentation.segmentNum - 1 && (!queue_index.empty()); i++) {
		cout << "apca_presentation.r[i]: " << apca_presentation.r[i] << ", queue.top.r: " << *queue_index.top().r << endl;
		if (apca_presentation.r[i] == *queue_index.top().r) {

			cout << "i:" << i << ", apca_presentation.r[i]: " << apca_presentation.r[i] << endl;
			if (i != 0) {
				apca_presentation.v[i] = (apca_presentation.v[i] * (apca_presentation.r[i] - apca_presentation.r[i - 1]) + apca_presentation.v[i + 1] * (apca_presentation.r[i + 1] - apca_presentation.r[i])) / (apca_presentation.r[i + 1] - apca_presentation.r[i - 1]);
				cout << apca_presentation.v[i] << endl;

			}
			else {
				apca_presentation.v[i] = (apca_presentation.v[i] * double(apca_presentation.r[i] + 1.0) + apca_presentation.v[i + 1] * (apca_presentation.r[i + 1] - apca_presentation.r[i])) / double(apca_presentation.r[i + 1] + 1);
				cout << "*********" << apca_presentation.v[i] << endl;
			}
			apca_presentation.r[i] = apca_presentation.r[i + 1];
			merge_start_index = i + 1;
			while (merge_start_index < apca_presentation.segmentNum) {
				apca_presentation.r[merge_start_index] = apca_presentation.r[merge_start_index + 1];
				apca_presentation.v[merge_start_index] = apca_presentation.v[merge_start_index + 1];
				merge_start_index++;
			}

			queue_index.pop();
			i--;
		}
	}

	apca_presentation.segmentNum -= merge_frequency;
	temp_merge_index.r = pointer_beginer_r;
	//temp_merge_index.v = pointer_beginer_v;
	delete[] temp_merge_index.r;
	//delete[] temp_merge_index.v;
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
			cout << apca_presentation.v[merge_start_index] << endl;

			apca_presentation.v[merge_start_index] = (apca_presentation.v[merge_start_index] * ((apca_presentation.r[merge_start_index] - apca_presentation.r[merge_start_index - 1]) / initial_segment_length) + apca_presentation.v[merge_start_index + 1] * ((apca_presentation.r[merge_start_index + 1] - apca_presentation.r[merge_start_index]) / initial_segment_length)) / ((apca_presentation.r[merge_start_index + 1] - apca_presentation.r[merge_start_index - 1]) / initial_segment_length);
			//cout << temp1 << " *************** " << apca_presentation.v[merge_start_index] << endl;
			//if (temp != apca_presentation.v[merge_start_index]) cout << "ERROROROROROOROROROROROO! " << endl;
		}
		else {
			//.v[merge_start_index] = getAve(orginal_time_series, DataType(apca_presentation.r[merge_start_index + 1]));
			//double temp = apca_presentation.v[merge_start_index];
			apca_presentation.v[merge_start_index] = (apca_presentation.v[merge_start_index] * ((apca_presentation.r[merge_start_index] + 1) / initial_segment_length) + apca_presentation.v[merge_start_index + 1] * ((apca_presentation.r[merge_start_index + 1] - apca_presentation.r[merge_start_index]) / initial_segment_length)) / (apca_presentation.r[merge_start_index + 1]);
			//if (temp != apca_presentation.v[merge_start_index]) cout << "ERROROROROROOROROROROROO!!!!! " << endl;
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
	cout << "mergeSegmentsIterator()" << endl;
	if (apca_presentation.segmentNum > M) {
		int merge_frequency = apca_presentation.segmentNum - M;

		//getMergeIndex(orginal_time_series, merge_frequency, apca_presentation);

		getMergeIndexQueue(merge_frequency, apca_presentation);

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
		cout << "apca_presentation.segmentNum == M" << endl;
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
	cout << "Total Resolution: " << log2(getNextPowerOf2(n)) << ",  Resolution: " << log2(retained_coeffs_length) << endl;

	truncateZero(n, length_power_of_2, general_apca_presentation);
	getExactMeanValue(orginal_time_series, general_apca_presentation);
	//mergeSegmentsIterator(orginal_time_series, general_apca_presentation, n, N, italicC);
	mergeSegmentsRecursivly(orginal_time_series, general_apca_presentation, n, N, italicC);

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
#define ARRAY_LENGTH 9
	DataType n = ARRAY_LENGTH, M = 3;
	DataType test_array[ARRAY_LENGTH] = { 7,5,5,3,3,3,4,6,6 };
	//DataType test_array[ARRAY_LENGTH] = { -1,0.50206,0.54216,0.72238,1.4289,2.1365,2.2811,1.9363,1.4689,1.0088,0.38028,-0.29678,-0.51393,-0.25564,-0.1072,-0.28783,-0.41801,-0.31916,-0.26038,-0.35036,-0.50549,-0.71089,-0.82392,-0.8997,-1.1539,-1.2298,-1.0441,-1.202,-1.3922,-1.1301,-1.1799,-1.6493,-1.7266,-1.6084,-1.6628,-1.6507,-1.6973,-1.8387,-1.8026,-1.7805,-1.8252,-1.6448,-1.4238,-1.3922,-1.3604,-1.2002,-0.91863,-0.68592,-0.66794,-0.51272,-0.10169,0.063954,0.082614,0.23761,0.17479,0.12321,0.50339,0.68387,0.47499,0.5328,0.72355,0.66442,0.64794,0.75705,0.73207,0.62021,0.6063,0.67795,0.68908,0.59881,0.54265,0.58181,0.63063,0.66442,0.68632,0.65738,0.65089,0.72174,0.73971,0.69148,0.69523,0.75011,0.85384,0.91602,0.82775,0.88091,1.0888,0.93104,0.6103,0.63889,0.68468,0.58324,0.64052,0.70859,0.70501,0.71382,0.43376 };
	DataType* PAA_time_series = new DataType[M];
	//DataType test_array[ARRAY_LENGTH] = { 0.80467, 0.36777, 0.24389, 0.026614, -0.2744, 0.096731, -0.74773, -1.6098 };
	//DataType* mp_HDWT_time_series = new DataType[8];

	APCA test_apca;
	test_apca.segmentNum = M;
	test_apca.r = new DataType[M];
	test_apca.v = new DataType[M];
	unsigned int after_padding_length = NULL;
	int retained_coeffs_length = NULL;

	//nChoosek(97,72);
	//cout << Combinations(96, 72) << endl;

	//getHDWT(test_array, getNextPowerOf2(8), M, mp_HDWT_time_series, retained_coeffs_length);
	//getPAA(test_array, n, M, PAA_time_series);
	getAPCAPoint(test_array, n, M, test_apca);
	cout << endl;
	/*for (int i = 0; i < M; i++) {
		cout << PAA_time_series[i] << " ";
	}
	cout << endl;*/
	/*for (int i = 0; i < M; i++) {
		cout << test_apca.r[i] << " " << test_apca.v[i] << endl;
	}*/
	cout << "M: " << M << endl;
	//distanceAE(test_array, test_apca);

	//getCombinationNumber(7, 2);
	//test_apca.r = nullptr;
	delete[] PAA_time_series;
	delete[] test_apca.r;
	//test_apca.v = nullptr;
	delete[] test_apca.v;
	system("pause");
	return 0;
}

