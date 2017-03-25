#pragma once
#include "Array2D.h"
#include <omp.h>

#define POSLABEL 1
#define NEGLABEL -1
#define LESSEQ 1
#define GREATER -1

template<typename T>
class DecisionStump {
private:
	//member
	Array2D<T> m_aFeature;
	vector<char> m_vLabel;
	Array2D<unsigned int> m_aSortedLabel;
	vector<double> m_vWeight;
	char m_cParity;   -1<=   +1>
	T m_tSplit;
	double m_dTrainErr;
	double m_dTestErr;
	//method

public:
	//member
	//method
	DecisionStump() = default;
	DecisionStump();
	Save();
	Train();
	Test();



};

#define NEW 0
#if NEW
template <typename T=double, typename U=char>
class DecisionStump {
private:
	//input: before train
	T * m_feature = nullptr;
	U * m_label = nullptr;
	unsigned int *m_sortMap = nullptr;
	double *m_weight = nullptr;
	unsigned int m_sampleCounts = 0u;
	unsigned int m_featureCounts = 0u;
	//output: after train
	unsigned int m_featureIndex = 0u;
	char m_parity = LESSEQ;
	T m_split = 0;
	double m_trainErr = 1.0;
	double m_initErr = 0.0;

public:
	DecisionStump() {};
	DecisionStump(const T *fptr, const U *lptr, const unsigned int *smapptr, const double *wptr, const unsigned int samCounts, const unsigned int feaCounts)
		:m_feature(fptr), m_label(lptr), m_sortMap(smapptr), m_weight(wptr), m_sampleCounts(samCounts), m_featurCounts(feaCounts);
	void Train();
	void Test();
};


template<typename T, typename U>
inline DecisionStump<T, U>::DecisionStump(const T * fptr, const U * lptr, const unsigned int * smapptr, const double * wptr, const unsigned int samCounts, const unsigned int feaCounts){
	m_initErr = 0.0;
	for (unsigned int i = 0; i != m_sampleCounts; ++i)
		m_initErr += m_label[i] == POSLABEL ? m_weight[i] : 0;
	m_initErr < 1.0 - m_initErr ? (m_parity = LESSEQ, m_trainErr = m_initErr) : (m_parity = GREATER, m_trainErr = 1 - m_initErr);
}

template<typename T, typename U>
inline void DecisionStump<T, U>::Train(){
	vector<T> split(m_featureCounts, 0);
	vector<double> err(m_featureCounts, 1);
#pragma omp parallel for
	for (unsigned int x = 0; x != m_featureCounts; ++x) {
		unsigned int *sortMapStart = m_sortMap + x*m_sampleCounts;
		T * featureStart = m_feature + x*m_sampleCounts;
		split[x] = featureStart[sortMapStart[0]];
		double err = m_initErr;
		char parity = m_parity;
		for (unsigned int y = 0; y != m_sampleCounts; ++y) {
			data_index = sortMapStart[y];
			label = m_label[sortMapStart[y]];
			data = m_feature[featureStart[sortMapStart[y]]];
			weight = m_weight[sortMapStart[y]];
			res data <= split ? -1 : 1;
			if parity==LESSEQ
				data
			err
			if label==POSLABEL
				err-weight


			split[x] = featureStart[sortMapStart[y]];
			featureStart[sortMapStart[y]] <= split[x] ? ;
			label_real= m_label[]
			featureStart[sortMapStart[y]]
			
		}
		x*m_sampleCounts+y
	}


}

template<typename T, typename U>
inline void DecisionStump<T, U>::Test(){


}

#endif