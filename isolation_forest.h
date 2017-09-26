/*M///////////////////////////////////////////////////////////////////////////////////////
//
//  IMPORTANT: READ BEFORE DOWNLOADING, COPYING, INSTALLING OR USING.
//
//  By downloading, copying, installing or using the software you agree to this license.
//  If you do not agree to this license, do not download, install,
//  copy or use the software.
//
//
//                           License Agreement
//                For Open Source Computer Vision Library
//
// Copyright 864051262@qq.com, all rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
//   * Redistribution's of source code must retain the above copyright notice,
//     this list of conditions and the following disclaimer.
//
//   * Redistribution's in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//
//   * The name of the copyright holders may not be used to endorse or promote products
//     derived from this software without specific prior written permission.
//
// This software is provided by the copyright holders and contributors "as is" and
// any express or implied warranties, including, but not limited to, the implied
// warranties of merchantability and fitness for a particular purpose are disclaimed.
// In no event shall the Intel Corporation or contributors be liable for any direct,
// indirect, incidental, special, exemplary, or consequential damages
// (including, but not limited to, procurement of substitute goods or services;
// loss of use, data, or profits; or business interruption) however caused
// and on any theory of liability, whether in contract, strict liability,
// or tort (including negligence or otherwise) arising in any way out of
// the use of this software, even if advised of the possibility of such damage.
//
//M*/

//this code was complete a year ago ,then I was learning single class classification,I could found one_class svm,svdd's code ,butcounld not found 
//isolation forest,so I write it my self 
//this code is complete without practice test but just simple sin function test
//if you have any question,please feel free to contact me :864051262@qq.com


#pragma once
#include <map>
#include <iostream>
#include <stdio.h>
#include <random>
using namespace std;

#define DELETE_PTR(p) { if(NULL != p) { delete p; p = NULL;} }
#define DELETE_ARR_PTR(p) { if(NULL != p) { delete []p; p = NULL;} }
typedef struct min_max{
	float* p_min;
	float* p_max;
}min_max;

typedef struct sample_str{
	int sample_num;//total number of samples
	int feature_length;//length of feature
	//int feature_index;
	int train_num;//training sample number
	int max_tree_depth;//max depth
	int node_id;//node id
	sample_str(const int& num,const int& length,const int& train_nums,const int& id)
	{
		sample_num = num;
		feature_length = length;
		train_num = train_nums;
		p_num = new int[train_nums]();
		p_index = new int[train_nums]();
		p_sub_index0 = new int[train_nums]();
		p_sub_index1 = new int[train_nums]();
		//p_min = new float[feature_length]();
		//p_max = new float[feature_length]();
		memset(subnode,0,sizeof(subnode));
		node_id = id;
	}

	void get_random_sample()
	{   
		//get random number,validate under vs2015,other version of vs is unknosn
		std::random_device rd0 ; 
		std::uniform_real_distribution<> u0(1,sample_num);
		for (int i=0;i<train_num;i++)
		{
			bool isthesame=false;
			do 
			{
				p_index[i] = u0(rd0);
				//the comment code makes the sampling become bootstrap
				//i.e,sampled with replacement
				/*for (int j=0;j<i;j++)
				{
					isthesame = (p_index[i]==p_index[j])?true:false;
					if (isthesame)
					{
						break;
					}
				}*/
			} while (isthesame);
		}
		
	}

	~sample_str()
	{
		DELETE_ARR_PTR(p_num);
		DELETE_ARR_PTR(p_index);
		DELETE_ARR_PTR(p_sub_index0);
		DELETE_ARR_PTR(p_sub_index1);
		//DELETE_ARR_PTR(p_min);
		//DELETE_ARR_PTR(p_max);
	}
	int* p_num;//save the left and right sub node number
	int* p_index;//index to the training sample
	//T threshold;
	int subnode[2];//sub node sample number
	int* p_sub_index0;//index to the left node sample
	int* p_sub_index1;
	//float* p_min;
	//float* p_max;
	min_max value_minmax; //the min and max of the current node
}sample_str;


typedef struct ret_result
{
	bool is_leaf;//is node or not
	float path_length;//path from the root to the leaf
	int child_id;//child node id,if there is any
}ret_result;
///a node in a tree
template<typename T>
class CBinaryNode {
public:
    CBinaryNode();
	bool& isLeaf() const{return m_isLeaf};
	//read model
    void deserialize(std::istream&);

	T threshold() const{return m_threshold};
    
    //feature:sample features,sample_atr:sample information
	int learn( T** feature, sample_str& sample_atr);   
	int predict(T* feature,ret_result& result);
	//save node information
    void serialize(std::ostream&);
    //feature index used by current node
    int m_featureIndex;
    //node depth in a tree 
	int m_depth_count;
	//total sample number 
	int m_sample_num; 
	//node id,m_node_id[0]stands for current node is ，1 represent left node id
	//2 represent right child node
	int m_node_id[3];
	//thresholds used for split
    T m_threshold;
    //is leaf or not
	bool m_isLeaf;
private:
	//compare the feature with thresholds,return 0 or 1,
	//0 represent the sample goto left node,1 represent
	//the sample goto right node
	short comparebyFeature(const T&);
};

//predict
template<typename T>
int CBinaryNode<T>::predict(T* feature,ret_result& result)
{
	//result.path_length++;
	if (m_isLeaf){
		result.is_leaf = true;
		result.child_id = -1;
		//result.size = 
	} 
	else{
		result.path_length++;
		result.child_id = comparebyFeature(feature[m_featureIndex]);
	}

	//this->m_node_id[result.child_id+1]
	return 0;
}
/// Serialization.
template<typename T>
void CBinaryNode<T>::serialize(std::ostream& s)
{
    s << " " << m_node_id[0];
    s << " " << m_node_id[1];
    s << " " << m_node_id[2];
    s << " " << m_featureIndex;
    s << " " << m_threshold;
    s << " " << m_isLeaf;
    s << " " << m_depth_count;
    s << " " << m_sample_num;
}

/// Deserialization.
template<typename T>
void CBinaryNode<T>::deserialize(std::istream& s)
{
    s >> m_node_id[0];
    s >> m_node_id[1];
	s >> m_node_id[2];
    s >> m_featureIndex;
    s >> m_threshold;
	s >> m_isLeaf;
	s >> m_depth_count;
	s >> m_sample_num;
}

template <typename T> 
int CBinaryNode<T>::learn(T ** feature, sample_str& sample_atr)
{
	
	std::random_device rd ; 
	//std::default_random_engine e(rd) ; 
	std::uniform_int_distribution<> u(0,sample_atr.feature_length-1);

	m_featureIndex = u(rd);//u(e);
	m_sample_num = sample_atr.train_num;
	//printf("m_featureIndex is %d\n",m_featureIndex);
	
	//use local sample's min max not global min max of all sample
	float min = feature[sample_atr.p_index[0]][m_featureIndex];
	float max = feature[sample_atr.p_index[0]][m_featureIndex];
	for(int i=0;i<m_sample_num;i++) {
		if (feature[sample_atr.p_index[i]][m_featureIndex] > max)
			max = feature[sample_atr.p_index[i]][m_featureIndex];

		if (feature[sample_atr.p_index[i]][m_featureIndex] < min)
			min = feature[sample_atr.p_index[i]][m_featureIndex];
	}

	std::random_device rd0 ; 
	//std::default_random_engine e0(rd0) ;
	//std::uniform_real_distribution<> u0(sample_atr.value_minmax.p_min[m_featureIndex],sample_atr.value_minmax.p_max[m_featureIndex]);
	std::uniform_real_distribution<> u0(min,max);

	m_threshold = u0(rd0);//u0(e0);
	

	/*if (abs(m_threshold)>2)
	{
		printf("m_threshold is %f\n",m_threshold);
		printf("sample_atr.p_min[m_featureIndex] is %f\n",sample_atr.p_min[m_featureIndex]);
		printf("sample_atr.p_max[m_featureIndex] is %f\n",sample_atr.p_max[m_featureIndex]);
	}*/

	memset(sample_atr.p_num,0,sizeof(int)*m_sample_num);
	int count = 0;
	for(int i=0;i<m_sample_num;i++)
	{
		sample_atr.p_num[i] = comparebyFeature(feature[sample_atr.p_index[i]][m_featureIndex]);
		count += sample_atr.p_num[i];
		if (!sample_atr.p_num[i]){
			sample_atr.p_sub_index0[sample_atr.subnode[0]++] = sample_atr.p_index[i];
		}
		else{
			sample_atr.p_sub_index1[sample_atr.subnode[1]++] = sample_atr.p_index[i];
		}
	}
	m_isLeaf = (count==0 || count==m_sample_num)?true:false;

	return 0;

}

template <typename T> 
short CBinaryNode<T>::comparebyFeature(const T& feature)
{
	return feature<m_threshold?0:1;
}

template<typename T> 
CBinaryNode<T>::CBinaryNode():m_isLeaf(false)
{
	m_node_id[0] = -1;
	m_node_id[1] = -1;
	m_node_id[2] = -1;
	m_depth_count =-1;
}

/// A binary tree.
template<typename T>
class CBinaryTree {
public:
    typedef CBinaryNode<T> BinaryNodeType;

    CBinaryTree();
	~CBinaryTree();
	void deserialize(std::istream&);  
	 int learn( T ** feature, sample_str& sample_atr);   
    void serialize(std::ostream&);
	int predict(T* feature,ret_result& result);
	int m_node_count;
	inline int get_depth(){ return m_depth_count;};
private:
	int comparebyFeature(const T&,short& subnode);
	int recursiveLearn(T ** feature, sample_str& sample_atr);

	int m_sample_num; 
	int m_count_leaf;
	bool m_isend;
	float ajustment(float pathlength);
	int m_train_num;
	int m_depth_count;
	std::map<int,BinaryNodeType> m_binarynode;
};

template <typename T> 
CBinaryTree<T>::CBinaryTree():m_node_count(0),m_depth_count(0),m_count_leaf(0)
{	
}

template <typename T>
CBinaryTree<T>::~CBinaryTree()
{

}

template <typename T>
float CBinaryTree<T>::ajustment(float pathlength)
{	
	if (pathlength<2) return 0;
	if (pathlength==2) return 1;
	if (pathlength>2) {		
		return 2.0*(log(1.0*pathlength-1.0)+0.5772156649 ) - 2.0*(1.0*pathlength-1.0)/(pathlength*1.0);		
	}	
}
template <typename T>
int CBinaryTree<T>::predict(T* feature,ret_result& result)
{
	
	result.child_id=0;
	result.is_leaf = false;
	result.path_length = 0;
	float size;
	if (m_node_count==0)
	{
		return 0;
	}
	m_binarynode[0].predict(feature,result);
	if (result.is_leaf)
	{
		//result.path_length = result.path_length+m_depth_count;
		//result.path_length--;
		result.path_length = ajustment(m_binarynode[0].m_sample_num);
		//result.path_length = ajustment(m_binarynode[0].m_depth_count);
		//result.path_length = ajustment(result.path_length);
		//result.path_length = 1;
		return 0;
	}
	int id = m_binarynode[0].m_node_id[result.child_id+1];
	for (int i=0;i<10;i++)
	{
		//id = m_binarynode[i].m_node_id[result.child_id+1];
		//if (id<0)
		//{
		//	break;
		//}
		m_binarynode[id].predict(feature,result);

		size = m_binarynode[id].m_sample_num;
		//size = m_binarynode[id].m_depth_count;
		if (result.is_leaf)
		{
			break;
		}
		id = m_binarynode[id].m_node_id[result.child_id+1];
		if (id==-1)
		{
			break;
		}
	}
	
	result.path_length = ajustment(size);
	//result.path_length = ajustment(result.path_length);
	//result.path_length = ajustment(result.path_length);
	return 0;
}

template <typename T> 
int CBinaryTree<T>::learn(T** feature, sample_str& sample_atr)
{
	m_depth_count = 0;
	m_node_count = 0;
	m_node_count++;
	m_binarynode.insert(std::map<int, BinaryNodeType> :: value_type(0, BinaryNodeType()));
	m_binarynode[0].m_node_id[0] = 0;
	sample_atr.subnode[0] = 0;
	sample_atr.subnode[1] = 0;
	m_binarynode[0].learn(feature,sample_atr);
	m_depth_count++;
	m_binarynode[0].m_depth_count = m_depth_count;

	if (!m_binarynode[0].m_isLeaf)
	{
		m_binarynode[0].m_node_id[1] = m_node_count++;
		m_binarynode[0].m_node_id[2] = m_node_count++;
		
		m_depth_count++;
		m_binarynode.insert(std::map<int, BinaryNodeType> :: value_type(m_binarynode[0].m_node_id[1], BinaryNodeType()));
		
		m_binarynode.insert(std::map<int, BinaryNodeType> :: value_type(m_binarynode[0].m_node_id[2], BinaryNodeType()));
		m_binarynode[1].m_node_id[0] = 1;
		m_binarynode[2].m_node_id[0] = 2;
		m_binarynode[1].m_depth_count = m_depth_count;
		m_binarynode[2].m_depth_count = m_depth_count;

		sample_str sample_atributs_left(sample_atr.sample_num,sample_atr.feature_length,sample_atr.subnode[0],1);
		sample_atributs_left.max_tree_depth = sample_atr.max_tree_depth;
		memcpy(sample_atributs_left.p_index,sample_atr.p_sub_index0,sizeof(int)*sample_atr.subnode[0]);
		sample_atributs_left.value_minmax = sample_atr.value_minmax;
		
		sample_str sample_atributs_right(sample_atr.sample_num,sample_atr.feature_length,sample_atr.subnode[1],2);
		sample_atributs_right.max_tree_depth = sample_atr.max_tree_depth;
		memcpy(sample_atributs_right.p_index,sample_atr.p_sub_index1,sizeof(int)*sample_atr.subnode[1]);
		sample_atributs_right.value_minmax = sample_atr.value_minmax;

		recursiveLearn(feature,sample_atributs_left);
		recursiveLearn(feature,sample_atributs_right);

		//m_depth_count++;
	}else
		m_count_leaf++;
	m_isend = true;
	return 0;
}

template <typename T> 
void CBinaryTree<T>::serialize(std::ostream& s)
{
	s << " " << m_node_count;
	s << " " << m_count_leaf;
	for(int i=0;i<m_node_count;i++)
	{
		m_binarynode[i].serialize(s);
	}
}

template <typename T> 
void CBinaryTree<T>::deserialize(std::istream& s)
{
	m_node_count=0;
	s >> m_node_count;
	s >> m_count_leaf;
	if (m_node_count!=0)
	{
		for(int i=0;i<m_node_count;i++)
		{
			m_binarynode[i].deserialize(s);
		}
	}
	
}

template <typename T> 
int CBinaryTree<T>::recursiveLearn(T** feature, sample_str& sample_atr)
{
	m_binarynode[sample_atr.node_id].learn(feature,sample_atr);
	if (!m_binarynode[sample_atr.node_id].m_isLeaf && (m_binarynode[sample_atr.node_id].m_depth_count<sample_atr.max_tree_depth))
	{
		m_binarynode[sample_atr.node_id].m_node_id[1] = m_node_count++;
		m_binarynode[sample_atr.node_id].m_node_id[2] = m_node_count++;
		
		int node_id = m_binarynode[sample_atr.node_id].m_node_id[1]; 
		m_binarynode.insert(std::map<int, BinaryNodeType>::value_type(node_id, BinaryNodeType()));
		m_binarynode[node_id].m_node_id[0] = node_id;
		m_binarynode[node_id].m_depth_count = m_binarynode[sample_atr.node_id].m_depth_count + 1;

		node_id = m_binarynode[sample_atr.node_id].m_node_id[2];
		m_binarynode.insert(std::map<int, BinaryNodeType>::value_type(node_id, BinaryNodeType()));
		m_binarynode[node_id].m_node_id[0] = node_id;
		m_binarynode[node_id].m_depth_count = m_binarynode[sample_atr.node_id].m_depth_count + 1;

		sample_str sample_atributs_left(sample_atr.sample_num,sample_atr.feature_length,sample_atr.subnode[0],m_binarynode[sample_atr.node_id].m_node_id[1]);
		sample_atributs_left.max_tree_depth = sample_atr.max_tree_depth;
		memcpy(sample_atributs_left.p_index,sample_atr.p_sub_index0,sizeof(int)*sample_atr.subnode[0]);
		sample_atributs_left.value_minmax = sample_atr.value_minmax;
		
		sample_str sample_atributs_right(sample_atr.sample_num,sample_atr.feature_length,sample_atr.subnode[1],m_binarynode[sample_atr.node_id].m_node_id[2]);
		sample_atributs_right.max_tree_depth = sample_atr.max_tree_depth;
		memcpy(sample_atributs_right.p_index,sample_atr.p_sub_index1,sizeof(int)*sample_atr.subnode[1]);
		sample_atributs_right.value_minmax = sample_atr.value_minmax;		

		recursiveLearn(feature,sample_atributs_left);
		recursiveLearn(feature,sample_atributs_right);
	}else
		m_count_leaf++;

	return 0;
}

/// A isolation forest.
template<typename T>
class CIsolationForest {
public:
    typedef CBinaryTree<T> BinaryTreeType;
	CIsolationForest();
	~CIsolationForest();
	 void deserialize(std::istream&);  
	int learn( T ** feature, sample_str& sample_atr,int tree_num=256);   
    void serialize(std::ostream&);
	int predict(T* feature,float& result);
private:
	int m_tree_num;
	int m_train_num;
	float m_cn;
	std::vector<BinaryTreeType> m_binarytree;
};

template <typename T>
CIsolationForest<T>::CIsolationForest()
{
}

template <typename T>
CIsolationForest<T>::~CIsolationForest()
{
}

template <typename T>
int CIsolationForest<T>::predict(T* feature,float& result)
{
	float total_path_length = 0;

	for (int i=0;i<m_tree_num;i++)
	{
		ret_result ret_results = {0,0,0};
		m_binarytree[i].predict(feature,ret_results);
		total_path_length += ret_results.path_length;
	}
	float avg_path = total_path_length/m_tree_num;
	//printf("avg path is %f\n",avg_path);
	result = std::pow(2,-avg_path/m_cn);
	return 0;
}

template <typename T>
int CIsolationForest<T>::learn( T ** feature, sample_str& sample_atr,int tree_num)
{
	min_max minmax;
	minmax.p_min = new float[sample_atr.feature_length]();
	minmax.p_max = new float[sample_atr.feature_length]();
	int mins,maxs;

	for (int i=0;i<sample_atr.feature_length;i++)
	{
		minmax.p_min[i] = feature[0][i];
		minmax.p_max[i] = feature[0][i];
		for (int j=0;j<sample_atr.sample_num;j++)
		{
			minmax.p_min[i] = minmax.p_min[i]<feature[j][i] ? minmax.p_min[i]:feature[j][i];
			minmax.p_max[i] = minmax.p_max[i]>feature[j][i] ? minmax.p_max[i]:feature[j][i];
		}
		minmax.p_min[i] = minmax.p_min[i];
		minmax.p_max[i] = minmax.p_max[i];
	}
	sample_atr.value_minmax = minmax;
	m_binarytree.clear();
	m_tree_num = tree_num;
	m_train_num = sample_atr.train_num;
	sample_atr.max_tree_depth = 0.9+log(1.0*m_train_num)/log(2.0);

	m_cn = 2.0*(log(1.0*(1.0*m_train_num-1)) + 0.5772156649 ) - 2.0*(1.0*m_train_num-1)/(1.0*m_train_num);
	
	int count_zero_node=0;
	for(int i=0;i<m_tree_num;i++)
	{
		BinaryTreeType binarytree;
		sample_atr.get_random_sample();
		binarytree.learn(feature,sample_atr);
		while (binarytree.m_node_count==0)
		{
			sample_atr.get_random_sample();
			binarytree.learn(feature,sample_atr);
			continue;
		}
		if (binarytree.m_node_count==0)
		{
			count_zero_node++;
			continue;
		}
		m_binarytree.push_back(binarytree);
	}

	m_tree_num = m_tree_num-count_zero_node;

	DELETE_ARR_PTR(minmax.p_min);
	DELETE_ARR_PTR(minmax.p_max);
	return 0;
}

template <typename T>
void CIsolationForest<T>::serialize(std::ostream& s)
{
	s << " " << m_tree_num;
	s << " " << m_train_num;
	s << " " << m_cn;
	for(int i=0;i<m_tree_num;i++)
	{
		m_binarytree[i].serialize(s);
	}
}

template <typename T>
void CIsolationForest<T>::deserialize(std::istream& s)
{
	s >> m_tree_num;
	s >> m_train_num;
	s >> m_cn;
	m_binarytree.clear();
	m_binarytree.resize(m_tree_num);
	for(int i=0;i<m_tree_num;i++)
	{
		m_binarytree[i].deserialize(s);
	}
	
}
