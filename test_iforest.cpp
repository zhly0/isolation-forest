#include <iostream>
#include <fstream>
#include <math.h>
#include <cstdlib>
#include "isolation_forest.h"
using namespace std;

int save_model()
{
    int sample_num = 256;
	int feature_length = 180;
    
    float** sins = new float*[sample_num]();
	for (int i=0;i<sample_num;i++)
	{
		sins[i] = new float[feature_length]();
	}
	std::random_device rd0 ; 
	std::default_random_engine e0(rd0) ;
	std::uniform_real_distribution<> u0(-0.5,0.5);
	for (int i=0;i<sample_num;i++)
	{
		for (int j=0;j<feature_length;j++)
		{
			sins[i][j] = u0(e0)+sin(j*3.14/180);
		}
	}
	int train_num = 128;
	sample_str str_sample(sample_num,feature_length,train_num,0);

	CIsolationForest<float> iforest;
	int tree_num = 100;

	iforest.learn(sins,str_sample,tree_num);
	ofstream s("model1.data",ios::out);
	iforest.serialize(s);
	s.close();
    return 0;
}

int test_model()
{
	int sample_num = 256;
	int feature_length = 180;
	
	CIsolationForest<float> iforest_predict;
	ifstream ss("model1.data",ios::in);
	iforest_predict.deserialize(ss);
	float* floats = new float[feature_length]();
	for (int i=0;i<feature_length;i++)
	{
		floats[i] =sin(i*3.14/180);//sins[0][i];//sin(i*3.14/180);
		//printf("result is %f\n",floats[i] );
	}
	float results;
	iforest_predict.predict(floats,results);
	printf("result is %f\n",results);
    return 0;
}

void main()
{
    //step one :save model:
    save_model();
    //step two : test model
    test_model();
}
