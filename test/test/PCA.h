#pragma once

#include "linalg.h"
#include "MatrixOp.h"
#include "FFT.h"
#include "math.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef PCAMODEL_EXPORTS    
#define PCAMODEL_EXPORTS   __declspec( dllexport )
#else  
#define PCAMODEL_EXPORTS __declspec(dllimport)   
#endif  

class PCAMODEL_EXPORTS CPCAModel
{
public:
	 CPCAModel(void);
	 ~CPCAModel(void);

	void Training(double** DataSet, int sensorNum, int SampleNum, double confidscore, bool saveModel, char modelName[]);
	short LoadModel(char modelName[]);
	short PCATransfer(double* Indata, double* OutData);
	short GetModelParameter(int* NumofRaw, int* NumofColum, int* numPC, double* confindSet);


private:
	short SetParameter(int sensorNum, int SampleNum, double confidscore);
	void SaveModel(char modelName[]);

	bool ModelLoaded;
	double** Model;
	int numObserve;
	int InDimension;
	double* DimemsionCenter;
	unsigned short PCnum;
	double* score;
	double confidence;

	// head information for save model
	struct headInf
	{
	    char myhead[13];
		char numhead[14]; 
		char dimhead[15];
		char centerhead[17];
		char modelhead[12];
		char pcsNumhead[24];
		char scorehead[6]; 
		char confidhead[11];
	};
	headInf inf;
	
};

#ifdef __cplusplus
}
#endif
