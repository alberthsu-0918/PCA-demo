// PCAModel.cpp : Defines the exported functions for the DLL application.
//

#include "StdAfx.h"
#include "PCA.h"


using namespace alglib;


CPCAModel::CPCAModel(void)
{
	this->ModelLoaded = false;
	this->numObserve = -1;
	this->InDimension = -1;
	this->Model = new double*[1];
	this->Model[0] = NULL;
	this->DimemsionCenter = NULL;
	this->confidence = 80.0;

	char head[] = "Albert's PCA";
	memcpy(this->inf.myhead, head, sizeof(head));
	char confidhead[] = "Confidence";
	memcpy(this->inf.confidhead, confidhead, sizeof(confidhead));
	char numhead[] = "NumOfObserver"; 
	memcpy(this->inf.numhead, numhead, sizeof(numhead));
	char dimhead[] = "NumOfDimension";
	memcpy(this->inf.dimhead, dimhead, sizeof(dimhead));
	char centerhead[] = "CoordinateCenter";
	memcpy(this->inf.centerhead, centerhead, sizeof(centerhead));
	char modelhead[] = "ModelMatrix";
	memcpy(this->inf.modelhead, modelhead, sizeof(modelhead));
	char pcsNumhead[] = "NumOfPrincipalComponent";
	memcpy(this->inf.pcsNumhead, pcsNumhead, sizeof(pcsNumhead));
	char scorehead[] = "Score";
	memcpy(this->inf.scorehead, scorehead, sizeof(scorehead));
}

CPCAModel::~CPCAModel(void)
{
	
}



short CPCAModel::SetParameter(int sensorNum, int SampleNum, double confidscore)
{
	if(confidscore<0 || confidscore>=100 )
	{
		printf("Range of Confidence has to be 0 < ~ < 100%");
		return -1;
	}

	if(sensorNum<1 || SampleNum<1)
	{
		printf("Be careful!");
		return -1;
	}

	//check if Model memory space is used
	if(this->Model[0] != NULL)
	{
		if(this->numObserve>0 && this->InDimension>0)
		{
			for(int i=0;i<numObserve;i++)
			{
				free(this->Model[i]);
			}
			free(this->Model);
			free(this->DimemsionCenter);
			free(this->score);
		}
		else
		{
			printf("Model is created, but dimension information <0\n");
			return -2;
		}
	}
	else
	{
		free(this->Model[0]);
		free(this->Model);
	}
	
	if(confidscore != NULL)
	{
		this->confidence = confidscore;
	}
	this->numObserve = sensorNum;
	this->InDimension = SampleNum;

	this->DimemsionCenter = new double [SampleNum];
	this->Model = new double* [SampleNum];
	for(int i=0;i< SampleNum;i++)
	{
		this->Model[i] = new double [SampleNum]; 
	}
	this->score = new double [SampleNum];

	return 0;
}

void CPCAModel::SaveModel(char modelName[])
{
	FILE* modelfile = fopen(modelName, "wb");
	
	fwrite(this->inf.myhead, 1, sizeof(this->inf.myhead)/sizeof(char) , modelfile);

	fwrite(this->inf.confidhead, 1, sizeof(this->inf.confidhead)/sizeof(char), modelfile);
	fwrite(&this->confidence, sizeof(double), 1, modelfile);

	fwrite(this->inf.numhead, 1, sizeof(this->inf.numhead)/sizeof(char), modelfile);
	fwrite(&this->numObserve, sizeof(int), 1, modelfile);

	fwrite(this->inf.dimhead, 1, sizeof(this->inf.dimhead)/sizeof(char), modelfile);
	fwrite(&this->InDimension, sizeof(int), 1, modelfile);

	fwrite(this->inf.centerhead, 1, sizeof(this->inf.centerhead)/sizeof(char), modelfile);
	fwrite(this->DimemsionCenter, sizeof(double), this->InDimension, modelfile);

	fwrite(this->inf.modelhead, 1, sizeof(this->inf.modelhead)/sizeof(char), modelfile);
	for(int i=0;i<this->InDimension;i++)
	{
		fwrite(this->Model[i], sizeof(double), this->InDimension, modelfile);
	}

	fwrite(this->inf.pcsNumhead, 1, sizeof(this->inf.pcsNumhead), modelfile);
	fwrite(&this->PCnum, sizeof(unsigned short), 1, modelfile);

	fwrite(this->inf.scorehead, 1, sizeof(this->inf.scorehead)/sizeof(char), modelfile);
	fwrite(this->score, sizeof(double), this->InDimension, modelfile);

	fclose(modelfile);
}

short CPCAModel::LoadModel(char modelName[])
{
	short err = -1;
	int sensorNum = -1, sampleNum = -1;
	long offset = 0;
	double confidence;
	double* dimCenter;
	FILE* svfile = fopen(modelName, "rb");
	if(svfile == NULL)
	{
		printf("invalid model\n");
		err = -1;
		return err;
	}

	offset = offset + sizeof(this->inf.myhead) + sizeof(this->inf.confidhead);
	fseek(svfile, offset, SEEK_SET);
	fread(&confidence, sizeof(double), 1, svfile);

	offset = offset + 1*sizeof(double) + sizeof(this->inf.numhead);
	fseek(svfile, offset, SEEK_SET);
	fread(&sensorNum, sizeof(int), 1, svfile);

	offset = offset + 1*sizeof(int) + sizeof(this->inf.dimhead);
	fseek(svfile, offset, SEEK_SET);
	fread(&sampleNum, sizeof(int), 1, svfile);

	err = SetParameter(sensorNum, sampleNum, confidence);
	if(err!=0)
	{
		printf("SetParameter fail\n");
		return err;
	}

	offset = offset + 1*sizeof(int) + sizeof(this->inf.centerhead);
	fseek(svfile, offset, SEEK_SET);
	fread(this->DimemsionCenter, sizeof(double), sampleNum, svfile);

	offset = offset + sampleNum*sizeof(double) + sizeof(this->inf.modelhead);
	fseek(svfile, offset, SEEK_SET);
	for(int i=0; i<sampleNum; i++)
	{
		fread(this->Model[i], sizeof(double), sampleNum, svfile);
	}

	offset = offset + sampleNum*sampleNum*sizeof(double) + sizeof(this->inf.pcsNumhead);
	fseek(svfile, offset, SEEK_SET);
	fread(&this->PCnum, sizeof(unsigned short), 1, svfile);

	offset = offset + 1*sizeof(unsigned short) + sizeof(this->inf.scorehead);
	fseek(svfile, offset, SEEK_SET);
	fread(this->score, sizeof(double), sampleNum, svfile);
	/*
	fseek(svfile, sizeof(head) + sizeof(numhead) + 1*sizeof(int) + sizeof(dimhead) + 1*sizeof(int) + sizeof(centerhead)+sizeof(double)*2048+sizeof(modelhead)+sizeof(double)*2048*2047 , SEEK_SET);
	double ob[2048];
	fread(ob, sizeof(double), 2048, testRead);
	 */

	fclose(svfile);
	this->ModelLoaded = true;
	
	return err;
}

//=================== Parameter of Training() ===================
// double** DataSet(In):   sensor_1: /data1_0 data1_1 data1_2 ........data1_n\
//                         sensor_2: |data2_0 data2_1 data2_2 ........data2_n|
//                         ...       | ...                                   |
//                         sensor_m  \datam_0 datam_1 datam_2 ........datam_n/
// int sennsorNum:         number of raw of your input matrix, DataSet. Take above DataSet as an example, int sennsorNum = m
// int SampleNum:          number of column of your input matrix, DataSet. Take above DataSet as an example, int sennsorNum = n. 
// double confidscore:     percentage of confidence you want. Default = 80.
// bool saveModel:         true or false. see if you want to save trained model to bin file, that you can load it in the future without training
// car ModelName[]:        filename, ex: "ModelS.bin".
void CPCAModel::Training(double** DataSet, int sensorNum, int SampleNum, double confidscore, bool saveModel, char modelName[])
{

	short err = SetParameter(sensorNum, SampleNum, confidscore);
	if(err!=0)
	{
		printf("Fail to set Parameter\n");
		return;
	}

	if( (this->Model[0] == NULL) || (this->DimemsionCenter == NULL) || this->numObserve<0 || this->InDimension<0)
	{
		printf("Error. Please reset parameters by SetParameter()\n");
		return;
	}
	else
	{
		int rank = 0;
		double*  eigval;
		double** Data = new double* [this->numObserve];
		double sum =0.0, *mean = new double [this->InDimension];
		for(int i=0;i<this->numObserve;i++)
		{
			Data[i] = new double [this->InDimension];
			memcpy(Data[i], DataSet[i], sizeof(double)*(this->InDimension));
		}


		//---- Centralize data. mean value of Data in each dimension has to be 0 ----
		for(int i=0; i<this->InDimension;i++)
		{
			sum = 0.0;
			for(int j=0; j < this->numObserve;j++)
			{
				sum = sum + DataSet[j][i];
			}
			mean[i] = sum / this->numObserve;
			this->DimemsionCenter[i] = mean[i];

			for(int j=0; j < this->numObserve;j++)
			{
				Data[j][i] = Data[j][i] - mean[i];
			}
		}

		//test save data before SVD
		/*
		FILE* testtttf = fopen("test0.bin", "wb");
		for(int i=0;i<this->numObserve;i++)
		{
			fwrite(Data[i], sizeof(double), this->InDimension, testtttf);
		}
		fclose(testtttf);
		*/

		//---- Using SVD to build transfer model, Vt' ---- 
		double* a = new double [this->numObserve*this->InDimension]; //input
		double* v; //V'

		for(int i=0;i<this->numObserve;i++)
		{
			memcpy(a+i*this->InDimension, Data[i], sizeof(double)*this->InDimension);
		}
		
		real_2d_array A, U, Vt;
		real_1d_array W;
		A.setcontent(this->numObserve, this->InDimension, a);
		U.setlength(this->numObserve, this->numObserve);
		Vt.setlength(this->InDimension, this->InDimension);
		int svdok = rmatrixsvd(A, this->numObserve, this->InDimension, 2, 2, 2, W, U, Vt);
		rank = W.length();
		eigval = W.getcontent();
		v = Vt[0];
		
		//transpose Vt and save to Model
		for(int i=0;i<this->InDimension;i++)
		{
			for(int j=0;j<this->InDimension;j++)
			{
				this->Model[j][i] = v[i*this->InDimension+j];
			}
		}
		
		//test: see model and singular value
		/*
		for(int i=0;i<20;i++)
		{
			for(int j=0;j<20;j++)
				printf("%.4lf ", this->Model[i][j]);
			printf("\n");
		}

		for(int i=0;i<20;i++)
		{
			printf("%.4lf \n", eigval[i]);
		}
		
		for(int i=0;i<20;i++)
		{
			for(int j=0;j<20;j++)
				printf("%.4lf ", U[i][j]);
			printf("\n");
		}
		*/
		//---- Record eigenvalues(variance) as score, and find number of principal componment ---- 
		double pcSum = 0.0, tempSum =0.0;
		int pcsNum=0;
		memset(this->score, 0, sizeof(double)*(this->InDimension));
		memcpy(this->score, eigval, sizeof(double)*rank);
		
		for(int i=0; i<this->InDimension; i++)
		{
			pcSum = pcSum + score[i]*score[i];
		}

		do{
			tempSum = tempSum + score[pcsNum]*score[pcsNum]; 
			pcsNum++;
		}while(tempSum < (this->confidence*pcSum/100));
		
		this->PCnum = pcsNum;
		

		//---- Save Model, if savemodel == true ----
		if(saveModel)
		{
			SaveModel(modelName);	
		}

		this->ModelLoaded = true;

		//---test--- see transOut, result after project to a new coordinate
		/*
		double** transOut = new double* [this->numObserve];
		for(int i=0;i<this->numObserve; i++)
		{
			transOut[i] = new double [this->InDimension];
		}
		int err = multiplyM(Data, this->numObserve, this->InDimension, this->Model, this->InDimension, this->InDimension, transOut);
		printf("\n\n Check final result:\n\n");
		for(int i=0;i<20;i++)
		{
			for(int j=0;j<20;j++)
			{
				printf("%.4lf ", transOut[i][j]);
			}
			printf("\n");
		}
		*/
		//---end of the test---

		

		//---- free memory ----
		for(int i=0; i<this->numObserve;i++)
		{
			delete [] Data[i];
		}

		delete [] Data;
		delete [] a;
		
		
	}
}

short CPCAModel::PCATransfer(double* indata, double* outdata)
{
	short err = -1;
	double** datain = new double* [1], **dataout = new double* [1];
	

	if(this->ModelLoaded)
	{
		datain[0] = new double [this->InDimension];
		dataout[0] = new double [this->PCnum];

		memcpy(datain[0], indata, sizeof(double)*this->InDimension );

		//move input data to model's central coordinate
		for(int i=0;i<this->InDimension;i++)
		{
			datain[0][i] = datain[0][i] - this->DimemsionCenter[i];
		}

		err = multiplyM(datain, 1, this->InDimension, this->Model, this->InDimension, this->PCnum, dataout);
		memcpy(outdata, dataout[0], sizeof(double)*this->PCnum);

		free(datain[0]);
		free(dataout[0]);
	}
	else
	{
		printf("No valid model is loaded\n");
		return err;
	}

	free(datain);
	free(dataout);

	return err;
}

short CPCAModel::GetModelParameter(int* NumofRaw, int* NumofColum, int* numPC, double* confindSet)
{
	short err = -1;

	if(this->ModelLoaded)
	{
		*NumofRaw = this->numObserve;
		*NumofColum = this->InDimension;
		*numPC = this->PCnum;
		*confindSet = this->confidence;
		err = 0;
	}
	else
	{
		printf("No Model is loaded\n");
		return -1;
	}

	return err;
}


