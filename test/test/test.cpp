// test.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "PCA.h"

double calculateSD(double data[], int dataLength);

int _tmain(int argc, _TCHAR* argv[])
{
	int tableNum = 20, goodNum = 10, badNum = 10; //number of good and bad should be close
	int dataNum = 4096;
	double** DataTable = new double*[tableNum];      //Input Data for building model
	double** DataTableW = new double*[tableNum];     //Input Data after pre-process
	double** spectrumPW = new double*[tableNum];     //power spectrum of Input Data 
	double** spectrumPWdB = new double*[tableNum];   //power spectrum (dB) of Input Data 
	double** spectrumPWdBN = new double*[tableNum];  //normalized power spectrum (dB) of Input Data 
	double** spectrumPha = new double*[tableNum];    //phase of Input Data
	double*  STD = new double [dataNum / 2];
	double** temp = new double*[dataNum / 2];
	char** filename = new char*[tableNum];
	FILE** goodFILE = new FILE* [goodNum];
	FILE** badFILE = new FILE* [badNum];
	CPCAModel ModelS; // I want one...
	CPCAModel Model3; 

	for(int i=0;i<tableNum;i++)
	{
		DataTable[i] = new double [dataNum];
		DataTableW[i] = new double [dataNum];
		spectrumPW[i] = new double [dataNum];
		spectrumPWdB[i] = new double [dataNum/2];
		spectrumPWdBN[i] = new double[dataNum / 2];
		spectrumPha[i] = new double [dataNum]; 
		filename[i] = new char [15];
	}

	for(int i=0;i<goodNum;i++)
	{
		sprintf(filename[i], "g%d.csv", i+1);
		goodFILE[i] = fopen(filename[i], "r");
	}

	for(int i=0;i<badNum;i++)
	{
		sprintf(filename[i+goodNum], "bfErr%d.csv", i+1);
		badFILE[i] = fopen(filename[i+goodNum], "r");
	}

	//----Load data from file----
	printf("Load Data from files. ");
	printf("1~10 are good sampes. 11~20 are bad samples\n");
	int index = 0;
	for(int i=0;i<goodNum;i++)
	{
		index = 0;
		while((fscanf(goodFILE[i],"%lf", DataTable[i]+index)) !=EOF){
			index++;
		}
	}

	for(int i=0;i<badNum;i++)
	{
		index = 0;
		while((fscanf(badFILE[i],"%lf", DataTable[i+goodNum]+index)) !=EOF){
			index++;
		}
	}

	//----Pre process Data. Transfer Data to freq domain in dB unit----
	for(int i=0; i<tableNum;i++)
	{
		hanning(dataNum, DataTable[i], DataTableW[i]);
		fft_spectrum_double(dataNum, DataTableW[i], 0, spectrumPW[i], spectrumPha[i]);
		for(int j=0;j<dataNum/2;j++)
		{
			spectrumPWdB[i][j] = 10*log10(spectrumPW[i][j]);
		}
	}


	//----training PCA Model without normalization----
	ModelS.Training(spectrumPWdB, tableNum, dataNum/2, NULL, true, "ModelS.bin"); //After training, the model is saved as "ModelS.bin" in your file. You can distribute the bin file to other PC. 
	printf("Done ModelS Training. \n");


	//----training PCA Model with normalization----
	for (int i = 0; i < dataNum / 2; i++)
	{
		temp[i] = new double[tableNum];
	}

	for (int i = 0; i < tableNum; i++)
	{
		memcpy(spectrumPWdBN[i], spectrumPWdB[i], (dataNum / 2)*sizeof(double));

		for (int j = 0; j < dataNum / 2; j++)
		{
			temp[j][i] = spectrumPWdBN[i][j];
		}
	}

	for (int i = 0; i < dataNum / 2; i++)
	{
		STD[i] = calculateSD(temp[i], tableNum);
	}

	for (int i = 0; i < dataNum / 2; i++)
	{
		for (int j = 0; j < tableNum; j++)
		{
			spectrumPWdBN[j][i] = spectrumPWdBN[j][i] / STD[i];
		}
	}

	Model3.Training(spectrumPWdBN, tableNum, dataNum / 2, NULL, true, "Model3.bin");
	printf("Done Model3 Training. \n");


	//----Test trainging Model by loading and transfer function----
	int a, b, c, sampleRead;
	double d;
	double* plotOut;

	ModelS.LoadModel("ModelS.bin"); //Model is loaded now
	printf("Done Model Loading\n\n");
	ModelS.GetModelParameter(&a, &b, &c, &d);
	plotOut = new double [c];

	sampleRead = a;
	printf("NewCord    1       2\n");
	for(int i=0;i<sampleRead;i++)
	{
	    ModelS.PCATransfer(spectrumPWdB[i], plotOut); //Input data you record to the model. The plotOut is output 
		if(i<10)
			printf("g%d    : ", i);
		else
			printf("bad%d : ", i);
		printf("%.4lf  %.4lf\n", plotOut[0], plotOut[1]);
	}


	//----Test trainging Model3 by loading and transfer function----
	int a3, b3, c3, sampleRead3;
	double d3;
	double* plotOut3;

	Model3.LoadModel("Model3.bin"); //Model is loaded now
	printf("Done Model3 Loading\n\n");
	Model3.GetModelParameter(&a3, &b3, &c3, &d3);
	plotOut3 = new double[c3];

	sampleRead3 = a3;
	printf("NewCord    1       2\n");
	for (int i = 0; i < sampleRead3; i++)
	{
		Model3.PCATransfer(spectrumPWdBN[i], plotOut3); //Input data you record to the model. The plotOut is output 
		if (i < 10)
			printf("g%d    : ", i);
		else
			printf("bad%d : ", i);
		printf("%.4lf  %.4lf\n", plotOut3[0], plotOut3[1]);
	}


	
	

	//delete memory and close file
	for(int i=0;i<tableNum;i++)
	{
		delete [] DataTable[i];
		delete [] DataTableW[i];
		delete [] spectrumPW[i];
		delete [] spectrumPha[i];
		delete [] spectrumPWdB[i];
		delete [] spectrumPWdBN[i];
	}

	for (int i = 0; i < dataNum / 2; i++)
	{
		delete[] temp[i];
	}

	delete[] STD;

	for(int i=0;i<goodNum;i++)
	{
		fclose(goodFILE[i]);
	}
	for(int i=0;i<badNum;i++)
	{
		fclose(badFILE[i]);
	}
	
	system("pause");
	
	return 0;
}

double calculateSD(double data[], int dataLength)
{
	float sum = 0.0, mean, standardDeviation = 0.0;
	int i;

	for (i = 0; i < dataLength; ++i)
	{
		sum += data[i];
	}

	mean = sum / dataLength;

	for (i = 0; i < dataLength; ++i)
		standardDeviation += pow(data[i] - mean, 2);

	return sqrt(standardDeviation / 10);
}

