// test.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "PCA.h"


int _tmain(int argc, _TCHAR* argv[])
{
	int tableNum = 20, goodNum = 10, badNum = 10; //number of good and bad should be close
	int dataNum = 4096;
	double** DataTable = new double*[tableNum];      //Input Data for building model
	double** DataTableW = new double*[tableNum];     //Input Data after pre-process
	double** spectrumPW = new double*[tableNum];     //power spectrum of Input Data 
	double** spectrumPWdB = new double*[tableNum/2]; //power spectrum in dB. Half only, no repeat data
	double** spectrumPha = new double*[tableNum];    //phase of Input Data
	char** filename = new char*[tableNum];
	FILE** goodFILE = new FILE* [goodNum];
	FILE** badFILE = new FILE* [badNum];
	CPCAModel ModelS; // I want one...
	 

	for(int i=0;i<tableNum;i++)
	{
		DataTable[i] = new double [dataNum];
		DataTableW[i] = new double [dataNum];
		spectrumPW[i] = new double [dataNum];
		spectrumPWdB[i] = new double [dataNum/2];
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


	//----training PCA Model----
	ModelS.Training(spectrumPWdB, tableNum, dataNum/2, NULL, true, "ModelS.bin"); //After training, the model is saved as "ModelS.bin" in your file. You can distribute the bin file to other PC. 
	printf("Done Model Training. ");


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


	
	

	//delete memory and close file
	for(int i=0;i<tableNum;i++)
	{
		delete [] DataTable[i];
		delete [] DataTableW[i];
		delete [] spectrumPW[i];
		delete [] spectrumPha[i];
	}

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

