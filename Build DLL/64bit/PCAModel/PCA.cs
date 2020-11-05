using System.Runtime.InteropServices;
using System;

public class PCA
{
    [DllImport("PCAModel.dll")]
    public static extern IntPtr CPCAModel_Create();
    [DllImport("PCAModel.dll")]
    public static extern void LoadModel_Training(IntPtr pCPCAModel, double[][] DataSet, int sensorNum, int SampleNum, double confidscore, bool saveModel, char[] modelName);
    [DllImport("PCAModel.dll")]
    public static extern short CPCAModel_LoadModel(IntPtr pCPCAModel, char[] modelName);
    [DllImport("PCAModel.dll")]
    public static extern short PCATransfer(double[] Indata, IntPtr OutData);
    [DllImport("PCAModel.dll")]
    public static extern short CPCAModel_GetParameters(IntPtr pCPCAModel, IntPtr NumofRaw, IntPtr NumofColum, IntPtr numPC, IntPtr confindSet);
    [DllImport("PCAModel.dll")]
    public static extern void CPCAModel_Delete(IntPtr pCPCAModel);
}
