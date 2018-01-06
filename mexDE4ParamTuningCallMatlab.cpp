/**
	MATLAB: [ Coefficients ] = mexDE4ParamTuning( nDim, nPop, ParamMatrix);
	ParamMatrix = [lambda_FC, lambda_FI, lambda_Alpha, lambda_Beta, saliency_ratio, FC_ratio]

	nDim is parameter's number
	nPop is iter's number
	ParamMatrix is [min, max; min, max; min, max]
*/

#include "mex.h"
#include "./DE/DESolver.h"
#include "./DE/DESolver.cpp"

using namespace std;

#define g_nDEMaxGeneration	50
#define g_dDEScaleF			0.9
#define g_dDECrossRatio		1.0

// Parameters training by DE4ParamTuningSolver
class DE4ParamTuningSolver : public DESolver
{
	 
public:
	DE4ParamTuningSolver(int dim, int pop, double fusing_th_) : DESolver(dim,pop,fusing_th_), count(0){;}
	double EnergyFunction(double *pTrialSolution, bool &bAtSolution);
    double EnergyFunction1(double *pTrialSolution);

private:
	int count;

};

/*
double DE4ParamTuningSolver::EnergyFunction(double *pTrialSolution, bool &bAtSolution)
{
	double para0 = pTrialSolution[0];
	double para1 = pTrialSolution[1];
	double para2 = pTrialSolution[2];
	//int para3 = (int)pTrialSolution[3];
	//int para4 = (int)pTrialSolution[4];
    //int para5 = (int)pTrialSolution[5];
	//int para6 = (int)pTrialSolution[6];
	   
	mxArray *parm[3];
	mxArray *lhs;

	parm[0] = mxCreateDoubleScalar(para0);
	parm[1] = mxCreateDoubleScalar(para1);
	parm[2] = mxCreateDoubleScalar(para2);
	//parm[3] = mxCreateDoubleScalar(para3);
	//parm[4] = mxCreateDoubleScalar(para4);
	//parm[5] = mxCreateDoubleScalar(para5);
	//parm[6] = mxCreateDoubleScalar(para6);

    mexCallMATLAB(1, &lhs, 3, parm, "EnergyFunction");
	
	double err = (double)(*mxGetPr(lhs));
	//mxFree(lhs);
    //mxFree(parm[6]);
    //mxFree(parm[5]);
    //mxFree(parm[4]);
    //mxFree(parm[3]);
    //mxFree(parm[2]);
    //mxFree(parm[1]);
    //mxFree(parm[0]);

	mexPrintf("%f \n", (double)err);

   	return (double)err;

}*/

double DE4ParamTuningSolver::EnergyFunction(double *pTrialSolution, bool &bAtSolution)
{
	//double para0 = pTrialSolution[0];
	//double para1 = pTrialSolution[1];
	//double para2 = pTrialSolution[2];
	//int para3 = (int)pTrialSolution[3];
	//int para4 = (int)pTrialSolution[4];
    //int para5 = (int)pTrialSolution[5];
	//int para6 = (int)pTrialSolution[6];
	   
	mxArray **parm=(mxArray**)malloc(nDim*sizeof(mxArray*));
	mxArray *lhs;

	// parm[0] = mxCreateDoubleScalar(para0);
	// parm[1] = mxCreateDoubleScalar(para1);
	// parm[2] = mxCreateDoubleScalar(para2);
	for(int i=0; i<nDim; i++){
		if(isInt[i]==0)
			parm[i]=mxCreateDoubleScalar(double(pTrialSolution[i]));
		else
			parm[i]=mxCreateDoubleScalar(double(int(pTrialSolution[i])));
	}

    mexCallMATLAB(1, &lhs, nDim, parm, "EnergyFunction");
	
	double err = (double)(*mxGetPr(lhs));
	//mxFree(lhs);
    //mxFree(parm[6]);
    for(int i=0; i<nDim; i++){
		mxDestroyArray(parm[i]);
	}
	free(parm);
	parm=NULL;
	//mxDestroyArray(parm);

	mexPrintf("%f \n", (double)err);
	

   	return (double)err;

}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
	// initialize
	double *pnDim = mxGetPr(prhs[0]);
	double *pnPop = mxGetPr(prhs[1]);
	double *p_fusing_th=mxGetPr(prhs[2]);
	int nD = (int)(*pnDim);
	int nP = (int)(*pnPop);
	double fusing_th=(double)(*p_fusing_th);

	DE4ParamTuningSolver deSolver(nD, nP, fusing_th);

	// setup DE
	double *pdMin = new double[nD];
	double *pdMax = new double[nD];
	double *pIsInt= new double[nD];

	double* pDatParamMatrix = mxGetPr(prhs[3]);
	for (int i=0;i<nD;i++)
	{
		pdMin[i] = (double)pDatParamMatrix[0*nD+i];
		pdMax[i] = (double)pDatParamMatrix[1*nD+i];
		pIsInt[i]= (double)pDatParamMatrix[2*nD+i];
	}

	deSolver.Setup(pdMin, pdMax, pIsInt, stBest1Exp, g_dDEScaleF, g_dDECrossRatio);
	// show DE param
	deSolver.ShowDEParameters();

	// setup DE constraints (no discrete constraints)
	bool *pbRange = new bool[nD];
	for (int i=0; i<nD; i++)
	{
		pbRange[i] = true;
	}

	deSolver.SetupRangeConstraints(pdMin, pdMax, pbRange);

	// run DE
	mexPrintf("Calculating...\n\n");
    double bestEner=deSolver.Solve(g_nDEMaxGeneration);
	mexPrintf("Calculating  end \n");
	double *solution = deSolver.Solution();
	mexPrintf("\nBest Coefficients:\n");

// 	//output Coefficients
// 	plhs[0] = mxCreateDoubleMatrix(1, nD, mxREAL);
// 	double *pMatCoefficients = (double*)mxGetPr(plhs[0]);
// 
	for (int i=0;i<nD;i++)
	{
		//pMatCoefficients[i*nD] = solution[i];
		mexPrintf("[%d]: %f, ",i, pIsInt[i] ? int(solution[i]) : solution[i]);
	}
	mexPrintf("\n");


    //double x=deSolver.EnergyFunction1(solution);


	delete [] pdMin;
	delete [] pdMax;
	delete [] pbRange;

}
