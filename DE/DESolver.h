// Differential Evolution Solver Class
// Based on algorithms developed by Dr. Rainer Storn & Kenneth Price
// Written By: Lester E. Godwin
//             PushCorp, Inc.
//             Dallas, Texas
//             972-840-0208 x102
//             godwin@pushcorp.com
// Created: 6/8/98
// Last Modified: 6/8/98
// Revision: 1.0
// -----------------------------------------------------------------------
// Extended to handle: (1) InRange constraint;
//                     (2) Integer variables
//                     (3) Discrete variables
// For constraints of any form, revise the objective function by regularization
// Extended by: Wei Feng
// Last Modified: Nov. 10, 2010
//

#if !defined(_DESOLVER_H)
#define _DESOLVER_H

#define stBest1Exp			0
#define stRand1Exp			1
#define stRandToBest1Exp	2
#define stBest2Exp			3
#define stRand2Exp			4
#define stBest1Bin			5
#define stRand1Bin			6
#define stRandToBest1Bin	7
#define stBest2Bin			8
#define stRand2Bin			9

class DESolver;

typedef void (DESolver::*StrategyFunction)(int);

class DESolver
{
public:
	DESolver(int dim,int popSize, double fusing_th_);
	~DESolver(void);
	
	// Setup() must be called before solve to set min, max, strategy etc.
	void Setup(double min[],double max[], double pIsInt[], int deStrategy,
							double diffScale,double crossoverProb);

	//void Setup2(double min[],double max[],int deStrategy,
	//	double diffScale,double crossoverProb,double*iniCandidate);

	// SetupXConstraints() must be called before solve to set necessary constraints
	void SetupRangeConstraints(double *min,double *max,bool *inrange);
	void SetupIntegerConstraints(bool *isint);
	// If a variable is Discrete, it must also be Integer and InRange.
	// So, when setting up Discrete constraints, the InRange and Integer
	//     constraints should also be set up correctly.
	void SetupDiscreteConstraints(bool*isdiscrete,double**discvaltab,int*discvalnum=0);

	// Solve() returns true if EnergyFunction() returns true.
	// Otherwise it runs maxGenerations generations and returns false.
    virtual double Solve(int maxGenerations);

	// EnergyFunction must be overridden for problem to solve
	// testSolution[] is nDim array for a candidate solution
	// setting bAtSolution = true indicates solution is found
	// and Solve() immediately returns true.
	virtual double EnergyFunction(double testSolution[],bool &bAtSolution) = 0;
	
	// Info
	int Dimension(void) { return(nDim); }
	int Population(void) { return(nPop); }
	void ShowDEParameters(void);

	// Call these functions after Solve() to get results.
	double Energy(void) { return(bestEnergy); }
	double *Solution(void) { return(bestSolution); }

	int Generations(void) { return(generations); }

protected:
	void SelectSamples(int candidate,int *r1,int *r2=0,int *r3=0,
									 int *r4=0,int *r5=0);
	double RandomUniform(double min,double max);

	int nDim;
	int nPop;
	int generations;

	double * initialCandidate;

	int strategy;
	StrategyFunction calcTrialSolution;
	double scale;
	double probability;
	double fusing_th;
	double trialEnergy;
	double bestEnergy;

	double *trialSolution;	
	double *bestSolution;
	double *popEnergy;
	double *population;

	// Constraints related (memory new/delete outside)
	double *isInt;
	bool *bInRange;
	double *dMin;
	double *dMax;
	bool *bIsInteger;
	bool *bIsDiscrete;	
	int    *nDiscreteValNum;
	double **dDiscreteValTab;
	bool bRangeConstrained;
	bool bIntegerConstrained;
	bool bDiscreteConstrained;

	void ForceSolutionInRange(double *solution);
	void GetEvaluteSolution(double *evasolution,double *insolution);

private:
	void Best1Exp(int candidate);
	void Rand1Exp(int candidate);
	void RandToBest1Exp(int candidate);
	void Best2Exp(int candidate);
	void Rand2Exp(int candidate);
	void Best1Bin(int candidate);
	void Rand1Bin(int candidate);
	void RandToBest1Bin(int candidate);
	void Best2Bin(int candidate);
	void Rand2Bin(int candidate);
};

/*  Calling Procedure

	Step 1. Define a new class MySolver inherited from DESolver, with EnergyFunction() defined;
	Step 2. Calling MySolver.Setup() to configure the solver;
	Step 3. Calling MySolver.SetupRangeConstraints() and/or
	                MySolver.SetupIntegerConstraints() and/or 
					MySolver.SetupDiscreteConstraints()
			to correctly setup necessary constriants;
	Step 4. Calling MySolver.ShowDEParameters() to show the configuration of DE [optional];
	Step 5. Calling MySolver.Solve() to run the DE algorithm;
	Step 6. Calling MySolver.Solution() to retrieve the resultant solution;

*/

#endif // _DESOLVER_H
