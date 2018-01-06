#include "mex.h"
#include <memory.h>
#include "DESolver.h"

#define Element(a,b,c)  a[b*nDim+c]
#define RowVector(a,b)  (&a[b*nDim])
#define CopyVector(a,b) memcpy((a),(b),nDim*sizeof(double))
#define IntVar(a)  ((int)((double)a+0.5))

DESolver::DESolver(int dim,int popSize, double fusing_th_=-9.5) :
					nDim(dim), nPop(popSize), fusing_th(fusing_th_),
					generations(0), strategy(stRand1Exp),
					scale(0.7), probability(0.5), bestEnergy(0.0),
					trialSolution(0), bestSolution(0),
					popEnergy(0), population(0)
{
	trialSolution = new double[nDim];
	bestSolution  = new double[nDim];
	popEnergy	  = new double[nPop];
	population	  = new double[nPop * nDim];
	//initialCandidate = NULL;

	bRangeConstrained = false;
	bIntegerConstrained = false;
	bDiscreteConstrained = false;

	return;
}

DESolver::~DESolver(void)
{
// 	if (trialSolution) delete trialSolution;
// 	if (bestSolution) delete bestSolution;
// 	if (popEnergy) delete popEnergy;
// 	if (population) delete population;
// 	if (initialCandidate) delete initialCandidate;

	//trialSolution = bestSolution = popEnergy = population = 0;
    delete []trialSolution;
    delete []bestSolution;
    delete []popEnergy;
    delete []population;
	delete []isInt;


    dMin = NULL;
    dMax = NULL;
    nDiscreteValNum = NULL;
    dDiscreteValTab = NULL;
    bInRange = NULL;
    bIsInteger = NULL;
    bIsDiscrete = NULL;

    //delete initialCandidate;
	return;
}

void DESolver::Setup(double *min,double *max, double *pIsInt,
						int deStrategy,double diffScale,double crossoverProb)
{
	int i;

	strategy	= deStrategy;
	scale		= diffScale;
	probability = crossoverProb;	
	
	isInt=new double[nDim];
	for(int i=0; i<nDim; i++){
		isInt[i]=pIsInt[i];
	}
	
	for (i=0; i < nPop; i++)
	{
		for (int j=0; j < nDim; j++)
			Element(population,i,j) = RandomUniform(min[j],max[j]);

		popEnergy[i] = 1.0E20;
	}

	for (i=0; i < nDim; i++)
		bestSolution[i] = 0.0;

	switch (strategy)
	{
		case stBest1Exp:
			//calcTrialSolution = Best1Exp;
			calcTrialSolution = &DESolver::Best1Exp;			
			break;

		case stRand1Exp:
			//calcTrialSolution = Rand1Exp;
			calcTrialSolution = &DESolver::Rand1Exp;
			break;

		case stRandToBest1Exp:
			//calcTrialSolution = RandToBest1Exp;
			calcTrialSolution = &DESolver::RandToBest1Exp;
			break;

		case stBest2Exp:
			//calcTrialSolution = Best2Exp;
			calcTrialSolution = &DESolver::Best2Exp;
			break;

		case stRand2Exp:
			//calcTrialSolution = Rand2Exp;
			calcTrialSolution = &DESolver::Rand2Exp;
			break;

		case stBest1Bin:
			//calcTrialSolution = Best1Bin;
			calcTrialSolution = &DESolver::Best1Bin;
			break;

		case stRand1Bin:
			//calcTrialSolution = Rand1Bin;
			calcTrialSolution = &DESolver::Rand1Bin;
			break;

		case stRandToBest1Bin:
			//calcTrialSolution = RandToBest1Bin;
			calcTrialSolution = &DESolver::RandToBest1Bin;
			break;

		case stBest2Bin:
			//calcTrialSolution = Best2Bin;
			calcTrialSolution = &DESolver::Best2Bin;
			break;

		case stRand2Bin:
			//calcTrialSolution = Rand2Bin;
			calcTrialSolution = &DESolver::Rand2Bin;
			break;
	}

	return;
}


void DESolver::SetupRangeConstraints(double *min,double *max,bool *inrange)
{
	bRangeConstrained = true;
	dMin = min;
	dMax = max;
	bInRange = inrange;	

	return;
}

void DESolver::SetupIntegerConstraints(bool *isint)
{	
	bIntegerConstrained = true;
	bIsInteger = isint;

	return;
}

void DESolver::SetupDiscreteConstraints(bool*isdiscrete,double**discvaltab,int*discvalnum)
{
	bDiscreteConstrained = true;
	bIsDiscrete = isdiscrete;	
	nDiscreteValNum = discvalnum;
	dDiscreteValTab = discvaltab;

	// If a variable is Discrete, it must also be Integer and InRange.
	// So, when setting up Discrete constraints, the InRange and Integer
	//     constraints should also be set up correctly.

	return;
}

void DESolver::ForceSolutionInRange(double *solution)
{
	for (int i=0; bRangeConstrained && (i < nDim); i++)
	{
		if (bInRange[i] && ((solution[i] < dMin[i]) || (solution[i] > dMax[i])))
		{
			double r = RandomUniform(0.0,1.0);
			solution[i] = r*(dMax[i]-dMin[i]) + dMin[i];
		}
	}

	return;
}

void DESolver::GetEvaluteSolution(double *evasolution,double *insolution)
{
	CopyVector(evasolution,insolution);
	
	for (int i=0; bIntegerConstrained && (i < nDim); i++)
	{
		if (bIsInteger[i])
			evasolution[i] = IntVar(insolution[i]);
	}
	
	for (int i=0; bDiscreteConstrained && (i < nDim); i++)
	{
		if (bIsDiscrete[i])
		{
			int id = IntVar(insolution[i]);
			evasolution[i] = dDiscreteValTab[i][id];
		}
	}

	return;
}

void DESolver::ShowDEParameters(void)
{
	printf("*** Differential Evolution (DE) Minimization ***\n");
	printf("Solution_Dim = %d\n", nDim);
	printf("Population_Num = %d\n", nPop);
	printf("Scale_Factor = %.4f\n", scale);
	printf("Crossover_Ratio = %.4f\n", probability);
	printf("DE_Strategy = ");
	switch (strategy)
	{
		case stBest1Exp:
			printf("Best1Exp");
			break;
		case stRand1Exp:
			printf("Rand1Exp");
			break;
		case stRandToBest1Exp:
			printf("RandToBest1Exp");
			break;
		case stBest2Exp:
			printf("Best2Exp");
			break;
		case stRand2Exp:
			printf("Rand2Exp");
			break;
		case stBest1Bin:
			printf("Best1Bin");
			break;
		case stRand1Bin:
			printf("Rand1Bin");
			break;
		case stRandToBest1Bin:
			printf("RandToBest1Bin");
			break;
		case stBest2Bin:
			printf("Best2Bin");
			break;
		case stRand2Bin:
			printf("Rand2Bin");
			break;
	}
	printf("\n************************************************\n");

	return;
}


double DESolver::Solve(int maxGenerations)
{
	int generation;
	int candidate;
	bool bAtSolution;
    int breakjudge=0;

	bestEnergy = 1.0E20;
	bAtSolution = false;
	double *evaluateSolution = new double[nDim];

	for (int generation = 0; generation < maxGenerations; generation++)
	{
        if(breakjudge==1)
        {
            break;
        }
		mexPrintf("Generation... %d\n", generation);
		for (candidate = 0; candidate < nPop; candidate++)
		{
			mexPrintf("%d -- %d 	",generation,candidate);
			(this->*calcTrialSolution)(candidate);
			ForceSolutionInRange(trialSolution);
			GetEvaluteSolution(evaluateSolution, trialSolution);
            /*if(candidate==0&&generation==0)
            {
                trialSolution[0]=5;
                trialSolution[1]=0.3;
                //trialSolution[2]=0.1;
                trialSolution[2]=4;
                trialSolution[3]=3;
                trialSolution[4]=21;
                trialSolution[5]=3;
                trialSolution[6]=0.5;
                //trialSolution[8]=300;

            }*/
            trialEnergy = EnergyFunction(trialSolution, bAtSolution);

			if (trialEnergy < popEnergy[candidate])
			{
				// New low for this candidate
				popEnergy[candidate] = trialEnergy;
				CopyVector(RowVector(population, candidate), trialSolution);

				// Check if all-time low
				if (trialEnergy < bestEnergy)
				{
					bestEnergy = trialEnergy;
					CopyVector(bestSolution, trialSolution);
				}
			}
            if(trialEnergy<=fusing_th)
            {
               breakjudge=1;
               break;
            }
		}
	}

    if (evaluateSolution)  delete evaluateSolution;
	generations = generation;
    GetEvaluteSolution(bestSolution,bestSolution);
    return(-bestEnergy);
}

void DESolver::Best1Exp(int candidate)
{
	int r1, r2;
	int n;

	SelectSamples(candidate,&r1,&r2);
	n = (int)RandomUniform(0.0,(double)nDim);

	CopyVector(trialSolution,RowVector(population,candidate));
	for (int i=0; (RandomUniform(0.0,1.0) < probability) && (i < nDim); i++) 
	{
		trialSolution[n] = bestSolution[n]
							+ scale * (Element(population,r1,n)
							- Element(population,r2,n));
		n = (n + 1) % nDim;
	}

	return;
}

void DESolver::Rand1Exp(int candidate)
{
	int r1, r2, r3;
	int n;

	SelectSamples(candidate,&r1,&r2,&r3);
	n = (int)RandomUniform(0.0,(double)nDim);

	CopyVector(trialSolution,RowVector(population,candidate));
	for (int i=0; (RandomUniform(0.0,1.0) < probability) && (i < nDim); i++) 
	{
		trialSolution[n] = Element(population,r1,n)
							+ scale * (Element(population,r2,n)
							- Element(population,r3,n));
		n = (n + 1) % nDim;
	}

	return;
}

void DESolver::RandToBest1Exp(int candidate)
{
	int r1, r2;
	int n;

	SelectSamples(candidate,&r1,&r2);
	n = (int)RandomUniform(0.0,(double)nDim);

	CopyVector(trialSolution,RowVector(population,candidate));
	for (int i=0; (RandomUniform(0.0,1.0) < probability) && (i < nDim); i++) 
	{
		trialSolution[n] += scale * (bestSolution[n] - trialSolution[n])
							 + scale * (Element(population,r1,n)
							 - Element(population,r2,n));
		n = (n + 1) % nDim;
	}

	return;
}

void DESolver::Best2Exp(int candidate)
{
	int r1, r2, r3, r4;
	int n;

	SelectSamples(candidate,&r1,&r2,&r3,&r4);
	n = (int)RandomUniform(0.0,(double)nDim);

	CopyVector(trialSolution,RowVector(population,candidate));
	for (int i=0; (RandomUniform(0.0,1.0) < probability) && (i < nDim); i++) 
	{
		trialSolution[n] = bestSolution[n] +
							scale * (Element(population,r1,n)
										+ Element(population,r2,n)
										- Element(population,r3,n)
										- Element(population,r4,n));
		n = (n + 1) % nDim;
	}

	return;
}

void DESolver::Rand2Exp(int candidate)
{
	int r1, r2, r3, r4, r5;
	int n;

	SelectSamples(candidate,&r1,&r2,&r3,&r4,&r5);
	n = (int)RandomUniform(0.0,(double)nDim);

	CopyVector(trialSolution,RowVector(population,candidate));
	for (int i=0; (RandomUniform(0.0,1.0) < probability) && (i < nDim); i++) 
	{
		trialSolution[n] = Element(population,r1,n)
							+ scale * (Element(population,r2,n)
										+ Element(population,r3,n)
										- Element(population,r4,n)
										- Element(population,r5,n));
		n = (n + 1) % nDim;
	}

	return;
}

void DESolver::Best1Bin(int candidate)
{
	int r1, r2;
	int n;

	SelectSamples(candidate,&r1,&r2);
	n = (int)RandomUniform(0.0,(double)nDim);

	CopyVector(trialSolution,RowVector(population,candidate));
	for (int i=0; i < nDim; i++) 
	{
		if ((RandomUniform(0.0,1.0) < probability) || (i == (nDim - 1)))
			trialSolution[n] = bestSolution[n]
								+ scale * (Element(population,r1,n)
											- Element(population,r2,n));
		n = (n + 1) % nDim;
	}

	return;
}

void DESolver::Rand1Bin(int candidate)
{
	int r1, r2, r3;
	int n;

	SelectSamples(candidate,&r1,&r2,&r3);
	n = (int)RandomUniform(0.0,(double)nDim);

	CopyVector(trialSolution,RowVector(population,candidate));
	for (int i=0; i < nDim; i++) 
	{
		if ((RandomUniform(0.0,1.0) < probability) || (i  == (nDim - 1)))
			trialSolution[n] = Element(population,r1,n)
								+ scale * (Element(population,r2,n)
												- Element(population,r3,n));
		n = (n + 1) % nDim;
	}

	return;
}

void DESolver::RandToBest1Bin(int candidate)
{
	int r1, r2;
	int n;

	SelectSamples(candidate,&r1,&r2);
	n = (int)RandomUniform(0.0,(double)nDim);

	CopyVector(trialSolution,RowVector(population,candidate));
	for (int i=0; i < nDim; i++) 
	{
		if ((RandomUniform(0.0,1.0) < probability) || (i  == (nDim - 1)))
			trialSolution[n] += scale * (bestSolution[n] - trialSolution[n])
									+ scale * (Element(population,r1,n)
												- Element(population,r2,n));
		n = (n + 1) % nDim;
	}

	return;
}

void DESolver::Best2Bin(int candidate)
{
	int r1, r2, r3, r4;
	int n;

	SelectSamples(candidate,&r1,&r2,&r3,&r4);
	n = (int)RandomUniform(0.0,(double)nDim);

	CopyVector(trialSolution,RowVector(population,candidate));
	for (int i=0; i < nDim; i++) 
	{
		if ((RandomUniform(0.0,1.0) < probability) || (i  == (nDim - 1)))
			trialSolution[n] = bestSolution[n]
								+ scale * (Element(population,r1,n)
											+ Element(population,r2,n)
											- Element(population,r3,n)
											- Element(population,r4,n));
		n = (n + 1) % nDim;
	}

	return;
}

void DESolver::Rand2Bin(int candidate)
{
	int r1, r2, r3, r4, r5;
	int n;

	SelectSamples(candidate,&r1,&r2,&r3,&r4,&r5);
	n = (int)RandomUniform(0.0,(double)nDim);

	CopyVector(trialSolution,RowVector(population,candidate));
	for (int i=0; i < nDim; i++) 
	{
		if ((RandomUniform(0.0,1.0) < probability) || (i  == (nDim - 1)))
			trialSolution[n] = Element(population,r1,n)
								+ scale * (Element(population,r2,n)
											+ Element(population,r3,n)
											- Element(population,r4,n)
											- Element(population,r5,n));
		n = (n + 1) % nDim;
	}

	return;
}

void DESolver::SelectSamples(int candidate,int *r1,int *r2,
										int *r3,int *r4,int *r5)
{
	if (r1)
	{
		do
		{
			*r1 = (int)RandomUniform(0.0,(double)nPop);
		}
		while (*r1 == candidate);
	}

	if (r2)
	{
		do
		{
			*r2 = (int)RandomUniform(0.0,(double)nPop);
		}
		while ((*r2 == candidate) || (*r2 == *r1));
	}

	if (r3)
	{
		do
		{
			*r3 = (int)RandomUniform(0.0,(double)nPop);
		}
		while ((*r3 == candidate) || (*r3 == *r2) || (*r3 == *r1));
	}

	if (r4)
	{
		do
		{
			*r4 = (int)RandomUniform(0.0,(double)nPop);
		}
		while ((*r4 == candidate) || (*r4 == *r3) || (*r4 == *r2) || (*r4 == *r1));
	}

	if (r5)
	{
		do
		{
			*r5 = (int)RandomUniform(0.0,(double)nPop);
		}
		while ((*r5 == candidate) || (*r5 == *r4) || (*r5 == *r3)
													|| (*r5 == *r2) || (*r5 == *r1));
	}

	return;
}

/*------Constants for RandomUniform()---------------------------------------*/
#define SEED 3
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double DESolver::RandomUniform(double minValue,double maxValue)
{
	long j;
	long k;
	static long idum;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	double result;

	if (iy == 0)
		idum = SEED;

	if (idum <= 0)
	{
		if (-idum < 1)
			idum = 1;
		else
			idum = -idum;

		idum2 = idum;

		for (j=NTAB+7; j>=0; j--)
		{
			k = idum / IQ1;
			idum = IA1 * (idum - k*IQ1) - k*IR1;
			if (idum < 0) idum += IM1;
			if (j < NTAB) iv[j] = idum;
		}

		iy = iv[0];
	}

	k = idum / IQ1;
	idum = IA1 * (idum - k*IQ1) - k*IR1;

	if (idum < 0)
		idum += IM1;

	k = idum2 / IQ2;
	idum2 = IA2 * (idum2 - k*IQ2) - k*IR2;

	if (idum2 < 0)
		idum2 += IM2;

	j = iy / NDIV;
	iy = iv[j] - idum2;
	iv[j] = idum;

	if (iy < 1)
		iy += IMM1;

	result = AM * iy;

	if (result > RNMX)
		result = RNMX;

	result = minValue + result * (maxValue - minValue);
	return(result);
}
