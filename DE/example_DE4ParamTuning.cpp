// Parameters training by DE-MaxF
class DEMaxFSolver : public DESolver
{
public:
	DEMaxFSolver(int dim, int pop) : DESolver(dim,pop), Counter(0) {;}

	ofstream ofResfile;
	string   strDEBestParamFn;

	void Prepare(const string &strFileList, 
				 map<string, string> &mapSyllDict, 
				 const struct Parameter &stPara) 
	{
		CopyParam(inParameter, stPara);
		(this->*calcTrialSolution)(0);
		ForceSolutionInRange(trialSolution);
		GetEvaluteSolution(bestSolution,trialSolution);
		Solution2Param(inParameter, bestSolution);
		Process_Prepare(vecAllDataSet,vecAllmapDict,vecAllSentences,vecAllBndId,
					    strFileList, mapSyllDict, inParameter);
		CopyParam(inParameter, stPara);
	};
	double EnergyFunction(double *pTrialSolution, bool &bAtSolution);

	void Param2Solution(double *pSolution, const struct Parameter &stPara);
	void Solution2Param(struct Parameter &stPara, double *pSolution);

	void Solution2File(const string strParamFn);
	void Solution2File() { if (!strDEBestParamFn.empty()) Solution2File(strDEBestParamFn); };
private:
	struct Parameter inParameter;
	vector< vector<string> > vecAllDataSet; 
	vector< map<string, int> > vecAllmapDict;
	vector< vector<Sentence> > vecAllSentences;
	vector< vector<int> > vecAllBndId;
	int Counter;
};

void DEMaxFSolver::Param2Solution(double *pSolution, const struct Parameter &stPara)
{
	pSolution[0] = (double)stPara.iBlockLen;
	pSolution[1] = (double)stPara.iCutOff;
	pSolution[2] = (double)stPara.iSmoothLen;
	pSolution[3] = stPara.dSmoothAlpha;
	pSolution[4] = stPara.dBeta;
}

void DEMaxFSolver::Solution2Param(struct Parameter &stPara, double *pSolution)
{
	CopyParam(stPara, inParameter);
	stPara.iBlockLen = (int)pSolution[0];
	stPara.iCutOff = (int)pSolution[1];
	stPara.iSmoothLen = (int)pSolution[2];
	stPara.dSmoothAlpha = pSolution[3];
	stPara.dBeta = pSolution[4];
}

double DEMaxFSolver::EnergyFunction(double *pTrialSolution, bool &bAtSolution)
{
	struct Parameter curPara;
	Solution2Param(curPara, pTrialSolution);
	double curEngy = -Process_silent_lite(curPara, vecAllDataSet, vecAllmapDict, 
		                                  vecAllSentences, vecAllBndId);
	if (Counter++ % nPop == 0)
	{
		double dFm = (Counter==1) ? -curEngy : -Energy();
		printf("#Iter %d,  F-measure: %.6f\n", Counter/nPop + 1, dFm);
#ifndef SILENCE_MODE_DE
		ofResfile << dFm << endl;
#endif
		if (Counter/nPop > 0)
			Solution2File();
	}
	return curEngy;
}

void DEMaxFSolver::Solution2File(const string strParamFn)
{
	struct Parameter bestPara;
	CopyParam(bestPara, inParameter);
	Solution2Param(bestPara, Solution());
	ParamToFile(strParamFn, bestPara);
}

void Train_by_DE_F(const struct Parameter &stPara, 
				   map<string, string> &mapSyllDict,
				   const string &strFileList, const string &strMethodNam)

{
	// initialize
	int nD   = 5;
	DEMaxFSolver deSolver(nD, 20);
	
	// setup DE
	double *pdMin = new double[nD];
	double *pdMax = new double[nD];
	pdMin[0] = (double)stPara.iBlockLen_mn;
	pdMax[0] = (double)stPara.iBlockLen_mx;
	pdMin[1] = (double)stPara.iCutOff_mn;
	pdMax[1] = (double)stPara.iCutOff_mx;
	pdMin[2] = (double)stPara.iSmoothLen_mn;
	pdMax[2] = (double)stPara.iSmoothLen_mx;
	pdMin[3] = stPara.dSmoothAlpha_mn;
	pdMax[3] = stPara.dSmoothAlpha_mx;
	pdMin[4] = stPara.dBeta_mn;
	pdMax[4] = stPara.dBeta_mx;
	deSolver.Setup(pdMin, pdMax, g_nDEStrategy, g_dDEScaleF, g_dDECrossRatio);

	// setup DE constraints (no discrete constraints)
	bool *pbRange = new bool[nD];
	bool *pbInt = new bool[nD];
	for (int i=0; i<nD; i++)
	{
		pbRange[i] = true;
		pbInt[i] = true;
	}
	pbInt[3] = false;
	pbInt[4] = false;
	deSolver.SetupRangeConstraints(pdMin, pdMax, pbRange);
	deSolver.SetupIntegerConstraints(pbInt);

	// show DE param
	deSolver.ShowDEParameters();

	// run DE
#ifndef SILENCE_MODE_DE
	string saResFn = strMethodNam;
	saResFn.append("_DE-MaxF_Res.txt");
	deSolver.ofResfile.open(saResFn.data(), ios_base::out);
#endif
	deSolver.strDEBestParamFn = strMethodNam;
	deSolver.strDEBestParamFn.append("_DE-MaxF_Best.ini");
	deSolver.Prepare(strFileList, mapSyllDict, stPara);
	deSolver.Solve(g_nDEMaxGeneration);
#ifndef SILENCE_MODE_DE
	deSolver.ofResfile.close();
#endif

	// output solution
	deSolver.Solution2File();

#ifndef SILENCE_MODE
	cout<<"--------------------------------------------------------"<<endl;
#endif
	cout << "DEBest_Fmeasure: " << -deSolver.Energy() << ",  DEBest_ParamFileName: " << deSolver.strDEBestParamFn << endl;
	
	delete [] pdMin;
	delete [] pdMax;
	delete [] pbRange;
	delete [] pbInt;
}