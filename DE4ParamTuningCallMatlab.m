%% author: Zezheng Wang, Shijie Zhang, Wei Feng
% A demo about "Differential Evolution Algorithm" with matlab and c++
% aim to find the minimal value of LossFunction()
% ParamMatrix(i,:)=[parai_min, parai_max, isInt];
% In this demo: loss=(x-1)^2 + (y-2)^2 + (z-3)^2 + (q-4)^2;

clear, clc, close all;

global maxscore;
global runnumber;
global besterr;
besterr=10000;
maxscore=0;
runnumber=0;

nDim = 4;
nPop = 50*nDim;
fusing_th=-11;

ParamMatrix = zeros(nDim, 3);

para0_min=-1.0;
para0_max=10.0; 
ParamMatrix(1,:) = [para0_min, para0_max, 0];

para1_min=-1.0;
para1_max=10.0;
ParamMatrix(2,:) = [para1_min, para1_max, 0];

para2_min=-1.0;
para2_max=10.0;
ParamMatrix(3,:) = [para2_min, para2_max, 0];

para3_min=-1.0;
para3_max=10.0;
ParamMatrix(4,:) = [para3_min, para3_max, 0];

mexDE4ParamTuningCallMatlab( nDim, nPop, fusing_th, ParamMatrix );
