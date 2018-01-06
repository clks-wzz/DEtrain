# DEtrain
Matlab and C++ implementation of "Differential Evolution Algorithm"

# Run
1. mex -O mexDE4ParamTuningCallMatlab.cpp
2. run DE4ParamTuningCallMatlab.m

# Make your own trainer
1. encapsulate your loss in the "LossFunction.m"
2. adapt your own params in the "DE4ParamTuningCallMatlab.m"
