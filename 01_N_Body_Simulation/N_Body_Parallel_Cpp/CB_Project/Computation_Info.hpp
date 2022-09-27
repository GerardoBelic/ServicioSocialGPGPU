#ifndef COMPUTATION_INFO
#define COMPUTATION_INFO

struct Computation_Info
{

	Computation_Info(float _timeStep, unsigned _numIterations, unsigned _numThreads = 1) :
					 timeStep(_timeStep), numIterations(_numIterations), numThreads(_numThreads) { }

	float timeStep;
	unsigned numIterations;
	unsigned numThreads;
};

#endif
