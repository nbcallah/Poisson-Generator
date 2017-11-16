#include <vector>

#ifndef POISSONGEN_H
#define POISSONGEN_H

class poissonGen {
	public:
	//Constructor
	poissonGen(double maxRate);
	~poissonGen();
	//Generate next realization of the poisson process
	double getNextEvent(double maxTime);
	//Generate events until the given time has elapsed
	std::vector<double> getEventsToTime(double length);
	//Set the time of the generator
	void setTime(double time);
	
	private:
	//Need to implement method to get a uniform number
	virtual double getUniformNumber() = 0;
	//Need to implement method to evaluate rate function
	virtual double evaluateFunction(double t) = 0;
	
	//Holds the expected max rate (for rejection sampling)
	double maxRate;
	//Holds the current time to evaluate rate function at
	double currentTime;
};

#endif