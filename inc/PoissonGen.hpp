#include <vector>

#ifndef POISSONGEN_H
#define POISSONGEN_H

class poissonGen {
	public:
	poissonGen(double maxRate);
	~poissonGen();
	double getNextEvent(double maxTime);
	std::vector<double> getEventsToTime(double length);
	void setTime(double time);
	
	private:
	virtual double getUniformNumber() = 0;
	virtual double evaluateFunction(double t) = 0;
	
	double maxRate;
	double currentTime;
};

#endif