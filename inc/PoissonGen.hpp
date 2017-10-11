#ifndef POISSONGEN_H
#define POISSONGEN_H

class poissonGen {
	public:
	poissonGen(double maxRate);
	~poissonGen();
	std::vector<double> getNEvents(unsigned int num);
	double getNextEvent();
	std::vector<double> getEventsToTime(double length);
	void setTime(double time);
	
	private:
	virtual double getUniformNumber();
	virtual double evaluateFunction(double t);
	
	double maxRate;
	double currentTime;
}

#endif