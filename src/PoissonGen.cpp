#include <cmath>
#include <vector>
#include "../inc/PoissonGen.hpp"

poissonGen::poissonGen(double maxRate) {
	this->maxRate = maxRate;
	this->currentTime = 0.0;
}

poissonGen::~poissonGen() {
}

double poissonGen::getNextEvent(double maxTime) {
	this->currentTime += -log(this->getUniformNumber())/(this->maxRate);
	while(this->getUniformNumber() > this->evaluateFunction(this->currentTime)/this->maxRate && this->currentTime < maxTime) {
		this->currentTime += -log(this->getUniformNumber())/(this->maxRate);
	}
	return this->currentTime < maxTime ? this->currentTime : std::numeric_limits<double>::infinity();
}

std::vector<double> poissonGen::getEventsToTime(double length) {
	std::vector<double> events;
	double event = this->getNextEvent(length);
	while(event < length) {
		events.push_back(event);
		event = this->getNextEvent(length);
	}

	return events;
}

void poissonGen::setTime(double time) {
	this->currentTime = time;
}
