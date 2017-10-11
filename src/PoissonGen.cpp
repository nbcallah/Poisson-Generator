#include <cmath>
#include <vector>
#include "../inc/PoissonGen.hpp"

poissonGen::poissonGen(double maxRate) {
	this->maxRate = maxRate;
	this->currentTime = 0.0;
}

poissonGen::~poissonGen() {
}

std::vector<double> poissonGen::getNEvents(unsigned int num) {
	std::vector<double> events;
	for(unsigned int i = 0; i < num; i++) {
		events.push_back(this->getNextEvent());
	}
	return events;
}

double poissonGen::getNextEvent() {
	currentTime += -log(this->getUniformNumber())/(it->maxRate);
	while(this->getUniformNumber() > this->evaluateFunction(this->currentTime)/this->maxRate) {
		currentTime += -log(this->getUniformNumber())/(it->maxRate);
	}
	
	return currentTime;
}

std::vector<double> poissonGen::getEventsToTime(double length) {
	std::vector<double> events;
	event = this->getNextEvent();
	while(event < length) {
		events.push_back(event);
		event = this->getNextEvent();
	}

	return events;
}

void poissonGen::setTime(double time) {
	this->currentTime = time;
}
