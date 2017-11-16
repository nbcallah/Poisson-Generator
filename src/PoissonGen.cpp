#include <cmath>
#include <vector>
#include "../inc/PoissonGen.hpp"

/* Constructor
	Arguments: maxRate - bounding rate of rate function.
*/
poissonGen::poissonGen(double maxRate) {
	this->maxRate = maxRate;
	this->currentTime = 0.0;
}

poissonGen::~poissonGen() {
}

/* Realize one event
	Arguments: maxTime - time when generator gives up and returns infinity
*/
double poissonGen::getNextEvent(double maxTime) {
	//Advance time by the bounding rate
	this->currentTime += -log(this->getUniformNumber())/(this->maxRate);
	//Reject samples based on R(t)/R_max
	while(this->getUniformNumber() > this->evaluateFunction(this->currentTime)/this->maxRate && this->currentTime < maxTime) {
		this->currentTime += -log(this->getUniformNumber())/(this->maxRate);
	}
	//When a sample is accepted, return the time
	return this->currentTime < maxTime ? this->currentTime : std::numeric_limits<double>::infinity();
}

/* Realize many events
	Arguments: length - time to generate events to
*/
std::vector<double> poissonGen::getEventsToTime(double length) {
	std::vector<double> events;
	//Generate a single event
	double event = this->getNextEvent(length);
	//While the arrival time is less than the length, append to events
	while(event < length) {
		events.push_back(event);
		event = this->getNextEvent(length);
	}

	return events;
}

/* Asynchronously set time
	Arguments: time - time
*/
void poissonGen::setTime(double time) {
	this->currentTime = time;
}
