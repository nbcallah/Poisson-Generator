#include "../inc/PoissonGen.hpp"
#include <cmath>
#include <stdio.h>
#include <random>

//Holds evaluation of a rate that goes ~exp(tau)
class exponentialCountingProcess : public poissonGen {
	public:
	exponentialCountingProcess(double num, double tau, std::mt19937_64* gen);
	
	private:
	virtual double getUniformNumber();
	virtual double evaluateFunction(double t);
	//Stores RNG
	std::mt19937_64* gen;
    std::uniform_real_distribution<double> randGen;
	//Parameters for function
	double tau;
	double num;
};

//Bounding rate on process is at t=0 with rate=(num/tau)*exp(0)=num/tau
exponentialCountingProcess::exponentialCountingProcess(double num, double tau, std::mt19937_64* gen) : poissonGen(num/tau) {
	this->gen = gen;
	this->randGen = std::uniform_real_distribution<double>(0.0, 1.0);
	this->tau = tau;
	this->num = num;
}

double exponentialCountingProcess::getUniformNumber() {
	return this->randGen(*(this->gen));
}

//Evaluate exponential rate
double exponentialCountingProcess::evaluateFunction(double t) {
	return (this->num/this->tau)*exp(-t/this->tau);
}

int main(int argc, char** argv) {
	//Create seed
	std::vector<uint64_t> seed;
	seed.push_back(2736687128);
	seed.push_back(234302120);
	seed.push_back(3355772407);
	seed.push_back(657836083);
	std::seed_seq sseq(seed.begin(), seed.end());
	std::mt19937_64 gen;
	gen.seed(sseq);
	
	//Create counting process object
	exponentialCountingProcess p(100.0, 10.0, &gen);
	
	//Get events out to 20s
	std::vector<double> realization = p.getEventsToTime(20);
	
	//Print number of events
	printf("Generated %lu Events!\n", realization.size());
	
	return 0;
}