#include "../inc/PoissonGen.hpp"
#include <cmath>
#include <stdio.h>
#include <random>

class exponentialCountingProcess : public poissonGen {
	public:
	exponentialCountingProcess(double num, double tau, std::mt19937_64* gen);
	
	private:
	virtual double getUniformNumber();
	virtual double evaluateFunction(double t);
	std::mt19937_64* gen;
    std::uniform_real_distribution<double> randGen;
	double tau;
	double num;
};

exponentialCountingProcess::exponentialCountingProcess(double num, double tau, std::mt19937_64* gen) : poissonGen(num/tau) {
	this->gen = gen;
	this->randGen = std::uniform_real_distribution<double>(0.0, 1.0);
	this->tau = tau;
	this->num = num;
}

double exponentialCountingProcess::getUniformNumber() {
	return this->randGen(*(this->gen));
}

double exponentialCountingProcess::evaluateFunction(double t) {
	return (this->num/this->tau)*exp(-t/this->tau);
}

class constantCountingProcess : public poissonGen {
	public:
	constantCountingProcess(double rate, std::mt19937_64* gen);
	
	private:
	virtual double getUniformNumber();
	virtual double evaluateFunction(double t);
	std::mt19937_64* gen;
    std::uniform_real_distribution<double> randGen;
	double rate;
};

constantCountingProcess::constantCountingProcess(double rate, std::mt19937_64* gen) : poissonGen(rate) {
	this->gen = gen;
	this->randGen = std::uniform_real_distribution<double>(0.0, 1.0);
	this->rate = rate;
}

double constantCountingProcess::getUniformNumber() {
	return this->randGen(*(this->gen));
}

double constantCountingProcess::evaluateFunction(double t) {
	return rate;
}

int main(int argc, char** argv) {
	std::vector<uint64_t> seed;
	seed.push_back(2736687128);
	seed.push_back(234302120);
	seed.push_back(3355772407);
	seed.push_back(657836083);
	std::seed_seq sseq(seed.begin(), seed.end());
	std::mt19937_64 gen;
	gen.seed(sseq);
	
	auto getMeanTime = [](std::vector<double> &sig, std::vector<double> &bkg)->double{
		double mean = 0.0;
		for(auto it = sig.begin(); it < sig.end(); it++) {
			mean += *it;
		}
		for(auto it = bkg.begin(); it < bkg.end(); it++) {
			mean += *it;
		}
		mean = mean / ((double)sig.size() + (double)bkg.size());
		return mean;
	};
	
	exponentialCountingProcess processLarge(20000, 7.0, &gen);
	exponentialCountingProcess processSmall(20000*exp(-1300.0/877.7), 7.0, &gen);
	constantCountingProcess processConst(0.1, &gen);
	
//	for(int j = 0; j < 100; j++) {
		std::vector<double> meanTimeLarge;
		std::vector<double> meanTimeSmall;
		for(int i = 0; i < 1000; i++) {
			processLarge.setTime(0.0);
			processConst.setTime(0.0);
			std::vector<double> largeCounts = processLarge.getEventsToTime(100.0);
			std::vector<double> constCountsLarge = processConst.getEventsToTime(100.0);
			meanTimeLarge.push_back(getMeanTime(largeCounts, constCountsLarge));
			processSmall.setTime(0.0);
			processConst.setTime(0.0);
			std::vector<double> smallCounts = processSmall.getEventsToTime(100.0);
			std::vector<double> constCountsSmall = processConst.getEventsToTime(100.0);
			meanTimeSmall.push_back(getMeanTime(smallCounts, constCountsSmall));
		}
		double avgDiff = 0.0;
		for(int i = 0; i < 1000; i++) {
			avgDiff += (meanTimeSmall[i]-meanTimeLarge[i]);
		}
		printf("Avg: %f\n", avgDiff/1000.0);
//	}
}