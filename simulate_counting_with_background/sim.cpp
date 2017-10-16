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

class dipStructure {
	public:
	dipStructure(int numDips, double* dipTCs, double* dipLengths, double* dipStrengths, double numHits, double hold_t, std::mt19937_64* gen);
	double getShortMean();
	double getLongMean();
	
	private:
	std::vector<double> dipLengths;
	std::vector<double> offsets;
	std::vector<exponentialCountingProcess> dipShort;
	std::vector<exponentialCountingProcess> dipLong;
	constantCountingProcess bkg;
};

double dipStructure::getShortMean() {
	double sumShort = 0.0;
	double ctsShort = 0.0;
	for(int i = 0; i < dipLengths.size(); i++) {
		dipShort[i].setTime(0.0);
		std::vector<double> cts = dipShort[i].getEventsToTime(dipLengths[i]);
		sumShort += std::accumulate(cts.begin(), cts.end(), 0.0);
		sumShort += offsets[i]*cts.size();
		ctsShort += cts.size();
		bkg.setTime(0.0);
		std::vector<double> bkgCts = bkg.getEventsToTime(dipLengths[i]);
		sumShort += std::accumulate(bkgCts.begin(), bkgCts.end(), 0.0);
		sumShort += offsets[i]*bkgCts.size();
		ctsShort += bkgCts.size();
	}
	return sumShort/ctsShort;
}

double dipStructure::getLongMean() {
	double sumLong = 0.0;
	double ctsLong = 0.0;
	for(int i = 0; i < dipLengths.size(); i++) {
		dipLong[i].setTime(0.0);
		std::vector<double> cts = dipLong[i].getEventsToTime(dipLengths[i]);
		sumLong += std::accumulate(cts.begin(), cts.end(), 0.0);
		sumLong += offsets[i]*cts.size();
		ctsLong += cts.size();
		bkg.setTime(0.0);
		std::vector<double> bkgCts = bkg.getEventsToTime(dipLengths[i]);
		sumLong += std::accumulate(bkgCts.begin(), bkgCts.end(), 0.0);
		sumLong += offsets[i]*bkgCts.size();
		ctsLong += bkgCts.size();
	}
	return sumLong/ctsLong;
}

dipStructure::dipStructure(int numDips, double* dipTCs, double* dipLengths, double* dipStrengths, double numHits, double hold_t, std::mt19937_64* gen) : bkg(0.15, gen) {
	double dip_norm = 0.0;
	for(int i = 0; i < numDips; i++) {
		dip_norm += dipStrengths[i]*(1-exp(-(dipLengths[i])/dipTCs[i]));
	}
	dip_norm = 1.0/dip_norm;
	
	for(int i = 0; i < numDips; i++) {
		offsets.push_back(std::accumulate(this->dipLengths.begin(), this->dipLengths.end(), 0.0));
		this->dipLengths.push_back(dipLengths[i]);
		dipShort.push_back(exponentialCountingProcess(numHits*dip_norm*dipStrengths[i], dipTCs[i], gen));
		dipLong.push_back(exponentialCountingProcess(numHits*exp(-hold_t/877.7)*dip_norm*dipStrengths[i], dipTCs[i], gen));
	}
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
	
//	double dip9_TCs[9] = {100, 65.981913, 23.739640, 19.241204, 18.700906, 17.779887, 19.758315, 14.361219, 8.065494};
//	double dip9_lengths[9] = {20, 40, 20, 20, 20, 20, 20, 20, 60};
//	double dip9_strengths[9] = {0.00000000001, 0.002746, 0.010865, 0.013028, 0.011499, 0.012766, 0.009688, 0.008515, 0.003804};
//	dipStructure dips(9, &dip9_TCs[0], &dip9_lengths[0], &dip9_strengths[0], 21700, 990, &gen);
	
	double dip1_TCs[1] = {8.0};
	double dip1_lengths[1] = {100.0};
	double dip1_strengths[1] = {1.0};
	dipStructure dips(1, &dip1_TCs[0], &dip1_lengths[0], &dip1_strengths[0], 22000, 990, &gen);
	
	std::vector<double> avgDifferences;
	const int numSamples = 200;
	const int numRuns = 55;
	for(int i = 0; i < numSamples; i++) {
		std::vector<double> meanTimeShort;
		std::vector<double> meanTimeLong;
		for(int j = 0; j < numRuns; j++) {
			meanTimeShort.push_back(dips.getShortMean());
		}
		for(int j = 0; j < numRuns; j++) {
			meanTimeLong.push_back(dips.getLongMean());
		}
		double avgDiff = 0.0;
		for(int j = 0; j < numRuns; j++) {
			avgDiff += (meanTimeLong[j]-meanTimeShort[j]);
		}
		avgDiff = avgDiff/numRuns;
		avgDifferences.push_back(avgDiff);
		double stddev = 0.0;
		for(int j = 0; j < numRuns; j++) {
			stddev += pow((meanTimeLong[j]-meanTimeShort[j])-avgDiff, 2);
		}
		stddev = stddev/(numRuns - 1);
		stddev = sqrt(stddev);
		printf("Avg: %f stddev: %f\n", avgDiff, stddev/sqrt(1));
	}
	
	double mean = std::accumulate(avgDifferences.begin(), avgDifferences.end(), 0.0) / avgDifferences.size();
	double meanstddev = sqrt(std::accumulate(avgDifferences.begin(), avgDifferences.end(), 0.0, [mean](double sum, double val){return sum + pow(val - mean, 2);})/(avgDifferences.size()-1));
	printf("Mean: %f stddev: %f\n", mean, meanstddev);
	
	
	return 0;
}