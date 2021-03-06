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
	std::mt19937_64* gen;
    std::uniform_real_distribution<double> randGen;
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

//Holds evaluation of a rate that is constant
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

//Bounding rate is just the rate (all generated events are accepted and there is no thinning)
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

//A class that holds one or more exponential processes plus a background process.
//Models a counting experiment where a detector is sensitive to different populations
//At discrete times. Each population drains for a while before moving to the next.
//Counting process is modeled to have population decay, and the population can be counted
//At varying times.
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

//Get mean arrival time for when population is counted during the beginning
double dipStructure::getShortMean() {
	double sumShort = 0.0;
	double ctsShort = 0.0;
	//For each population that is counted
	for(int i = 0; i < dipLengths.size(); i++) {
		//Set time of poisson process to 0
		dipShort[i].setTime(0.0);
		//Generate events out to a the length of the counting window
		std::vector<double> cts = dipShort[i].getEventsToTime(dipLengths[i]);
		//Sum the arrival times
		sumShort += std::accumulate(cts.begin(), cts.end(), 0.0);
		sumShort += offsets[i]*cts.size();
		ctsShort += cts.size();
		//Generate a constant background process
		bkg.setTime(0.0);
		std::vector<double> bkgCts = bkg.getEventsToTime(dipLengths[i]);
		//Add in arrival times
		sumShort += std::accumulate(bkgCts.begin(), bkgCts.end(), 0.0);
		sumShort += offsets[i]*bkgCts.size();
		ctsShort += bkgCts.size();
	}
	//return mean time
	return sumShort/ctsShort;
}

//Get mean arrival time for when population is counted after storage and decay for a length of time
double dipStructure::getLongMean() {
	double sumLong = 0.0;
	double ctsLong = 0.0;
	//Do everything the same, but with alternate timings
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

//Constructor for counting experiment simulation.
dipStructure::dipStructure(int numDips, double* dipTCs, double* dipLengths, double* dipStrengths, double numHits, double hold_t, std::mt19937_64* gen) : bkg(0.15, gen) {
	//Normalize the given structure so that the parameters can be multiplied by the desired population size.
	double dip_norm = 0.0;
	for(int i = 0; i < numDips; i++) {
		dip_norm += dipStrengths[i]*(1-exp(-(dipLengths[i])/dipTCs[i]));
	}
	dip_norm = 1.0/dip_norm;
	
	for(int i = 0; i < numDips; i++) {
		//Construct a list of counting windows for the experiment
		//First, find the offset for the ith counting window based on previous entries
		offsets.push_back(std::accumulate(this->dipLengths.begin(), this->dipLengths.end(), 0.0));
		//Add the end of the current counting process.
		this->dipLengths.push_back(dipLengths[i]);
		//Create 2 processes: one for holding the decaying population for no time, and the other for a given time
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
	
	//These arrays hold parameters so that different counting experiments can be tested.
//	double dip9_TCs[9] = {100, 65.981913, 23.739640, 19.241204, 18.700906, 17.779887, 19.758315, 14.361219, 8.065494};
//	double dip9_lengths[9] = {20, 40, 20, 20, 20, 20, 20, 20, 60};
//	double dip9_strengths[9] = {0.0000000000001, 0.002746, 0.010865, 0.013028, 0.011499, 0.012766, 0.009688, 0.008515, 0.003804};
//	dipStructure dips(9, &dip9_TCs[0], &dip9_lengths[0], &dip9_strengths[0], 21700, 990, &gen);
	
	double dip1_TCs[1] = {8.0};
	double dip1_lengths[1] = {100.0};
	double dip1_strengths[1] = {1.0};
	dipStructure dips(1, &dip1_TCs[0], &dip1_lengths[0], &dip1_strengths[0], 22000, 990, &gen);
	
//	double dip3_TCs[3] = {100, 52, 7.8};
//	double dip3_lengths[3] = {20, 20, 100};
//	double dip3_strengths[3] = {0.0000000000001, 0.0027, 0.083};
//	dipStructure dips(3, &dip3_TCs[0], &dip3_lengths[0], &dip3_strengths[0], 28500, 1370, &gen);
	
	//We want to see how much the mean arrival time changes when the signal to noise of the process changes.
	//We'll simulate lots of poisson processes and measure the mean arrival time of the process under high
	//and low signal to noise.
	std::vector<double> avgDifferences;
	//Number of averages to compute
	const int numSamples = 5;
	//Number of simulation runs per average
	const int numRuns = 200;
	for(int i = 0; i < numSamples; i++) {
		std::vector<double> meanTimeShort;
		std::vector<double> meanTimeLong;
		//Measure the mean arrival time for numRuns experiments (both high and low signal)
		for(int j = 0; j < numRuns; j++) {
			meanTimeShort.push_back(dips.getShortMean());
		}
		for(int j = 0; j < numRuns; j++) {
			meanTimeLong.push_back(dips.getLongMean());
		}
		//Compute average differences in mean arrival time
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
		printf("Avg Arrival Time Shift: %f stddev: %f\n", avgDiff, stddev/sqrt(1));
	}
	
	//Calculate grand average
	double mean = std::accumulate(avgDifferences.begin(), avgDifferences.end(), 0.0) / avgDifferences.size();
	double meanstddev = sqrt(std::accumulate(avgDifferences.begin(), avgDifferences.end(), 0.0, [mean](double sum, double val){return sum + pow(val - mean, 2);})/(avgDifferences.size()-1));
	printf("Grand Avg Arrival Time Shift: %f stddev: %f\n", mean, meanstddev);
	
	
	return 0;
}