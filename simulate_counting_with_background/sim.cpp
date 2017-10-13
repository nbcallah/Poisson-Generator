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
	
//	auto getMeanTime = [](std::vector<double> &sig, std::vector<double> &bkg)->double{
//		double mean = 0.0;
//		for(auto it = sig.begin(); it < sig.end(); it++) {
//			mean += *it;
//		}
//		for(auto it = bkg.begin(); it < bkg.end(); it++) {
//			mean += *it;
//		}
//		mean = mean / ((double)sig.size() + (double)bkg.size());
//		return mean;
//	};
	auto getSumTime = [](std::vector<double> &sig)->double{
		double sum = 0.0;
		for(auto it = sig.begin(); it < sig.end(); it++) {
			sum += *it;
		}
		return sum;
	};
	
//	dip dip1_arr[1] = {{7.8, 0.0, 1.0}};
//	const int numDips = 1;
//	dip dip3_arr[2] = {{52, 0.0, 0.0027}, {7.8, 20.0, 0.083}};
//  const int numDips = 2;
	double dip9_arr[8][3] = {{65.981913, 0.000000, 0.002746}, {23.739640, 40.000000, 0.010865}, {19.241204, 60.000000, 0.013028}, {18.700906, 80.000000, 0.011499}, {17.779887, 100.000000, 0.012766}, {19.758315, 120.000000, 0.009688}, {14.361219, 140.000000, 0.008515}, {8.065494, 160.000000, 0.003804}};
	const int numDips = 8;
	
	double dip_norm = 0.0;
		for(int i = 0; i < 7; i++) {
		dip_norm += dip9_arr[i][2]*(1-exp(-(dip9_arr[i+1][1]-dip9_arr[i][1])/dip9_arr[i][0]));
	}
	dip_norm += dip9_arr[7][2]*(1-exp(-100/dip9_arr[7][0]));
	dip_norm = 1.0/dip_norm;
	
	constantCountingProcess bkg(0.1176, &gen);
	std::vector<exponentialCountingProcess> dip9Short;
	std::vector<exponentialCountingProcess> dip9Long;
	for(int i = 0; i < 8; i++) {
		dip9Short.push_back(exponentialCountingProcess(15000*dip_norm*dip9_arr[i][2], dip9_arr[i][0], &gen));
		dip9Long.push_back(exponentialCountingProcess(15000*exp(-1020.0/877.7)*dip_norm*dip9_arr[i][2], dip9_arr[i][0], &gen));
	}
	const int numSamples = 100;
	const int numRuns = 70;
	for(int i = 0; i < numSamples; i++) {
		std::vector<double> meanTimeShort;
		std::vector<double> meanTimeLong;
		for(int j = 0; j < numRuns; j++) {
			double sumShort = 0.0;
			double ctsShort = 0;
			for(auto it = dip9Short.begin(); it < dip9Short.end(); it++) {
				it->setTime(0.0);
			}
			for(int dipNum = 0; dipNum < 7; dipNum++) {
				std::vector<double> cts = dip9Short[dipNum].getEventsToTime(dip9_arr[dipNum+1][1]-dip9_arr[dipNum][1]);
				sumShort += getSumTime(cts);
				sumShort += (dip9_arr[dipNum][1])*cts.size();
				ctsShort += cts.size();
				bkg.setTime(0.0);
				std::vector<double> bkgCts = bkg.getEventsToTime(dip9_arr[dipNum+1][1]-dip9_arr[dipNum][1]);
				sumShort += getSumTime(bkgCts);
				sumShort += (dip9_arr[dipNum][1])*bkgCts.size();
				ctsShort += bkgCts.size();
			}
			std::vector<double> cts = dip9Short[7].getEventsToTime(100);
			sumShort += getSumTime(cts);
			sumShort += (dip9_arr[7][1])*cts.size();
			ctsShort += cts.size();
			bkg.setTime(0.0);
			std::vector<double> bkgCts = bkg.getEventsToTime(100);
			sumShort += getSumTime(bkgCts);
			sumShort += (dip9_arr[7][1])*bkgCts.size();
			ctsShort += bkgCts.size();
			meanTimeShort.push_back(sumShort/ctsShort);
		}
		for(int j = 0; j < numRuns; j++) {
			double sumLong = 0.0;
			double ctsLong = 0;
			for(auto it = dip9Long.begin(); it < dip9Long.end(); it++) {
				it->setTime(0.0);
			}
			for(int dipNum = 0; dipNum < 7; dipNum++) {
				std::vector<double> cts = dip9Long[dipNum].getEventsToTime(dip9_arr[dipNum+1][1]-dip9_arr[dipNum][1]);
				sumLong += getSumTime(cts);
				sumLong += (dip9_arr[dipNum][1])*cts.size();
				ctsLong += cts.size();
				bkg.setTime(0.0);
				std::vector<double> bkgCts = bkg.getEventsToTime(dip9_arr[dipNum+1][1]-dip9_arr[dipNum][1]);
				sumLong += getSumTime(bkgCts);
				sumLong += (dip9_arr[dipNum][1])*bkgCts.size();
				ctsLong += bkgCts.size();
			}
			std::vector<double> cts = dip9Long[7].getEventsToTime(100);
			sumLong += getSumTime(cts);
			sumLong += (dip9_arr[7][1])*cts.size();
			ctsLong += cts.size();
			bkg.setTime(0.0);
			std::vector<double> bkgCts = bkg.getEventsToTime(100);
			sumLong += getSumTime(bkgCts);
			sumLong += (dip9_arr[7][1])*bkgCts.size();
			ctsLong += bkgCts.size();
			meanTimeLong.push_back(sumLong/ctsLong);
		}
		double avgDiff = 0.0;
		for(int j = 0; j < numRuns; j++) {
			avgDiff += (meanTimeLong[j]-meanTimeShort[j]);
		}
		avgDiff = avgDiff/numRuns;
		double stddev = 0.0;
		for(int j = 0; j < numRuns; j++) {
			stddev += pow((meanTimeLong[j]-meanTimeShort[j])-avgDiff, 2);
		}
		stddev = stddev/(numRuns - 1);
		stddev = sqrt(stddev);
		printf("Avg: %f stddev: %f\n", avgDiff, stddev);
	}
	
	return 0;
}