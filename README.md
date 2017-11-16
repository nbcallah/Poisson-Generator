# Nonhomogeneous Poisson Process Realizer

## Synopsis

Poisson Generator is a library to realize a nonhomogeneous poisson process. A nonhomogeneous poisson process is a poisson process where the rate is variable in time.

## Code Example

C code:

```
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

int main(int argc, char** argv) {
	std::vector<uint64_t> seed;
	seed.push_back(2736687128);
	seed.push_back(234302120);
	seed.push_back(3355772407);
	seed.push_back(657836083);
	std::seed_seq sseq(seed.begin(), seed.end());
	std::mt19937_64 gen;
	gen.seed(sseq);
	
	exponentialCountingProcess p(100.0, 10.0, &gen);
	
	std::vector<double> realization = p.getEventsToTime(20);
	
	printf("Generated %lu Events!\n", realization.size());
	
	return 0;
}

```

## Motivation

Generating a realization of a poisson process can be useful in generating synthetic datasets or in simulations. Classical poisson processes can be easily generated either by generating arrival times, or uniformly distributing poisson distributed numbers of events on an interval. When the rate is time variable, other techniques need to be employed such as thinning.

## Algorithm

This code is based off of "Generating Nonhomogeneous Poisson Processes" by Raghu Pasupathy. The method employed is the simple thinning algorithm as described. A set of arrival times is generated from a bounding rate. At each arrival time, the process is "thinned" by the ratio of the variable rate function / bounding rate. This works well for processes that have a rate comparable to the bounding rate, but fails to efficiently generate events on tails. However, it is a general method that needs only a bounding rate and a way to evaluate the time-dependent rate of interest.

## Installation

An example program (and makefile) is provided in the simulate_counting_with_background directory. The library is a single header file and a single source file; compilation for use in other programs is left as an exercise to the reader.

## API Reference (C++)

The library is supplied as a class that needs to be inherited and whose virtual methods need to be implemented.

### `poissonGen::poissonGen(double maxRate)`

Constructor that takes the bounding rate of the process as an argument.

### `double poissonGen::getNextEvent(double maxTime)`

Method that generates the next realization of the poisson process, or gives up if it does not find an event at t < maxTime (model time, not evalutation time).

### `std::vector<double> poissonGen::getEventsToTime(double length)`

Method that generates a vector of events from the poisson process. Goes to length in model time (model time, not evaluation time).

### `virtual double getUniformNumber()`

Virtual method that generates uniform random numbers between 0 and 1.

### `virtual double evaluateFunction(double t)`

Virtual method that evaluates the rate function at time t.

## Tests

Compile the source in the simulate_counting_with_background directory and run minimal. Output should be number of events generated in the first 20s of a process with a rate that goes ~exp(-t/10) with 100 events expected out to infinity.

The program sim contains a slightly more complicated scenario: A number of exponential counting processes (between 1-9) are concatenated. A uniform background is also added. The difference in mean arrival time between a process with high signal to background ratio and a process with a low signal to background. The average difference is calculated over many simulated runs.

## Contributors

As this is a small project, I will not actively maintain it. If you find this code, use it, and find a bug, feel free to let me know.

## Sources

I have based the algorithms on "Generating Nonhomogeneous Poisson Processes" by Raghu Pasupathy.

## Authors

Nathan Callahan

## License

This project is licensed under the MIT License - see the LICENSE.md file for details
