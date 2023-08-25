#pragma once
#include <vector>
#include "./Bird.h"
#include <random>
#ifndef FLOCK_H_
#define FLOCK_H_
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <iostream>
#include <fstream>

class Flock
{
public:
	std::vector<Bird> bandada;
	double maxSpeed, inertia, friction, temperature, coupling;
	double numNeigh;
	double c0, c1, c2, exi;
	double siga, sigv, rho;
	double boxSize, cells;
	unsigned long int seed;
	double fil;
	Pvector field;
	gsl_rng *m_mt;

	Flock() = default;
	Flock(unsigned long int seed, double boxSize, double dt, double maxs, double iner, double fric, double T, double J, double N,double f);
    void writeRng(char *c);
    void saveFile(char *c,std::string envs, std::string birds,double dt);
	void loadFile(FILE *rng, std::string envs, std::string birds);
	void randVecs(Bird &b);
	void boundary();
	void circularBoundary(double,double);
	std::vector<Bird> flocking(int index);
	void updateFlock(double dt, double fieldStrength);
	double totalSpin();
	double totalEnergy();
	double totalTangSpeed();
	double avgRad();
	void setConstants(double dt);
	std::vector<double> measurements();
	std::vector<double> measurementsFerro();

	unsigned int test();
};

#endif
