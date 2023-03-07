#pragma once
#include "Pvector.h"
#include <vector>

#ifndef BOID_H_
#define BOID_H_

class Bird
{
public:
	// Attributes
	Pvector position, velocity, acc, force, deltaV;
	int idx, boxX, boxY;
	double lambda, mu;
	Pvector randA, randV;
	bool fixed;

	Bird() = default;
	Bird(bool fix, int idx, Pvector position, Pvector velocity, Pvector acc, double inertia, double maxSpeed, double length);

	void calcLambda(double c1, double c2, double dt, double maxSpeed);
	Pvector polarToScreen(double width, double height, double lengthX, double lengthY);
	void partialUpdate(double dt, double maxSp, double c0, double c1, double c2);
	void boundary(double lengthX, double lengthY);
	double getDistance(Bird vecino, double length);
	void updateVelocity(double dt, double c1, double c2);
	void updatePosition(double dt);
	void calcForce(double maxSpeed, double coupling, std::vector<Bird> &vecinos);
	void calcMu(double c2, double dt, double maxSpeed);
	void finalUpdate(double c2, double maxs, double dt, double inertia, double coupling, std::vector<Bird> &vecinos);
	double getSpin(double inertia, double maxSpeed);
};

#endif
