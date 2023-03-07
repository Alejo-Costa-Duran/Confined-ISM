#pragma once
#include <iostream>
using namespace std;
#ifndef Pvec_H
#define Pvec_H
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

class Pvector
{

public:
	double x, y;
	//Constructors
	Pvector() = default;
	Pvector(double x, double y)
	{
		this->x = x;
		this->y = y;
	}

	//Scalar functions
	void mulScalar(double s);
	void addScalar(double s);

	//Vector functions
	void addVector(Pvector v);
	void subVector(Pvector v);
	double dotProd(Pvector v);
	double crossProd(Pvector v);

	void normalize(double norma);
	double getNorm();
};

#endif
