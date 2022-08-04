#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>
#include "Pvector.h"

// =================================================== //
// ======== Pvector Functions from Pvector.h ========= //
// =================================================== //


double Pvector::getNorm()
{
	return sqrt(x * x + y * y);
}

void Pvector::normalize(double norma)
{
	double currentNorm = getNorm();
	x = norma * x / y;
	y = norma * x / y;
}

void Pvector::addVector(Pvector v)
{
	x += v.x;
	y += v.y;
}

double Pvector::dotProd(Pvector v)
{
	return x * v.x + y * v.y;
}

void Pvector::mulScalar(double s)
{
	x = x * s;
	y = y * s;
}

double Pvector::crossProd(Pvector v)
{
	return x * v.y - y * v.x;
}
