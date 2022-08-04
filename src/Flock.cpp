#define _USE_MATH_DEFINES
#include "Bird.h"
#include "Flock.h"
#include <algorithm>
#include <math.h>
#include <random>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>



// =============================================== //
// ======== Flock Functions from Flock.h ========= //
// =============================================== //

Flock::Flock(unsigned long int seed,double boxSize,double dt, double maxSp, double iner, double fric, double T, double J, double N)
{
    inertia = iner;
    this->seed = seed;
    friction = fric;
    temperature = T;
    coupling = J;
    numNeigh = N;
    cells = boxSize*2/floor(2*boxSize);
    maxSpeed = maxSp;
    setConstants(dt);
    this->m_mt = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(m_mt, seed);
    this->boxSize = boxSize;
}


void Flock::randVecs(Bird& b)
{
    gsl_ran_bivariate_gaussian(m_mt, siga, sigv, rho, &b.randA.x, &b.randV.x);
    gsl_ran_bivariate_gaussian(m_mt, siga, sigv, rho, &b.randA.y, &b.randV.y);
}

void Flock::setConstants(double dt)
{
    double etasv0 = friction / (maxSpeed * maxSpeed);
    double xi = etasv0 / (inertia);
    double xidt = xi * dt;
    this->exi = exp(-xidt);
    if (xidt < 1e-3)
    {
        c0 = 1 - xidt + xidt * xidt / 2 - xidt * xidt * xidt / 6;
        c1 = 1 - xidt / 2 + xidt * xidt / 6 - xidt * xidt * xidt / 24;
        c2 = 0.5 - xidt / 6 + xidt * xidt / 24;
        rho = sqrt(3.) * (0.5 - xidt / 16. - (17. / 1280.) * xidt * xidt + (17. / 6144) * xidt * xidt * xidt);
    }
    else
    {
        c0 = exi;
        c1 = (1 - c0) / xidt;
        c2 = (1 - c1) / xidt;
        rho = (1 - exi) * (1 - exi) / sqrt((1 - exi * exi) * (3 * xidt - 3 + 4 * exi - exi * exi));
    }
    siga = (temperature / inertia) * (1 - exi * exi);
    siga = sqrt(siga);
    if (etasv0 < 1e-3)
    {
        sigv = temperature * dt * dt * dt * etasv0 * (2. / 3. - 0.5 * xidt) / (inertia * inertia);
    }
    else
    {
        sigv = (temperature / etasv0) * (2 * dt - (3 - 4 * exi + exi * exi) / xi);
    }
    sigv = sqrt(sigv);
}

void Flock::boundary()
{
    for (unsigned int l = 0; l < bandada.size(); l++)
    {
        bandada[l].boundary(boxSize,boxSize);
    }
}

double Flock::totalTangSpeed()
{
    double tang = 0;
    for (unsigned int index = 0; index < bandada.size(); index++)
    {
        double angle = atan2(bandada[index].position.y, bandada[index].position.x);
        tang += -bandada[index].velocity.x * sin(angle) + bandada[index].velocity.y * cos(angle);
    }
    return tang;
}

vector<Bird> Flock::flocking(int index)
{
    Bird pajaro = bandada[index];
    /*
    std::vector<Bird> tmp;
    copy(bandada.begin(), bandada.end(), back_inserter(tmp));
    sort(tmp.begin(), tmp.end(), [&](const Bird a, const Bird b) {return bandada[index].getDistance(a) < bandada[index].getDistance(b); });
    std::vector<Bird> vecinos;
    copy(tmp.begin() + 1, tmp.begin() + 1 + numNeigh, back_inserter(vecinos));

    */
    std::vector<Bird> vecinos;
    int box = 2*boxSize-1;
    for (auto& b : bandada)
    {

        if(((abs(b.boxX - pajaro.boxX)%box)<3) && ((abs(b.boxY - pajaro.boxY)%box)<3))
        { 
            if (pajaro.getDistance(b, boxSize) < numNeigh && pajaro.getDistance(b, boxSize) != 0)
            {
                vecinos.push_back(b);
            }
        }
    }

    return vecinos;
}

double Flock::avgRad()
{
    double rad = 0;
    for (Bird b : bandada)
    {
        rad += sqrt(pow(b.position.x, 2) + pow(b.position.y, 2));
    }
    return rad / bandada.size();
}

void Flock::updateFlock(double dt)
{
    for (Bird& b : bandada)
    {
        //Calc  noise
        randVecs(b);

        //Update position
        b.position.x += b.velocity.x * dt;
        b.position.y += b.velocity.y * dt;

        //
        b.boundary(boxSize,boxSize);

        //Update box
	    b.boxX = floor((b.position.x+boxSize)/cells);
	    b.boxY = floor((b.position.y+boxSize)/cells);

        //Velocity without constraint
        b.deltaV.x = c1 * dt * b.acc.x + c2 * dt * dt * b.force.x + b.randV.x;
        b.deltaV.y = c1 * dt * b.acc.y + c2 * dt * dt * b.force.y + b.randV.y;

        //Lagrange multiplier
        double b0 = 2 * (b.velocity.x * b.deltaV.x + b.deltaV.y * b.velocity.y);
        double a = maxSpeed * maxSpeed;
        double c = b.deltaV.x * b.deltaV.x + b.deltaV.y * b.deltaV.y - a;
        double w;
        if (b0 == 0) { w = sqrt(-c / a); }
        else
        {
            double sgnb = b0 > 0 ? 1.0 : -1.0;
            double q = -0.5 * (b0 + sgnb * sqrt(b0 * b0 - 4 * a * c));
            w = b0 > 0 ? c / q : q / a;
        }
        b.lambda = (w - 1) / (c2 * dt * dt);

        //Partial update of acc
        b.acc.x = c0 * b.acc.x + (c1 - c2) * dt * (b.force.x + b.lambda * b.velocity.x) + b.randA.x;
        b.acc.y = c0 * b.acc.y + (c1 - c2) * dt * (b.force.y + b.lambda * b.velocity.y) + b.randA.y;

        //Update v
        b.velocity.x = (1+c2 * dt * dt * b.lambda) * b.velocity.x + b.deltaV.x;
        b.velocity.y = (1+c2 * dt * dt * b.lambda) * b.velocity.y + b.deltaV.y;

        //Update forces
        std::vector<Bird> vec = flocking(b.idx);
        b.calcForce(inertia, coupling, vec);

        //Update a
        b.acc.x += c2 * dt * b.force.x;
        b.acc.y += c2 * dt * b.force.y;

        //Calc mu
        b.mu = -b.velocity.dotProd(b.acc) / (maxSpeed * maxSpeed);

        //Final update a
        b.acc.x += b.velocity.x * b.mu;
        b.acc.y += b.velocity.y * b.mu;
    }
}

double Flock::totalSpin()
{
    double totSpin = 0;
    for (auto& b : bandada)
    {
        totSpin += b.getSpin(inertia, maxSpeed);
    }
    return totSpin;
}

std::vector<double> Flock::measurementsFerro()
{
    std::vector<double> meanValues;
    int N = bandada.size();
    double spin = 0, energy = 0, avgVecx = 0, avgVecy = 0, squareSpin = 0;
    double squareAvgVecX = 0, squareAvgVecY = 0, squareEnergy = 0;
    for (int idx = 0; idx < N; idx++)
    {
        Bird pajaro = bandada[idx];
        double e = 0;
        double s = pajaro.getSpin(inertia, maxSpeed);
        spin += s;
        squareSpin += s * s;
        avgVecx += pajaro.velocity.x;
        avgVecy += pajaro.velocity.y;
        squareAvgVecX += pajaro.velocity.x * pajaro.velocity.x;
        squareAvgVecY += pajaro.velocity.y * pajaro.velocity.y;
        e += s * s / (2 * inertia);
        std::vector<Bird> vec = flocking(idx);
        for (Bird b : vec) { e  += -(coupling / (2 * maxSpeed * maxSpeed)) * (pajaro.velocity.dotProd(b.velocity)); }
        energy += e;
        squareEnergy += e * e;
    }
    std::vector<double> measures = { spin/N, squareSpin/N, avgVecx/N,squareAvgVecX/N, avgVecy/N, squareAvgVecY/N, energy/N, squareEnergy/N };
    return measures;
}

std::vector<double> Flock::measurements()
{
    std::vector<double> meanValues;
    std::vector<double> spins;
    std::vector<double> tangSpeeds;
    std::vector<double> energies;
    std::vector<double> vecinos;
    Pvector centerMass(0,0);
    int N = bandada.size();
    for (Bird b : bandada) { centerMass.x += b.position.x; centerMass.y += b.position.y; }
    centerMass.x = centerMass.x / N;
    centerMass.y = centerMass.y / N;
    double spin = 0, energy = 0, tangSpeed = 0, avgVec = 0;
    for (int idx=0; idx<N;idx++)
    {
        Bird pajaro = bandada[idx];
        spins.push_back(pajaro.getSpin(inertia, maxSpeed));
        spin += spins[idx];
        Pvector dist(centerMass.x - pajaro.position.x, centerMass.y - pajaro.position.y);
        tangSpeeds.push_back(pajaro.velocity.crossProd(dist) / dist.getNorm());
        tangSpeed += tangSpeeds[idx];
        
        energies.push_back(spins[idx] * spins[idx] / (2 * inertia));
        std::vector<Bird> vec = flocking(idx);
        vecinos.push_back(vec.size());
        avgVec += vecinos[idx];
        for (Bird b : vec) { energies[idx] += -(coupling / (2 * maxSpeed * maxSpeed)) * (pajaro.velocity.dotProd(b.velocity)); }
        energy += energies[idx];
    }
    std::vector<double> vecs;
    for (int i = 1; i < 11; i++)
    {
        vecs.push_back(0.0);
        for (int j = 10*(i-1)*(i); j < 10 * (i)*(i+1); j++)
        {
            vecs[i - 1] += vecinos[j]/(20*(i));
        }
    }
    spin = spin / N;
    tangSpeed = tangSpeed / N;
    energy = energy / N;
    avgVec = avgVec / N;

    double sigmaEnergy = 0 , sigmaSpin = 0, sigmaTang = 0, sigmaVec=0;
    for (int idx = 0; idx < N; idx++)
    {
        sigmaSpin += (spins[idx] - spin) * (spins[idx] - spin);
        sigmaTang += (tangSpeeds[idx] - tangSpeed) * (tangSpeeds[idx] - tangSpeed);
        sigmaEnergy += (energies[idx] - energy) * (energies[idx] - energy);
        sigmaVec += (vecinos[idx] - avgVec) * (vecinos[idx] - avgVec);
    }
    sigmaSpin = sqrt(sigmaSpin / N);
    sigmaTang = sqrt(sigmaTang / N);
    sigmaEnergy = sqrt(sigmaEnergy / N);
    sigmaVec = sqrt(sigmaVec / N);
    std::vector<double> measures = { spin, sigmaSpin, tangSpeed, sigmaTang, energy, sigmaEnergy};
    for (unsigned int idx = 0; idx < vecs.size(); idx++) { measures.push_back(vecs[idx]); }
    return measures;
}

double Flock::totalEnergy()
{
    double energy = 0;
    for (unsigned int index = 0; index < bandada.size(); index++)
    {
        energy += bandada[index].getSpin(inertia, maxSpeed) * bandada[index].getSpin(inertia, maxSpeed) / (2 * inertia);
        std::vector<Bird> vec = flocking(index);
        for (Bird b : vec)
        {
            energy += -coupling / (2 * maxSpeed * maxSpeed) * (bandada[index].velocity.x * b.velocity.x + bandada[index].velocity.y * b.velocity.y);
        }
    }
    return energy;
}
