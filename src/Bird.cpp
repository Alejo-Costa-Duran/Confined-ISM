#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include "Bird.h"
#include "Pvector.h"
#include "Flock.h"

Bird::Bird(int idx, Pvector position, Pvector velocity, Pvector acc, double inertia, double maxSpeed, double length)
{
	this->position = position;
	this->velocity = velocity;
	boundary(length,length);
	boxX = floor((position.x+length)*floor(2*length)/(2*length));
	boxY = floor((position.y+length)*floor(2*length)/(2*length));
	this->acc = acc;
	this->idx = idx;
	deltaV = Pvector(0, 0);
	randA = Pvector(0, 0);
	randV = Pvector(0, 0);
}

Pvector Bird::polarToScreen(double width, double height, double lengthX, double lengthY)
{
	Pvector screenPosition;
	screenPosition.x = (width / (2 * lengthX)) * (lengthX + position.x);
	screenPosition.y = (height / (2 * lengthY)) * (lengthY - position.y);
	return screenPosition;
}

void Bird::calcLambda(double c1, double c2, double dt, double maxSpeed)
{
	Pvector delV = Pvector(dt * c1 * acc.x + c2 * dt * dt * force.x + randV.x, dt * c1 * acc.y + c2 * dt * dt * force.y + randV.y);
	double b = 2 * (velocity.x * delV.x + delV.y * velocity.y);
	double a = maxSpeed * maxSpeed;
	double c = delV.x * delV.x + delV.y * delV.y - a;
	double w;
	if (b == 0) { w = sqrt(-c / a); }
	else
	{
		double sgnb = b > 0 ? 1.0 : -1.0;
		double q = -0.5 * (b + sgnb * sqrt(b * b - 4 * a * c));
		w = b > 0 ? c / q : q / a;
	}
	lambda = (w - 1) / (c2 * dt * dt);
}

void Bird::updatePosition(double dt)
{
	position.x += dt * velocity.x;
	position.y += dt * velocity.y;
}

void Bird::calcMu(double c2, double dt, double maxSpeed)
{
	mu = -(velocity.x * acc.x + velocity.y * acc.y) / (maxSpeed * maxSpeed * dt * c2);
}

void Bird::updateVelocity(double dt, double c1, double c2)
{
	velocity.x += c1 * dt * acc.x + dt * dt * c2 * force.x + dt * dt * c2 * lambda * velocity.x + randV.x;
	velocity.y += c1 * dt * acc.y + dt * dt * c2 * force.y + dt * dt * c2 * lambda * velocity.y + randV.y;
}

void Bird::calcForce(double inertia, double coupling, std::vector<Bird>& vecinos)
{
	double currentFx = 0;
	double currentFy = 0;
	for (auto& bird : vecinos)
	{
		currentFx += bird.velocity.x;
		currentFy += bird.velocity.y;
	}
	force.x = currentFx * coupling / inertia;
	force.y = currentFy * coupling / inertia;
}

void Bird::partialUpdate(double dt, double maxS, double c0, double c1, double c2)
{
	updatePosition(dt);
	calcLambda(c1, c2, dt, maxS);
	double tempx = (c1 - c2) * dt * (force.x + lambda * velocity.x);
	double tempy = (c1 - c2) * dt * (force.y + lambda * velocity.y);
	updateVelocity(dt, c1, c2);
	acc.x = c0 * acc.x + tempx;
	acc.y = c0 * acc.y + tempy;
}

double Bird::getSpin(double inertia, double maxSpeed)
{
	return (inertia / (maxSpeed * maxSpeed)) * (velocity.x * acc.y - velocity.y * acc.x);
}

void Bird::finalUpdate(double c2, double maxS, double dt, double inertia, double coupling, std::vector<Bird>& vecinos)
{
	calcForce(inertia, coupling, vecinos);
	acc.x += c2 * dt * force.x + randA.x;
	acc.y += c2 * dt * force.y + randA.y;
	calcMu(c2, dt, maxS);
	acc.x += mu * velocity.x;
	acc.y += mu * velocity.y;
}

void Bird::boundary(double lengthX, double lengthY)
{
	if (position.x > lengthX) { position.x -= 2 * lengthX; }
	if (position.y > lengthY) { position.y -= 2 * lengthY; }
	if (position.x < -lengthX) { position.x += 2 * lengthX; }
	if (position.y < -lengthX) { position.y += 2 * lengthY; }
}

double Bird::getDistance(Bird vecino, double length)
{
	double deltaX = abs(position.x - vecino.position.x);
	if (deltaX > length) { deltaX = deltaX - 2 * length; }
	double deltaY = abs(position.y - vecino.position.y);
	if (deltaY > length) { deltaY = deltaY - 2 * length; }
	return sqrt(deltaX * deltaX + deltaY * deltaY);
}
