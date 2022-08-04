#include "manager.h"
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <SFML/Graphics.hpp>
#include <fstream>

#define M_PI 3.14159926535897932385

void simRand(unsigned long int seed,int numSteps,int measure, int N, double boxSize,
	double inert, double coupling, double maxs, double timestep,
	double eta, double temp, std::string filename)
{
	ofstream file("../Data/Ferro/" + filename);
	Flock fl = Flock(seed, boxSize, timestep, maxs, inert, eta, temp, coupling, 1.0);
	for (int i = 0; i < N; i++)
	{
		double theta = gsl_rng_uniform(fl.m_mt) * M_PI * 2.0, posx = 2.0 * boxSize * gsl_rng_uniform(fl.m_mt)-boxSize, posy = 2.0 * boxSize * gsl_rng_uniform(fl.m_mt) - boxSize;
		Pvector pos(posx, posy);
		Pvector vel(maxs * cos(theta), maxs * sin(theta));
		Pvector acc(0, 0);
		fl.bandada.push_back(Bird(i, pos, vel, acc, inert, maxs,boxSize));
	}
	
	for (unsigned int idx = 0; idx < fl.bandada.size(); idx++)
	{
		std::vector<Bird> vec = fl.flocking(idx);
		fl.bandada[idx].calcForce(inert, coupling, vec);
	}
	
	int c = measure;
	if (file.is_open())
	{
		std::string header = "Spin,SigmaSpin,vecx,sigmavecx,vecy,sigmavecy,Energy,SigmaEnergy,"+std::to_string(seed);
		file << header + "\n";
		/*for (int therm = 0; therm < 100000; therm++)
		{
			fl.updateFlock(timestep);
			fl.boundary();
		}*/
		for (int step = 0; step < numSteps; step++)
		{
			fl.updateFlock(timestep);
			//fl.boundary();
			if (c == measure)
			{
				cout.precision(17);
				std::vector<double> medidas = fl.measurementsFerro();
				for (unsigned int l = 0; l < medidas.size() - 1; l++) { file << std::fixed << medidas[l]<< ","; }
				file << medidas[medidas.size() - 1]<< std::endl;
				c = 0;
				std::cout<<temp<<"\t"<< step << "\n";
			}
			c++;
		}
	}
}

void simCircles(int numRadi, double boxSize,double separation, int numAngles,double inert, double coupling, double maxs, double timestep, double eta, double temp, std::string filename)
{
	ofstream file("../Data/Ferro"+filename);
	Flock fl = Flock(0,boxSize,timestep, maxs, inert, eta, temp, coupling, 1.0);
	int l = 0;
	for (int i = 1; i < numRadi + 1; i++)
		for (int j = 0; j < numAngles * i; j++)
		{
			double theta = 2.0 * M_PI * j / (1.0*numAngles * i);
			double radius = separation * i + 2;
			Pvector pos(radius * cos(theta), radius * sin(theta));
			Pvector vel(-maxs * sin(theta), maxs * cos(theta));
			Pvector acc(-maxs * maxs * cos(theta) / radius, -maxs * maxs * sin(theta) / radius);
			fl.bandada.push_back(Bird(l, pos, vel, acc, inert, maxs,boxSize));
			++l;
		}
	for (unsigned int idx = 0; idx < fl.bandada.size(); idx++)
	{
		std::vector<Bird> vec = fl.flocking(idx);
		fl.bandada[idx].calcForce(inert, coupling, vec);
	}
	int c = 100;
	if (file.is_open())
	{
		char s[8] = "Spin   ";
		char ss[8] = "SigSpin";
		char e[8] = "Ener   ";
		char ee[8] = "SigEner";
		char ts[8] = "TangSpd";
		char sts[8] = "SigmSpd";
		char vec[80] = "Vec1,Vec2,Vec3,Vec4,Vec5,Vec6,Vec7,Vec8,Vec9,Vec10  ";

		file << s << "," << ss << "," << ts << "," << sts << "," << e << "," << ee << "," << vec << std::endl;
		for (int i = 0; i < 100000; i++)
		{
			fl.updateFlock(timestep);
			fl.boundary();
			if (c == 100)
			{
				cout.precision(17);
				std::vector<double> medidas = fl.measurements();
				for (unsigned int l = 0; l < medidas.size() - 1;l++) { file << std::fixed << medidas[l] << ","; }
				file << medidas[medidas.size()-1] << std::endl;
				c = 0;
				//std::cout << i << "\n";
			}
			c++;
		}
		file.close();
	}
}

int main()
{
	/*
	Manager man;
	man.Run(1.0);
*/
	
	for (int temp = 0; temp < 11; temp++)
	{
		double t = 1.0 * temp+0.01;
		int seed = rand();
		simRand(seed, 700000, 100, 1100, 9.586, 5.0, 1.0, 1.0,0.01 , 1.0, t, "aDataT"+std::to_string(t) + ".csv");
		//simCircles(10,15,0.8,20, 5.0, 1.0, 1.0,0.01 , 1.0, t, "DataasT"+std::to_string(t) + ".csv");
	}

}