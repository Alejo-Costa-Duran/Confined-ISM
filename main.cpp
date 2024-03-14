#include "manager.h"
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <SFML/Graphics.hpp>
#include <fstream>
#include <iomanip>
#include <Simulation.h>
#include <chrono>
#include <algorithm>

void simCaja(int sampling, int seed, double numBirds, double boxSize, double temperature, std::string filename)
{
    ofstream filePos("../../Data/CampoExterno/Caja/Pos"+filename);
    ofstream fileVel("../../Data/CampoExterno/Caja/Vel"+filename);
    Flock fl(seed,boxSize,0.001,1.0,1.0,1.0,temperature,1.0,1.0,0.0);

    for(int idx = 0; idx<numBirds; idx++)
    {
        double posx=2*boxSize*(-0.5+gsl_rng_uniform(fl.m_mt)), posy = posx=2*boxSize*(-0.5+gsl_rng_uniform(fl.m_mt));
        double dir0=2*M_PI*gsl_rng_uniform(fl.m_mt);
        Pvector Pos(posx,posy);
        Pvector Vel(-sin(dir0),cos(dir0));
        Pvector acc(0,0);
        fl.bandada.push_back(Bird(false,idx,Pos,Vel,acc,1.0,1.0,boxSize));
    }
    std::cout<<filename<<"\n";
    for (unsigned int idx = 0; idx < fl.bandada.size(); idx++)
	{
		std::vector<Bird> vec = fl.flocking(idx);
		fl.bandada[idx].calcForce(1,1, vec,0);
	}
    int samp = 10*sampling;
    for(int step = 0; step<200001;step++)
    {
        if(samp == 10*sampling)
        {
            for(int idx=0; idx<numBirds-1; idx++)
            {
            cout.precision(17);
            filePos << std::setprecision(10) << fl.bandada[idx].position.x << "," <<fl.bandada[idx].position.y <<",";
            fileVel << std::setprecision(10) << fl.bandada[idx].velocity.x << "," <<fl.bandada[idx].velocity.y <<",";
            samp = 0;
            }
            cout.precision(17);
            filePos << std::setprecision(10) << fl.bandada[numBirds-1].position.x << "," <<fl.bandada[numBirds-1].position.y <<"\n";
            fileVel << std::setprecision(10) << fl.bandada[numBirds-1].velocity.x << "," <<fl.bandada[numBirds-1].velocity.y <<"\n";
            samp = 0;
        }
        fl.updateFlock(0.001,0);
        samp++;
    }
    for(int step = 0; step<100000;step++)
    {
        
        if(samp == sampling)
        {
            for(int idx=0; idx<numBirds-1; idx++)
            {
            cout.precision(17);
            filePos << std::setprecision(10) << fl.bandada[idx].position.x << "," <<fl.bandada[idx].position.y <<",";
            fileVel << std::setprecision(10) << fl.bandada[idx].velocity.x << "," <<fl.bandada[idx].velocity.y <<",";
            samp = 0;
            }
            cout.precision(17);
            filePos << std::setprecision(10) << fl.bandada[numBirds-1].position.x << "," <<fl.bandada[numBirds-1].position.y <<"\n";
            fileVel << std::setprecision(10) << fl.bandada[numBirds-1].velocity.x << "," <<fl.bandada[numBirds-1].velocity.y <<"\n";
            samp = 0;
        }
        fl.updateFlock(0.001,0);
        samp++;
    }
    filePos.close();
    fileVel.close();
}

int main(int arg_count, char* args[])
{
   if(arg_count>1)
   {
    std::string folder = args[1];
    std::string run_id = args[2];
    std::string data_file = args[3];
    std::cout << "Los archivos se van a guardar en la carpeta: " + folder + "\n";
    std::cout << "Con el sufijo: " + run_id<<"\n";
    std::cout << "Los parámetros de la simulación se tomaron del archivo" + data_file<<"\n";
    Simulation s(data_file);
    int counter = 200;
    std::vector<double> inertias ={0.0003,0.005,0.1,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,5,10,15,20,30,50,100};
    std::vector<double> fields = {0.0,0.001,0.005,0.01,0.015,0.02,0.05,0.08,0.1,0.25,0.5,0.75,1.0,1.05,1.1,1.25,1.5,1.75,2.0};
    for(auto chi: inertias)
    {
        for(auto f: fields)
        {        
            std::cout << "Amount of runs: " << 1 <<"\n";
            int seed = rand();
            s.set_par("Inertia",chi);
            s.set_par("Field",f);
            s.set_par("Seed",seed);
            run_id = "Field"+std::to_string(counter);
            counter +=1;
            auto start = std::chrono::high_resolution_clock::now();
            s.simMosquito(s.parameters, folder, run_id,"/Info");
            auto stop = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
            std::cout << "\n Finished simulation on file id: " + run_id;
            std::cout << "\n It took: " << duration.count() << " seconds \n";
        }
    }
   }
   else
   {
    std::cout << "Error: Es necesario proveer el nombre de la carpeta y el sufijo de los archivos \n";    return 1;
   }
}
