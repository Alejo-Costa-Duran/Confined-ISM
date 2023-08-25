#include "Simulation.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <Flock.h>
#include <bits/stdc++.h> 

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

void printProgress(double percentage, std::string run_id) {
    std::cout<< "\t Progress on file id: " << run_id;
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush(stdout);
}


Simulation::Simulation()
{
}

Simulation::Simulation(std::string sim_file)
{
    simulation_inputs_file = sim_file;
    readSimFile(simulation_inputs_file);
}

void Simulation::simMosquito(std::map<std::string,double> par, std::string folder, std::string run_id, std::string fData)
{
    std::ofstream filePos(folder + "/Pos" + run_id + ".csv");
    std::ofstream fileVel(folder + "/Vel" + run_id + ".csv");
    std::ofstream fileSpin(folder + "/Spin" + run_id + ".csv");
    std::string rngFile = folder + "/CheckPoints/Rng"+run_id;
    char *cc =  const_cast<char*>(rngFile.c_str());
    std::string enviroment = folder + "/CheckPoints/env" + run_id + ".csv";
    std::string checkpointFile = folder + "/CheckPoints/checkFile" + run_id + ".csv";
    std::ofstream fileData;
    fileData.open(folder + "/" + fData + ".csv",ios::app);
    for(auto key : parameters)
    {
        fileData << key.second << ",";
    }
    fileData << "\n";
    fileData.close();
    int seed = par["Seed"];
    double boxSize = par["BoxSize"];
    double timestep = par["TimeStep"];
    double inert = par["Inertia"], fric = par["Friction"], temperature = par["Temperature"], interaccion = par["Interaction"];
    int numBirds = par["NumBirds"];
    double field = par["Field"];
    int sampling = par["Sampling"];
    int totSteps = par["TotalSteps"];
    int checkpoint = par["Checkpoint"];
    
    Flock fl(seed,boxSize,timestep,1.0,inert,fric,temperature,interaccion,1.0,0.0);
    
    std::string headerPos = "PosX0,PosY0";
    std::string headerVel = "VelX0,VelY0";
    std::string headerSpin = "Spin0";
    for(int idx=1; idx<numBirds; idx++)
    {
        headerPos = headerPos + ",PosX" + std::to_string(idx) + ",PosY" + std::to_string(idx);
        headerVel = headerVel + ",VelX" + std::to_string(idx) + ",VelY" + std::to_string(idx);
        headerSpin = headerSpin + ",Spin" + std::to_string(idx);
    
    }
    filePos<< headerPos <<"\n";
    fileVel << headerVel <<"\n";
    fileSpin << headerSpin <<"\n";
    for(int idx = 0; idx<numBirds; idx++)
    {
        Pvector Pos(1,0);
        Pvector Vel(0,1);
        Pvector acc(0,0);
        fl.bandada.push_back(Bird(false,idx,Pos,Vel,acc,inert,1.0,boxSize));
    }
    for (unsigned int idx = 0; idx < fl.bandada.size(); idx++)
	{
		std::vector<Bird> vec = fl.flocking(idx);
		fl.bandada[idx].calcForce(inert, interaccion, vec,field);
	}
    int samp = sampling;
    int checkFlag = checkpoint;
    double progress = 0;
    for(int step = 0; step<totSteps;step++)
    {
        if(samp == sampling)
        {
            for(int idx=0; idx<numBirds-1; idx++)
            {
            cout.precision(17);
            filePos << std::setprecision(10) << fl.bandada[idx].position.x << "," <<fl.bandada[idx].position.y <<",";
            fileVel << std::setprecision(10) << fl.bandada[idx].velocity.x << "," <<fl.bandada[idx].velocity.y <<",";
            fileSpin << std::setprecision(10) << fl.bandada[idx].getSpin(inert,1.0) <<",";
            samp = 0;
            }
            cout.precision(17);
            filePos << std::setprecision(10) << fl.bandada[numBirds-1].position.x << "," <<fl.bandada[numBirds-1].position.y <<"\n";
            fileVel << std::setprecision(10) << fl.bandada[numBirds-1].velocity.x << "," <<fl.bandada[numBirds-1].velocity.y <<"\n";
            fileSpin << std::setprecision(10) << fl.bandada[numBirds-1].getSpin(inert,1.0) <<"\n";
            samp = 0;
            printProgress(1.0*progress/totSteps, run_id);
            
        }
        
        if(checkFlag == checkpoint)
        {
            fl.saveFile(cc,enviroment,checkpointFile,timestep);
            checkFlag = 0;
        }
        fl.updateFlock(timestep,field);
        samp++;
        checkFlag++;
        progress += 1.0;
    }

    filePos.close();
    fileVel.close();
    fileSpin.close();
}

void Simulation::simAnnulus(std::map<std::string,double> par, std::string folder, std::string run_id, std::string fData)
{
    std::ofstream filePos(folder + "/Pos" + run_id + ".csv");
    std::ofstream fileVel(folder + "/Vel" + run_id + ".csv");
    std::ofstream fileSpin(folder + "/Spin" + run_id + ".csv");
    std::ofstream fileAng(folder + "/Ang" + run_id + ".csv");
    std::string rngFile = folder + "/CheckPoints/Rng"+run_id;
    char *cc =  const_cast<char*>(rngFile.c_str());
    std::string enviroment = folder + "/CheckPoints/env" + run_id + ".csv";
    std::string checkpointFile = folder + "/CheckPoints/checkFile" + run_id + ".csv";

    std::ofstream fileData;
    fileData.open(folder + "/" + fData + ".csv",ios::app);
    for(auto key : parameters)
    {
        cout.precision(17);
        fileData << std::setprecision(17)<< key.first << ",";
    }
    fileData << "Run Id" <<"\n";
    for(auto key : parameters)
    {
        fileData << key.second << ",";
    }
    fileData << run_id <<"\n";
    fileData.close();
    int seed = par["Seed"];
    double boxSize = par["BoxSize"];
    double timestep = par["TimeStep"];
    double inert = par["Inertia"], fric = par["Friction"], temperature = par["Temperature"], interaccion = par["Interaction"];
    int numBirds = par["NumBirds"];
    double field = par["Field"];
    int sampling = par["Sampling"];
    int totSteps = par["TotalSteps"];
    int checkpoint = par["Checkpoint"];

    Flock fl(seed,boxSize,timestep,1.0,inert,fric,temperature,interaccion,1.0,0.0);


    std::string headerPos = "PosX,PosY";
    std::string headerVel = "VelX,VelY";
    std::string headerSpin = "Spin";
    filePos<< headerPos <<"\n";
    fileVel << headerVel <<"\n";
    fileSpin << headerSpin <<"\n";
    fileAng << "Ang" <<"\n";
    for(int idx = 0; idx<numBirds; idx++)
    {
        double r0=7.5+4*(-0.5+gsl_rng_uniform(fl.m_mt)), theta0=2*M_PI*gsl_rng_uniform(fl.m_mt);
        double dir0=2*M_PI*gsl_rng_uniform(fl.m_mt);
        Pvector Pos(r0*cos(theta0),r0*sin(theta0));
        Pvector Vel(-sin(theta0),cos(theta0));
        Pvector acc(-cos(theta0)/r0,-sin(theta0)/r0);
        fl.bandada.push_back(Bird(false,idx,Pos,Vel,acc,inert,1.0,boxSize));
    }
    for (unsigned int idx = 0; idx < fl.bandada.size(); idx++)
	{
		std::vector<Bird> vec = fl.flocking(idx);
		fl.bandada[idx].calcForce(inert, interaccion, vec,field);
	}


    int samp = sampling;
    int checkFlag = checkpoint;
    double progress = 0;
    for(int step = 0; step<totSteps;step++)
    {
        if(samp == sampling)
        {
            for(int idx=0; idx<numBirds-1; idx++)
            {
            cout.precision(17);
            filePos << std::setprecision(10) << fl.bandada[idx].position.x << "," <<fl.bandada[idx].position.y <<",";
            fileVel << std::setprecision(10) << fl.bandada[idx].velocity.x << "," <<fl.bandada[idx].velocity.y <<",";
            fileSpin << std::setprecision(10) << fl.bandada[idx].getSpin(inert,1.0) <<",";
            double L = fl.bandada[idx].position.x*fl.bandada[idx].velocity.y-fl.bandada[idx].position.y*fl.bandada[idx].velocity.x;
            fileAng<< std::setprecision(10) << L<<",";
            samp = 0;
            }
            cout.precision(17);
            filePos << std::setprecision(10) << fl.bandada[numBirds-1].position.x << "," <<fl.bandada[numBirds-1].position.y <<"\n";
            fileVel << std::setprecision(10) << fl.bandada[numBirds-1].velocity.x << "," <<fl.bandada[numBirds-1].velocity.y <<"\n";
            fileSpin << std::setprecision(10) << fl.bandada[numBirds-1].getSpin(inert,1.0) <<"\n";
            double L = fl.bandada[numBirds-1].position.x*fl.bandada[numBirds-1].velocity.y-fl.bandada[numBirds-1].position.y*fl.bandada[numBirds-1].velocity.x;
            fileAng<< std::setprecision(10) << L<<"\n";
            samp = 0;
            printProgress(progress/totSteps, run_id);
        }
        if(checkFlag == checkpoint)
        {
            fl.saveFile(cc,enviroment,checkpointFile,timestep);
            checkFlag = 0;
        }
        fl.updateFlock(timestep,field);
        samp++;
        checkFlag++;
        progress += 1.0;
    }

    filePos.close();
    fileVel.close();
    fileSpin.close();
    fileAng.close();
}

void Simulation::print_usage()
{
    for(auto &key : parameters)
    {
        std::cout << key.first << ":" << key.second <<"\n";
    }
}
void Simulation::set_par(std::string par, double value)
{
    parameters[par] = value;
}
void Simulation::run(std::string type,std::string folder, std::string run_id,std::string infoFile)
{

    if(type=="Field")
    {
        simAnnulus(parameters, folder, run_id,infoFile);
    }
    else if(type=="Mosquito")
    {
        std::cout<<"ok linea run 248"<<"\n";
        simMosquito(parameters, folder, run_id,infoFile);   
    }
    else if(type=="Caja")
    {
        std::cout << "Not yet \n";
        exit(1);
    }
    else
    {
        std::cout << "Simulación de tipo '" + type + "' no implementada \n";
        exit(1);
    }
}

void Simulation::readSimFile(std::string file)
{
    std::ifstream in;
    in.open(file);

    if(!in.is_open())
    {
        std::cout << "Error loading the file" << std::endl;
        in.close();
        return;
    }

    std::string line;
    std::vector<std::string> row;
    std::string buffer;
    while(std::getline(in,line))
    {
        double num;
        std::stringstream st(line);
        st>>buffer;
        st>>num;
        parameters[buffer] = num;
    }
}